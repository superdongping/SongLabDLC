"""Single-job, preemptible idle processing queue, persisted in each session."""
import threading
import time
from pathlib import Path
from media_quality import finalize_video, ProcessingCancelled
from video_output import mp4_destination

class BackgroundSupport:
    def init_background(self):
        self.processing_cancel=threading.Event();self.processing_stop=threading.Event()
        self.processing_thread=None;self.processing_id=None;self.processing_force=False
        self.last_activity=time.monotonic();self.background_error=''
        self.scheduler=threading.Thread(target=self.background_loop,daemon=True);self.scheduler.start()

    def touch(self):self.last_activity=time.monotonic()

    def preempt_background(self):
        self.touch();self.processing_force=False;self.processing_cancel.set()
        t=self.processing_thread
        if t and t.is_alive():
            t.join(5)
            if t.is_alive():raise ValueError('Background process is stopping. Try again in a moment.')

    def queue_action(self,data):
        with self.command_lock:
            if not self.project_file:raise ValueError('Open a project first.')
            action=data.get('action')
            config=self.db.setdefault('processing',{'idle_seconds':60,'paused':False})
            if action=='pause':
                config['paused']=True;self.preempt_background()
            elif action=='resume':config['paused']=False;self.touch()
            elif action=='process':
                if self.active:raise ValueError('Close preview or calibration before processing videos.')
                config['paused']=False;self.processing_force=True
            elif action=='retry':
                self.preempt_background()
                for s in self.db['sessions']:
                    q=s['phases'].get('recording',{})
                    if q.get('mp4_status')=='failed':q.update(mp4_status='queued',qc='PENDING')
            elif action=='delay':
                seconds=float(data['seconds'])
                if not 5<=seconds<=3600:raise ValueError('Idle delay must be 5-3600 seconds.')
                config['idle_seconds']=seconds
            else:raise ValueError('Unknown processing action.')
            with self.lock:self.persist()

    def queue_state(self):
        jobs=[{'id':s['id'],'status':s['phases'].get('recording',{}).get('mp4_status','not queued')} for s in self.db['sessions'] if not s.get('deleted_at')]
        return dict(jobs=jobs,active=self.processing_id,error=self.background_error,
                    config=self.db.get('processing',{'idle_seconds':60,'paused':False}))

    def background_loop(self):
        while not self.processing_stop.wait(.25):
            if not self.command_lock.acquire(blocking=False):continue
            try:
                if self.active or not self.project_file or self.processing_id:continue
                config=self.db.get('processing',{})
                if config.get('paused'):continue
                if not self.processing_force and time.monotonic()-self.last_activity<config.get('idle_seconds',60):continue
                if self.processing_thread and self.processing_thread.is_alive():continue
                s=next((s for s in self.db['sessions'] if not s.get('deleted_at') and s['status'] in ('completed','interrupted')
                        and s['phases'].get('recording',{}).get('mp4_status')=='queued'),None)
                if not s:continue
                self.processing_cancel.clear();self.processing_id=s['id']
                with self.lock:
                    s['phases']['recording']['mp4_status']='processing';self.persist()
                self.processing_thread=threading.Thread(target=self.process_one,args=(s,),daemon=True);self.processing_thread.start()
            except Exception as exc:
                self.background_error=str(exc);self.processing_id=None
            finally:self.command_lock.release()

    def process_one(self,s):
        try:
            q=s['phases']['recording'];source=Path(s['folder'])/q['file']
            if not source.is_file():raise ValueError('Original MKV is missing. Restore the project data folder, then retry.')
            destination = mp4_destination(self.project_file.parent,q['mp4_relative']) if q.get('mp4_relative') else None
            result=finalize_video(self.ffmpeg_path(),source,s['capture']['fps'],s['recording_seconds'],self.process_flags,self.processing_cancel,destination=destination)
            with self.lock:
                q.update(result)
                if s['status']!='completed':q['qc']='REVIEW'
                q['processed_at']=__import__('datetime').datetime.now().astimezone().isoformat()
                self.save_session(s);self.persist()
        except ProcessingCancelled:
            with self.lock:
                s['phases']['recording'].update(mp4_status='queued',qc='PENDING');self.persist()
        except Exception as exc:
            with self.lock:
                s['phases']['recording'].update(mp4_status='failed',qc='REVIEW',mp4_error=str(exc));self.persist()
        finally:self.processing_id=None

    def close_services(self):
        self.processing_stop.set();self.preempt_background()
        if self.scheduler:self.scheduler.join(3)
        if self.project_handle:self.project_handle.close();self.project_handle=None
