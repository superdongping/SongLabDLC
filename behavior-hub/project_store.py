"""Portable project manifests. Recording files stay inside the selected project root."""
import copy
import datetime as dt
import json
import os
from pathlib import Path
import shutil
import uuid


def atomic(path, data):
    path=Path(path); tmp=path.with_name(path.name+'.writing')
    with tmp.open('w',encoding='utf-8') as f:
        json.dump(data,f,indent=2,ensure_ascii=False,allow_nan=False);f.flush();os.fsync(f.fileno())
    os.replace(tmp,path)


def inside(root, relative):
    path=Path(relative)
    if path.is_absolute(): raise ValueError('Project paths must be relative.')
    result=(root/path).resolve()
    if not result.is_relative_to(root.resolve()): raise ValueError('Project path escapes its folder.')
    return result


class ProjectSupport:
    def init_projects(self):
        self.project_file=None;self.project_handle=None
        self.legacy_db=copy.deepcopy(self.db)
        self.preferences_path=self.data_dir/'projects.json'
        self.recent_projects=[]
        if self.preferences_path.exists():
            self.recent_projects=json.loads(self.preferences_path.read_text(encoding='utf-8')).get('recent',[])

    def serialize_project(self):
        data=copy.deepcopy(self.db)
        data['output']='data'
        for s in data['sessions']:
            s['folder']=Path(s['folder']).resolve().relative_to(self.project_file.parent).as_posix()
        return data

    def project_persist(self):
        self.db['project']['saved_at']=dt.datetime.now().astimezone().isoformat()
        atomic(self.project_file,self.serialize_project())

    def remember_project(self):
        item={'path':str(self.project_file),'name':self.db['project']['name']}
        self.recent_projects=[item]+[x for x in self.recent_projects if x['path']!=item['path']][:9]
        atomic(self.preferences_path,{'recent':self.recent_projects})

    def project_state(self):
        return dict(project=self.db.get('project'),project_file=str(self.project_file) if self.project_file else None,
                    recent_projects=self.recent_projects,legacy_count=len(self.legacy_db['sessions']))

    def project_action(self,data):
        with self.command_lock:
            if self.active and self.active!='preview': raise ValueError('Finish the current recording before changing projects.')
            kind=data.get('action')
            if kind=='save':
                if not self.project_file: raise ValueError('Create or open a project first.')
                with self.lock:
                    if 'settings' in data:self.db['settings']=copy.deepcopy(data['settings'])
                    self.persist()
                return str(self.project_file)
            self.stop_preview();self.preempt_background()
            if kind=='import_legacy':
                if not self.project_file: raise ValueError('Open a project first.')
                imported=[]
                for old in self.legacy_db['sessions']:
                    if old.get('deleted_at'):continue
                    if any(s.get('legacy_id')==old['id'] for s in self.db['sessions']):continue
                    source=Path(old['folder'])
                    if not source.is_dir():continue
                    target=inside(self.project_file.parent,'data/IMPORTED/'+uuid.uuid4().hex)
                    shutil.copytree(source,target)
                    s=copy.deepcopy(old);s['folder']=str(target);s['legacy_id']=old['id']
                    if any(x['id']==s['id'] for x in self.db['sessions']):s['id']+='_'+uuid.uuid4().hex[:6]
                    s['project_id']=self.db['project']['id'];self.db['sessions'].append(s);self.save_session(s);imported.append(s['id'])
                with self.lock:self.persist()
                return imported
            if kind=='new':
                parent=Path(data['folder']).expanduser().resolve()
                if not parent.is_dir():raise ValueError('Select an existing parent folder.')
                name=str(data.get('name','')).strip()
                if not name or len(name)>100 or any(ord(c)<32 for c in name):raise ValueError('Enter a project name (1-100 characters).')
                root=parent/('BehaviorHub_'+dt.datetime.now().strftime('%Y%m%d_%H%M%S')+'_'+uuid.uuid4().hex[:6]);root.mkdir()
                file=root/'BehaviorHub.project.json'
                (root/'data').mkdir()
                atomic(file,dict(schema_version=1,project=dict(id=uuid.uuid4().hex,name=name,created_at=dt.datetime.now().astimezone().isoformat()),
                                 output='data',mice=[],sessions=[],settings={},calibrations={},processing={'idle_seconds':60,'paused':False}))
            elif kind=='open':
                file=Path(data['path']).expanduser().resolve()
            elif kind=='legacy':
                with self.lock:
                    self.persist()
                    if self.project_handle:self.project_handle.close()
                    self.project_handle=None;self.project_file=None;self.db=copy.deepcopy(self.legacy_db)
                return None
            else:raise ValueError('Unknown project action.')
            if file==self.project_file:return str(file)
            # Acquire the new project lock before replacing the current project.
            if not file.is_file():raise ValueError('Project file not found. Copy the complete project folder first.')
            handle=file.with_suffix('.lock').open('a+b')
            try:
                if os.name=='nt':
                    import msvcrt
                    handle.seek(0);msvcrt.locking(handle.fileno(),msvcrt.LK_NBLCK,1)
                db=json.loads(file.read_text(encoding='utf-8'))
                if db.get('schema_version')!=1 or not isinstance(db.get('project'),dict) or not isinstance(db.get('sessions'),list):
                    raise ValueError('Unsupported or invalid Behavior Hub project.')
                root=file.parent
                db['output']=str(inside(root,db.get('output','data')))
                if not Path(db['output']).is_dir():raise ValueError('Project data folder missing. Copy the complete folder, not just the project file.')
                ids=set()
                for s in db['sessions']:
                    if s['id'] in ids:raise ValueError('Duplicate recording IDs in project.')
                    ids.add(s['id']);s['folder']=str(inside(root,s['folder']))
                    for q in s['phases'].values():
                        if Path(q['file']).name!=q['file']:raise ValueError('Invalid video filename.')
                        if q.get('mp4_status')=='processing':q['mp4_status']='queued'
                    if s['status'] in ('starting','recording','finalizing'):
                        s['status']='interrupted';s['error']='Recording did not finish before the service stopped.'
                        q=s['phases'].get('recording',{});q.update(mp4_status='queued',qc='PENDING')
                    s['files_missing']=not Path(s['folder']).is_dir()
                for c in db.get('calibrations',{}).values():c['verification']='Needs verification'
                with self.lock:
                    self.persist()
                    if self.project_handle:self.project_handle.close()
                    self.project_file=file;self.project_handle=handle;self.db=db;self.persist();self.remember_project()
            except Exception:
                if handle is not self.project_handle:handle.close()
                raise
            return str(file)
