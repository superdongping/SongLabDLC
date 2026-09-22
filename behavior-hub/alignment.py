"""V2-style manual geometry selection + KLT features + robust affine tracking."""
import base64
import copy
import datetime as dt
import threading
import time
from collections import deque
import cv2
import numpy as np
cv2.setNumThreads(1)

class AlignmentEngine:
    def __init__(self):
        self.lock=threading.RLock();self.stop_event=threading.Event();self.wake=threading.Event()
        self.enabled=False;self.latest=None;self.serial=0;self.consumed=0;self.frozen=None
        self.gray=None;self.features=None;self.geometry=None;self.kind='rectangle';self.context={}
        self.capture_times=deque(maxlen=100);self.tracking_times=deque(maxlen=100)
        self.info={'status':'Off','verification':'Needs verification'}
        self.thread=threading.Thread(target=self.loop,daemon=True);self.thread.start()
    def begin(self,context):
        with self.lock:
            self.enabled=True;self.context=dict(context);self.features=None;self.geometry=None;self.frozen=None
            self.capture_times.clear();self.tracking_times.clear();self.info={'status':'Live preview — freeze to select arena','verification':'Needs verification'}
    def stop(self):
        with self.lock:self.enabled=False;self.features=None;self.info['status']='Off'
    def close(self):self.stop_event.set();self.wake.set();self.thread.join(3)
    def offer(self,jpeg):
        if not self.enabled:return
        with self.lock:self.latest=jpeg;self.serial+=1;self.capture_times.append(time.monotonic())
        self.wake.set()
    @staticmethod
    def rate(times):return (len(times)-1)/(times[-1]-times[0]) if len(times)>1 and times[-1]>times[0] else 0
    def state(self):
        with self.lock:
            data=copy.deepcopy(self.info);data.update(enabled=self.enabled,target_fps=25,
                capture_fps=round(self.rate(self.capture_times),1) if self.capture_times and time.monotonic()-self.capture_times[-1]<2 else 0,
                tracking_fps=round(self.rate(self.tracking_times),1) if self.tracking_times and time.monotonic()-self.tracking_times[-1]<2 else 0)
            return data
    def freeze(self):
        with self.lock:
            if not self.enabled or not self.latest:raise ValueError('Start calibration and wait for a camera frame.')
            self.frozen=cv2.imdecode(np.frombuffer(self.latest,np.uint8),cv2.IMREAD_COLOR)
            self.features=None;self.info={'status':'Frozen — select arena geometry','verification':'Needs verification'}
            return {'image':'data:image/jpeg;base64,'+base64.b64encode(self.latest).decode(),
                    'width':self.frozen.shape[1],'height':self.frozen.shape[0]}
    def seed(self,gray,geometry):
        mask=np.zeros_like(gray);x,y,w,h=cv2.boundingRect(geometry.astype(np.float32))
        pad=35;cv2.rectangle(mask,(max(0,x-pad),max(0,y-pad)),(min(gray.shape[1]-1,x+w+pad),min(gray.shape[0]-1,y+h+pad)),255,-1)
        features=cv2.goodFeaturesToTrack(gray,300,.01,8,mask=mask)
        if features is None or len(features)<12:raise ValueError('Too few visual features. Improve light/contrast or add fixed markers near the arena.')
        return features
    def select(self,data):
        with self.lock:
            if self.frozen is None:raise ValueError('Freeze a frame before selecting geometry.')
            h,w=self.frozen.shape[:2];kind=data.get('kind')
            if kind=='rectangle':
                points=np.array(data['points'],dtype=np.float32)
                if points.shape!=(4,2) or not np.isfinite(points).all():raise ValueError('Select four finite corner points.')
                if not cv2.isContourConvex(points) or abs(cv2.contourArea(points))<100:raise ValueError('Select corners in clockwise order: top-left, top-right, bottom-right, bottom-left.')
            elif kind=='circle':
                e=data['ellipse'];cx,cy,rx,ry,angle=[float(e[k]) for k in ('cx','cy','rx','ry','angle')]
                if not np.isfinite([cx,cy,rx,ry,angle]).all() or min(rx,ry)<5:raise ValueError('Draw a valid ellipse.')
                t=np.linspace(0,2*np.pi,80,endpoint=False);a=np.deg2rad(angle)
                points=np.column_stack((cx+rx*np.cos(t)*np.cos(a)-ry*np.sin(t)*np.sin(a),cy+rx*np.cos(t)*np.sin(a)+ry*np.sin(t)*np.cos(a))).astype(np.float32)
            else:raise ValueError('Choose rectangle or circle.')
            if np.any(points<0) or np.any(points[:,0]>=w) or np.any(points[:,1]>=h):raise ValueError('Arena geometry must be inside the frame.')
            gray=cv2.cvtColor(self.frozen,cv2.COLOR_BGR2GRAY);features=self.seed(gray,points)
            self.gray=gray;self.features=features;self.geometry=points;self.kind=kind;self.steps=0
            self.info={'status':'Tracking','verification':'Needs verification','points':points.tolist(),'width':w,'height':h}
            self.tracking_times.clear()
            self.profile={'kind':kind,'points':points.tolist(),'width':w,'height':h,'context':self.context,
                          'date':dt.date.today().isoformat(),'verification':'Needs verification',
                          'reference_jpeg':base64.b64encode(cv2.imencode('.jpg',self.frozen)[1]).decode()}
            return self.profile
    def metrics(self,points,w,h):
        center=points.mean(axis=0);offset=(center-np.array([w/2,h/2]))/np.array([w,h])
        if self.kind=='rectangle':
            a,b,c,d=points;top=b-a;bottom=c-d;left=d-a;right=c-b
            angle=float(np.mean([np.degrees(np.arctan2(v[1],v[0])) for v in (top,bottom)]))
            lr=float((np.linalg.norm(left)-np.linalg.norm(right))/max(np.linalg.norm(left),np.linalg.norm(right),1))
            tb=float((np.linalg.norm(top)-np.linalg.norm(bottom))/max(np.linalg.norm(top),np.linalg.norm(bottom),1))
            good=abs(angle)<1.5 and abs(lr)<.08 and abs(tb)<.08 and np.max(np.abs(offset))<.03
            values=dict(rotation_deg=round(angle,2),left_right_asymmetry=round(lr,3),top_bottom_asymmetry=round(tb,3))
        else:
            ellipse=cv2.fitEllipse(points.reshape(-1,1,2));axes=ellipse[1];ratio=min(axes)/max(axes)
            good=ratio>=.97 and np.max(np.abs(offset))<.03;values=dict(axis_ratio=round(ratio,4),ellipse_angle_deg=round(ellipse[2],2))
        values.update(center_offset_percent=[round(float(x)*100,2) for x in offset],alignment_ok=bool(good))
        return values
    def loop(self):
        while not self.stop_event.is_set():
            self.wake.wait(.1);self.wake.clear()
            with self.lock:
                if not self.enabled or self.features is None or self.serial==self.consumed:continue
                self.consumed=self.serial;jpeg=self.latest
                try:
                    frame=cv2.imdecode(np.frombuffer(jpeg,np.uint8),cv2.IMREAD_GRAYSCALE)
                    if frame is None or frame.shape!=self.gray.shape:raise ValueError('Frame size changed. Recalibrate.')
                    cur,valid,_=cv2.calcOpticalFlowPyrLK(self.gray,frame,self.features,None,winSize=(25,25),maxLevel=3)
                    back,valid_back,_=cv2.calcOpticalFlowPyrLK(frame,self.gray,cur,None,winSize=(25,25),maxLevel=3)
                    ok=(valid.ravel()>0)&(valid_back.ravel()>0)&(np.linalg.norm(self.features-back,axis=2).ravel()<3)
                    if np.sum(ok)<12:raise ValueError('Tracking weak — freeze and recalibrate.')
                    transform,inliers=cv2.estimateAffine2D(self.features[ok],cur[ok],method=cv2.RANSAC,ransacReprojThreshold=3,maxIters=200,confidence=.95)
                    if transform is None or inliers is None or np.sum(inliers)<12:raise ValueError('Transform unreliable — recalibrate.')
                    points=cv2.transform(self.geometry[None,:,:],transform)[0]
                    if not np.isfinite(points).all():raise ValueError('Invalid transform — recalibrate.')
                    self.geometry=points;self.gray=frame;self.features=cur[ok][inliers.ravel()>0];self.steps+=1
                    if self.steps%60==0:self.features=self.seed(frame,points)
                    h,w=frame.shape;metrics=self.metrics(points,w,h)
                    self.info.update(status='Alignment OK' if metrics['alignment_ok'] else 'Adjust camera',points=points.tolist(),metrics=metrics,
                                     features=len(self.features),width=w,height=h)
                    if not metrics['alignment_ok']:self.info['verification']='Needs verification'
                    self.tracking_times.append(time.monotonic())
                except Exception as exc:
                    self.features=None;self.info.update(status=str(exc),verification='Needs verification');self.info.pop('metrics',None)
    def confirm(self):
        with self.lock:
            if not self.enabled or self.features is None or not self.tracking_times or time.monotonic()-self.tracking_times[-1]>1 or not self.info.get('metrics',{}).get('alignment_ok'):
                raise ValueError('Obtain a reliable Alignment OK result before confirming.')
            self.info['verification']='Verified this session'
            self.profile.update(reference_jpeg=base64.b64encode(cv2.imencode('.jpg',self.gray)[1]).decode(),points=self.geometry.tolist(),date=dt.date.today().isoformat(),verification='Verified this session')
            return copy.deepcopy(self.profile)
