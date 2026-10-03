"""Validated source framing shared by preview and the recorded stream."""
import math

def transform_settings(value=None):
    value = {} if value is None else value
    if not isinstance(value, dict):
        raise ValueError('Invalid video framing settings.')
    result = {}
    for key, default, low, high in [('zoom',1,1,4),('pan_x',0,-1,1),('pan_y',0,-1,1)]:
        v=value.get(key,default)
        if isinstance(v,bool) or not isinstance(v,(float,int)) or not math.isfinite(v) or not low<=v<=high:
            raise ValueError('Invalid '+key+' setting.')
        result[key]=float(v)
    for key in ['flip_h','flip_v']:
        v=value.get(key,False)
        if not isinstance(v,bool):raise ValueError('Invalid '+key+' setting.')
        result[key]=v
    if result['zoom']==1:result.update(pan_x=0.,pan_y=0.)
    return result

def framing_filter(size, value=None):
    t=transform_settings(value)
    width,height=map(int,size.split('x'))
    filters=[]
    if t['zoom']>1:
        w=max(2,int(width/t['zoom']/2)*2);h=max(2,int(height/t['zoom']/2)*2)
        x=int((width-w)*(t['pan_x']+1)/4)*2
        y=int((height-h)*(t['pan_y']+1)/4)*2
        filters.extend([f'crop={w}:{h}:{x}:{y}',f'scale={width}:{height}:flags=bilinear'])
    if t['flip_h']:filters.append('hflip')
    if t['flip_v']:filters.append('vflip')
    return ','.join(filters) or 'null'
