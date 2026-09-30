import {JSDOM,VirtualConsole} from 'jsdom';
import {readFileSync} from 'node:fs';
import assert from 'node:assert/strict';
let html=readFileSync('index.html','utf8').replace('__TOKEN__','test').replace('__FOCUS_UI__',readFileSync('focus_ui.js','utf8'));
let applyDelay=0, rejectApply=false;
const errors=[],calls=[],vc=new VirtualConsole();vc.on('jsdomError',e=>errors.push(String(e)));
const st={project:{id:'p1',name:'Test'},project_file:'test.project.json',settings:{},sessions:[],output:'test',
  recent_projects:[],processing:{jobs:[],config:{idle_seconds:60}},active:null,preview:false,synthetic:false,focus_profiles:{}};
const caps={minimum:0,maximum:250,step:5,capabilities:3,value:30,flags:1,verified:false,camera_id:'@device_one'};
const dom=new JSDOM(html,{runScripts:'dangerously',url:'http://127.0.0.1:43836',virtualConsole:vc,beforeParse(w){
  w.fetch=async(u,o={})=>{
    let path=String(u),body=o.body?JSON.parse(o.body):{};calls.push([path,body]);let data;
    if(path==='/api/state') data=structuredClone(st);
    else if(path==='/api/devices') data={devices:['C920','C920'],details:[{name:'C920',id:'@device_one'},{name:'C920',id:'@device_two'}]};
    else if(path==='/api/presets') data=JSON.parse(readFileSync('presets.json'));
    else if(path==='/api/focus'){
      if(body.action==='read'){st.focus_status={...caps};st.live_capture=null;st.active=null;st.preview=false;}
      if(body.action==='apply'){
        if(applyDelay)await new Promise(r=>setTimeout(r,applyDelay));
        if(rejectApply)return {ok:false,status:400,json:async()=>({error:'readback mismatch'})};
        st.focus_status={...caps,verified:true,flags:body.focus.mode==='manual'?2:1,value:body.focus.value??30};
        st.live_capture={camera:body.camera,camera_id:body.camera_id,size:body.size,fps:body.fps,input_format:body.input_format,focus:body.focus};
        st.active='preview';st.preview=true;
      }
      if(body.action==='save')st.focus_profiles[body.camera_id]={focus:body.focus};
      data={ok:true,result:st.focus_status};
    }else if(path==='/api/preview')return {status:204,ok:true};
    else data={ok:true};
    return {ok:true,status:200,json:async()=>data};
  };
  w.HTMLDialogElement.prototype.showModal=function(){this.open=true};w.HTMLDialogElement.prototype.close=function(){this.open=false};
  w.HTMLCanvasElement.prototype.getContext=()=>({clearRect(){},beginPath(){},moveTo(){},lineTo(){},rect(){},arc(){},stroke(){}});
  w.URL.createObjectURL=()=> 'blob:test';w.URL.revokeObjectURL=()=>{};
}});
const w=dom.window,$=id=>w.document.getElementById(id);
async function wait(fn){const end=Date.now()+5000;while(Date.now()<end){if(fn())return;await new Promise(r=>setTimeout(r,20))}throw Error('UI timeout: '+$('error').textContent)}
async function idle(){await wait(()=>w.eval('!busy&&!focusPending&&!focusInFlight'))}
function change(id,value,event='change'){$(id).value=value;$(id).dispatchEvent(new w.Event(event))}
try{
  await wait(()=>$('camera').options.length===2);await idle();
  for(const id of ['focusRead','focusApply','focusSave','focusNumber','focusZoom','focusZoomOpen'])assert.equal($(id),null,id);
  assert.equal($('focusMode').value,'auto');assert.equal($('focusSlider').disabled,true);
  // Choosing Manual automatically discovers capabilities, applies, opens preview and saves.
  change('focusMode','manual');await idle();
  assert.equal($('focusSlider').disabled,false);assert.equal($('focusSlider').step,'5');assert.equal($('focusSlider').max,'250');
  assert.equal(st.live_capture.focus.mode,'manual');assert.equal(st.live_capture.focus.value,30);
  assert.equal($('focusStatus').hidden,true);
  // A burst of dragging applies only its last value; Start is unavailable until done.
  let n=calls.filter(([p,b])=>p==='/api/focus'&&b.action==='apply').length;
  for(const value of ['35','40','45'])change('focusSlider',value,'input');
  assert.equal($('start').disabled,true);await idle();
  let applied=calls.filter(([p,b])=>p==='/api/focus'&&b.action==='apply');
  assert.equal(applied.length,n+1);assert.equal(applied.at(-1)[1].focus.value,45);
  assert.equal(st.focus_profiles['@device_one'].focus.value,45);
  assert.equal($('start').disabled,false);
  // Dragging during a slow request must eventually apply/save the newest choice.
  applyDelay=500;change('focusSlider','50','input');
  await wait(()=>w.eval('focusInFlight'));
  change('focusSlider','55','input');change('focusSlider','60','input');await idle();applyDelay=0;
  assert.equal(st.live_capture.focus.value,60);assert.equal(st.focus_profiles['@device_one'].focus.value,60);
  assert.equal($('focusSlider').value,'60');
  // Failed application must not persist the rejected setting.
  rejectApply=true;change('focusSlider','65','input');await idle();
  assert.match($('focusStatus').textContent,/Could not adjust/);assert.equal(st.focus_profiles['@device_one'].focus.value,60);
  rejectApply=false;change('focusSlider','70','input');await idle();assert.equal(st.focus_profiles['@device_one'].focus.value,70);
  st.active='recording-test';st.preview=false;await w.eval('poll()');
  for(const id of ['focusMode','focusSlider'])assert.equal($(id).disabled,true,id);
  st.active=null;st.live_capture=null;await w.eval('poll()');change('camera','@device_two');
  assert.equal($('focusMode').value,'auto');
  change('camera','@device_one');assert.equal($('focusMode').value,'manual');
  assert.equal($('focusSlider').value,'70');
  assert.equal(errors.length,0,errors.join('\n'));
  console.log('FOCUS_DOM_OK: simplified UI, automatic capability lookup/apply/save, debounce, latest-value wins, failure handling, recording lock, profile restore');
}finally{w.eval('clearInterval(timer);clearInterval(imageTimer);clearTimeout(focusTimer)');w.close()}
