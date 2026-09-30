import {JSDOM,VirtualConsole} from 'jsdom';
import assert from 'node:assert/strict';
const base='http://127.0.0.1:43836',errors=[],vc=new VirtualConsole();vc.on('jsdomError',e=>errors.push(String(e)));
const dom=new JSDOM(await(await fetch(base)).text(),{url:base,runScripts:'dangerously',virtualConsole:vc,beforeParse(w){w.fetch=(u,o)=>String(u).endsWith('/api/project-folder')?Promise.resolve({ok:true,json:async()=>({result:process.env.BEHAVIORHUB_TEST_OUTPUT})}):fetch(new URL(u,base),o);w.prompt=()=>process.env.BEHAVIORHUB_TEST_OUTPUT;w.HTMLDialogElement.prototype.showModal=function(){this.open=true};w.HTMLDialogElement.prototype.close=function(){this.open=false};w.HTMLCanvasElement.prototype.getContext=()=>({save(){},restore(){},scale(){},closePath(){},clearRect(){},beginPath(){},moveTo(){},lineTo(){},rect(){},arc(){},stroke(){}});w.URL.createObjectURL=()=> 'blob:test';w.URL.revokeObjectURL=()=>{}}});
let d=dom.window.document,$=id=>d.getElementById(id);async function wait(fn){let until=Date.now()+30000;while(Date.now()<until){if(fn())return;await new Promise(r=>setTimeout(r,100))}throw Error('Timeout '+$('error').textContent)}
try{
 await wait(()=>$('presetNote').textContent.includes('imported'));
 await wait(()=>dom.window.eval('busy===false'));
 $('userGuideOpen').click();assert.equal($('userGuide').open,true);assert.match($('userGuide').textContent,/choose Manual and drag the slider/);$('userGuideClose').click();assert.equal($('userGuide').open,false);
 assert.equal($('size').value,'1280x720');assert.equal($('mouse_id').value,'');
 for(let [assay,min] of Object.entries({OFT:6,NPR:6,ZERO_MAZE:6,Y_MAZE:8,FST:5,TST:6})){ $('assay').value=assay;$('assay').dispatchEvent(new dom.window.Event('change'));assert.equal(Number($('duration').value),min)}
 $('newProject').click();await wait(()=>$('projectName').textContent!== 'Legacy Library');
 $('assay').value='OFT';$('test').click();$('duration').value=2;
 await new Promise(r=>setTimeout(r,500));$('start').click();
 await wait(()=>$('status').textContent.includes('Recording saved'));
 assert.match($('qc').textContent,/MP4: queued/);await wait(()=>dom.window.eval('busy===false'));$('processNow').click();await wait(()=>$('qc').textContent.includes('MP4: verified'));assert.match($('qc').textContent,/MP4: verified/);assert.match($('qc').textContent,/Measured FPS: 25.000/);
 $('shape').value='circle';$('shape').dispatchEvent(new dom.window.Event('change'));assert.equal($('calStart'),null);
 $('previewStart').click();await wait(()=>dom.window.eval('state.preview'));await new Promise(r=>setTimeout(r,1000));$('referenceCapture').click();await wait(()=>!$('reference').hidden);$('referenceClear').click();assert.equal($('reference').hidden,true);$('previewStop').click();await wait(()=>!dom.window.eval('state.active'));await wait(()=>dom.window.eval('busy===false'));
 let button=name=>[...d.querySelectorAll('#history button')].find(b=>b.textContent===name);
 button('Edit').click();$('edit_mouse_id').value='OptionalMouse';$('editReason').value='UI test';$('editForm').dispatchEvent(new dom.window.Event('submit',{cancelable:true}));
 await wait(()=>!$('editDialog').open);await wait(()=>$('details').textContent.includes('OptionalMouse'));
 await new Promise(r=>setTimeout(r,300));button('Delete').click();await wait(()=>d.querySelectorAll('#history tr').length===0);
 $('deleted').checked=true;$('deleted').dispatchEvent(new dom.window.Event('input'));await wait(()=>!!button('Restore'));
 await new Promise(r=>setTimeout(r,300));button('Restore').click();await wait(()=>!!button('Edit'));
 assert.doesNotMatch(d.body.textContent,/[\u4e00-\u9fff]/);assert.equal(errors.length,0,errors.join('\n'));
 console.log('BEHAVIORHUB_DOM_OK: presets, blank metadata, 720p, timed recording, MP4/QC, guide control, edit/delete/restore, English');
}finally{dom.window.eval('clearInterval(timer);clearInterval(imageTimer)');await new Promise(r=>setTimeout(r,500));dom.window.close()}
