// DOM emulation, not a real-browser visual test. Uses only an isolated synthetic service.
import {JSDOM,VirtualConsole} from 'jsdom';
import assert from 'node:assert/strict';
import path from 'node:path';
const base=process.env.VAME_TEST_URL||'http://127.0.0.1:43826', errors=[];
const vc=new VirtualConsole();vc.on('jsdomError',e=>errors.push(String(e)));
async function open(){const html=await (await fetch(base)).text();return new JSDOM(html,{url:base,runScripts:'dangerously',virtualConsole:vc,beforeParse(w){w.HTMLDialogElement.prototype.showModal=function(){this.open=true};w.HTMLDialogElement.prototype.close=function(){this.open=false};w.fetch=(u,opts)=>fetch(new URL(u,base),opts);w.prompt=()=>process.env.VAME_TEST_OUTPUT||path.resolve('validation');w.URL.createObjectURL=()=> 'blob:mock';w.URL.revokeObjectURL=()=>{};}})}
let dom=await open();let d=dom.window.document;
async function closeDom(){dom.window.eval('clearInterval(timer); clearInterval(imageTimer)');await new Promise(r=>setTimeout(r,900));dom.window.close()}
async function wait(fn,ms=20000){const until=Date.now()+ms;while(Date.now()<until){if(fn())return;await new Promise(r=>setTimeout(r,100))}throw Error('Timeout: '+d.getElementById('error').textContent)}
function val(id,v){d.getElementById(id).value=v}
function click(id){d.getElementById(id).click()}
try{
 await wait(()=>d.getElementById('camera').options.length>0);
 await new Promise(r=>setTimeout(r,600));
 assert.equal(d.getElementById('size').value,'1280x720');
 for(let [k,v] of Object.entries({mouse_id:'DOM_TEST',cage:'NO_ANIMAL',weight_g:'25',sex:'F',group:'Saline'}))val(k,v);
 d.getElementById('weight_g').dispatchEvent(new dom.window.Event('input'));
 assert.match(d.getElementById('dose').textContent,/0.3125 mL/);
 d.getElementById('mouseForm').dispatchEvent(new dom.window.Event('submit',{cancelable:true}));
 await wait(()=>d.getElementById('saved').textContent.includes('DOM_TEST'));
 await new Promise(r=>setTimeout(r,200));click('manualFolder');
 await wait(()=>d.getElementById('folder').textContent.includes(process.env.VAME_TEST_OUTPUT||'validation'));
 click('test');val('baseline','5');val('post','3');
 assert.match(d.getElementById('baselineLabel').textContent,/seconds/);
 await new Promise(r=>setTimeout(r,300));click('startBaseline');
 await wait(()=>d.getElementById('status').textContent.includes('DOM_TEST'));
 await closeDom();dom=await open();d=dom.window.document;
 await wait(()=>!d.getElementById('injectionPanel').hidden);
 assert.match(d.getElementById('injectionPrompt').textContent,/Saline/);
 assert.match(d.getElementById('injectionPrompt').textContent,/0.3125/);
 val('actual_volume_ml','.3125');val('injectionNote','SIMULATED DOM test');click('injection');
 await wait(()=>!d.getElementById('postPanel').hidden);
 await new Promise(r=>setTimeout(r,300));click('startPost');
 await wait(()=>d.getElementById('status').textContent.includes('Both phases saved'));
 assert.doesNotMatch(d.body.textContent,/[\u4e00-\u9fff]/);
 assert.equal(d.getElementById('stop').disabled,true);
 assert.match(d.getElementById('qc').textContent,/MP4: queued/);click('processNow');await wait(()=>d.getElementById('qc').textContent.split('MP4: verified').length===3);
 assert.match(d.getElementById('qc').textContent,/MP4: verified/);
 assert.match(d.getElementById('qc').textContent,/Observed FPS: 25.000/);
 assert.match(d.getElementById('qc').textContent,/MKV retained/);
 let button=label=>[...d.querySelectorAll('#history button')].find(x=>x.textContent===label);
 button('Edit').click();
 assert.equal(d.getElementById('editDialog').open,true);
 val('edit_weight_g','24');val('edit_notes','Corrected in UI test');val('edit_reason','Test correction');
 d.getElementById('editForm').dispatchEvent(new dom.window.Event('submit',{cancelable:true}));
 await wait(()=>!d.getElementById('editDialog').open);
 await wait(()=>d.getElementById('details').textContent.includes('Corrected in UI test'));
 await new Promise(r=>setTimeout(r,300));button('Delete').click();
 await wait(()=>d.querySelectorAll('#history tr').length===0);
 click('showDeleted');await wait(()=>!!button('Restore'));
 await new Promise(r=>setTimeout(r,300));button('Restore').click();
 await wait(()=>!!button('Edit'));
 assert.equal(errors.length,0,errors.join('\n'));
 console.log('DOM_OK: form, dose, folder path, test units, two stages, injection prompt, reload recovery, controls');
}finally{await closeDom()}
