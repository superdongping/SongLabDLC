// Serial, debounced focus updates: keep only the latest slider choice.
let focusContext='', focusCaps=null, focusValue=null;
let focusPending=false, focusInFlight=false, focusTimer=null, focusRevision=0, focusMessage='';
function focusSelection(){
  const mode=$('focusMode').value;
  return mode==='manual'?{mode,value:focusValue??0}:{mode};
}
function resetFocusControls(){
  clearTimeout(focusTimer);focusPending=false;focusRevision++;focusCaps=null;focusMessage='';
  const saved=state.focus_profiles?.[$('camera').value]?.focus;
  const name=$('camera').selectedOptions[0]?.dataset.name||'';
  $('focusMode').value=saved?.mode||(/c920/i.test(name)?'auto':'device');
  focusValue=saved?.value??null;
}
function acceptFocusCaps(report){
  if(report?.camera_id!==$('camera').value||!Number.isInteger(report.step))return;
  focusCaps=report;
  const slider=$('focusSlider');slider.min=report.minimum;slider.max=report.maximum;slider.step=report.step;
  if(focusValue===null)focusValue=report.value;
  slider.value=focusValue;
}
function renderFocus(){
  const context=(state.project?.id||'legacy')+'|'+$('camera').value;
  if(focusContext!==context){focusContext=context;resetFocusControls()}
  acceptFocusCaps(state.focus_status);
  const recording=!!state.active&&state.active!=='preview';
  $('focusMode').disabled=recording||!!state.synthetic||!$('camera').value;
  $('focusSlider').disabled=recording||!!state.synthetic||$('focusMode').value!=='manual'||!focusCaps||!(focusCaps.capabilities&2);
  $('focusStatus').textContent=focusInFlight||focusPending?'Adjusting focus...':focusMessage;
  $('focusStatus').hidden=!$('focusStatus').textContent;
}
function scheduleFocus(){
  if(state.active&&state.active!=='preview')return;
  focusRevision++;focusPending=true;focusMessage='';
  clearTimeout(focusTimer);focusTimer=setTimeout(flushFocus,300);render();
}
async function flushFocus(){
  if(!focusPending)return;
  if(busy||focusInFlight){focusTimer=setTimeout(flushFocus,100);return}
  const revision=focusRevision,context=focusContext;
  focusPending=false;focusInFlight=true;render();
  await action(async()=>{
    const current=()=>revision===focusRevision&&context===focusContext;
    try{
      if($('focusMode').value==='manual'&&!focusCaps){
        const report=(await api('focus',{...capture(),focus:{mode:'device'},action:'read'})).result;
        if(!current())return;
        acceptFocusCaps(report);
        if(!focusCaps||(focusCaps.capabilities&2)===0)throw Error('Manual focus is not available for this camera.');
      }
      if(!current())return;
      const request=capture();clearReference();
      if(request.focus.mode==='device')await api('preview',request);
      else{
        const report=(await api('focus',{...request,action:'apply'})).result;
        if(!current())return;
        acceptFocusCaps(report);
        if(state.project)await api('focus',{...request,action:'save'});
      }
      if(current())focusMessage='';
    }catch(e){
      if(current())focusMessage='Could not adjust focus. Try again.';
      // Detailed diagnostics remain available through the service; keep the UI concise.
      throw Error('Could not adjust focus. Try again.');
    }
  });
  focusInFlight=false;render();
  if(focusPending){clearTimeout(focusTimer);focusTimer=setTimeout(flushFocus,100)}
}
$('focusMode').onchange=scheduleFocus;
$('focusSlider').oninput=()=>{focusValue=Number($('focusSlider').value);scheduleFocus()};
