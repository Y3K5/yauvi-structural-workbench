/* Public entry point; the record renderer is shared with the full Research Home. */
(async function(){
 'use strict';
 const view=document.getElementById('function-view'),message=document.getElementById('function-message');
 const esc=v=>String(v??'').replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
 let token;
 const error=text=>{message.textContent=text;message.setAttribute('role','alert');};
 async function api(path,body){
  const options=body===undefined?{}:{method:'POST',headers:{'Content-Type':'application/json','X-Yauvi-Token':token},body:JSON.stringify(body)};
  const response=await fetch(path,options),value=await response.json();
  if(!response.ok)throw Error(value.error||'Local request unavailable');return value;
 }
 const records=window.WorkbenchProteinFunctions({esc,api,error});
 async function open(key,study=null){
  const url=new URL(location.href);url.pathname='/function-records/';url.search='';url.searchParams.set('protein',key);if(study)url.searchParams.set('study',study);history.replaceState(null,'',url);
  message.textContent='Source-derived protein record. Save proposals explicitly; no method executes here.';
  await records.show(key,'#function-view');
 }
 async function browse(){
  history.replaceState(null,'','/function-records/');document.getElementById('function-history').hidden=true;
  view.innerHTML=records.explorer();message.textContent='Reading the locally configured catalogue.';await records.mountExplorer(open);
 }
 document.getElementById('browse-functions').onclick=browse;
 try{
  token=(await api('/api/session')).token;
  const q=new URLSearchParams(location.search),study=q.get('study');
  if(study){
   const loaded=await api('/api/protein-functions/investigations/'+encodeURIComponent(study));
   await open(loaded.packet.function_record.identity.key,study);
   const saved=view.querySelector('[data-open-investigation="'+study+'"]');if(saved)await saved.onclick();
   const host=document.getElementById('function-history');host.hidden=false;
   host.innerHTML='<h2>Saved native investigation</h2><p>'+esc(loaded.study.title)+' · revision '+loaded.study.revision+'</p><p>'+esc(loaded.study.question)+'</p><p>Exact revision '+esc(loaded.study.revision_sha256)+'</p><p>Source state: '+esc(loaded.source_state)+'. Earlier revisions remain available by their exact hashes in the native API.</p><a href="/api/studies/'+encodeURIComponent(study)+'/bundle?revision='+loaded.study.revision_sha256+'">Export this exact native dossier</a><details><summary>Exact native record and retained snapshot</summary><pre>'+esc(JSON.stringify(loaded,null,2))+'</pre></details>';
  }else if(q.get('protein'))await open(q.get('protein'));
  else await browse();
 }catch(e){error(e.message+'. See the function-record configuration guide. No scientific work was started.');}
})();
