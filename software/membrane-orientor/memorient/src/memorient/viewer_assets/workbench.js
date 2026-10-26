'use strict';
// Same-origin Workbench custody/jobs. No filesystem paths or executable fields.
function configuration(){
  if(!scene||molecularLocked)throw Error('Place a receptor on a geometric surface first.');
  return {surface:scene.surface_spec,placement:scene.placement_spec,organism:$('organism').value,
    membrane_layer:$('envelope').value,sources:currentSource?[currentSource]:[],synthetic:bank[activeBank]?.synthetic===true};
}
async function restoreConfiguration(pdb,spec){
  if(!spec?.surface||!spec?.placement)throw Error('Missing surface or placement settings');
  await addProtein(pdb,'Reopened local scene','receptor',{synthetic:spec.synthetic===true});
  $('shape').value=spec.surface.kind; molecularLocked=false; $('sizeLabel').textContent={plane:'Width, height (Å)',sphere:'Radius (Å)',cylinder:'Radius, length (Å)',ellipsoid:'Semiaxes a, b, c (Å)'}[spec.surface.kind];
  $('size').value=spec.surface.size.join(', ');$('thickness').value=spec.surface.thickness;
  const p=spec.placement;
  for(const key of ['u','v','side','depth'])$(key).value=p[key]??(key==='side'?1:0);
  $('axis').value=p.axis.join(', ');$('anchor').value=(p.anchor||[]).join(', ');$('spin').value=(p.spin||0)*180/Math.PI;
  $('organism').value=spec.organism||'';$('envelope').value=spec.membrane_layer||'';
  revision++;await place(true);
}
$('jsonDownload').onclick=()=>{try{const doc={schema_version:'memorient.project.v1',pdb:source,configuration:configuration(),
  pins:pins.map(p=>({pdb:p.source_pdb,surface:p.surface_spec,placement:p.placement_spec})),camera:viewer.getView()};
  save('memorient-project.json',JSON.stringify(doc,null,2),'application/json');}catch(e){status(e.message,true);}};
$('projectImport').onchange=async()=>{try{const file=$('projectImport').files[0];if(!file)return;if(file.size>10000000)throw Error('Scene limit: 10 MB');
  const doc=JSON.parse(await file.text());if(doc.schema_version!=='memorient.project.v1'||!Array.isArray(doc.pins)||doc.pins.length>20)throw Error('Unsupported saved scene');
  const restored=[];for(const p of doc.pins){const response=await localFetch('/place',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify(p)});const r=await response.json();if(!response.ok)throw Error(r.error);restored.push({...r,source_pdb:p.pdb});}
  await restoreConfiguration(doc.pdb,doc.configuration);pins=restored;render(scene,false);$('pinInfo').textContent=`${pins.length} reopened placements.`;
  if(Array.isArray(doc.camera)&&doc.camera.length===8&&doc.camera.every(Number.isFinite)){viewer.setView(doc.camera);viewer.render();}
}catch(e){status(e.message,true);}};

let wbBusy=false;
let portableModels=null;
async function wbApi(path,data,method='POST'){
  const session=await localSession;if(!session)throw Error('Open Organism Orient inside Workbench for saved jobs.');
  const raw=data instanceof ArrayBuffer;
  const response=await fetch(path,{method:data===undefined?'GET':method,headers:data===undefined?{}:{'X-Yauvi-Token':session.token,'Content-Type':raw?'application/octet-stream':'application/json'},body:data===undefined?undefined:raw?data:JSON.stringify(data)});
  const result=await response.json();if(!response.ok)throw Error(result.error||'Workbench request failed');return result;
}
async function ingest(id,role,file){
  if(!file||file.size>250*1024*1024||!file.size)throw Error('Choose a nonempty local file under 250 MiB.');
  const raw=await file.arrayBuffer(),hash=Array.from(new Uint8Array(await crypto.subtle.digest('SHA-256',raw)),x=>x.toString(16).padStart(2,'0')).join('');
  const upload=await wbApi(`/api/analyses/${id}/ingest`,{role,file_name:file.name,size:raw.byteLength,sha256:hash});
  for(let offset=0,index=0;offset<raw.byteLength;offset+=4*1024*1024,index++)await wbApi(`/api/ingests/${upload.upload_id}/${index}`,raw.slice(offset,offset+4*1024*1024),'PUT');
  await wbApi(`/api/ingests/${upload.upload_id}/finalize`,{});
}
async function submit(kind,files,parameters={}){
  if(wbBusy)return;wbBusy=true;$('wbStatus').textContent='Saving local inputs and checking admission…';
  try{
    const id='membrane-'+crypto.randomUUID();
    await wbApi('/api/analyses',{analysis_id:id,analysis_type:kind,question:kind==='membrane_preparation'?'Prepare and minimize reviewed molecular inputs; stop before equilibration.':'Inspect a declared membrane model and its source records.',subject_id:$('organism').value||'unknown'});
    for(const [role,file] of files)await ingest(id,role,file);
    if(Object.keys(parameters).length)await wbApi(`/api/analyses/${id}/parameters`,{parameters});
    const check=await wbApi(`/api/analyses/${id}/validate`,{});
    if(!check.valid)throw Error('Held: '+check.checks.filter(c=>!c.ok).map(c=>c.detail).join(' · '));
    await wbApi(`/api/analyses/${id}/run`,{experimental_acknowledged:true});
    $('wbStatus').textContent='Queued locally. Equilibration will not run.';await refreshWorkbench();
  }catch(e){$('wbStatus').textContent=e.message;}finally{wbBusy=false;}
}
$('wbSave').onclick=()=>{try{const spec=configuration();submit('membrane_surface_placement',[
  ['structure',new File([source],'receptor.pdb')],['placement_spec',new File([JSON.stringify(spec)],'placement.json')]]);}catch(e){$('wbStatus').textContent=e.message;}};
$('wbPrepare').onclick=()=>{if(!$('wbAck').checked){$('wbStatus').textContent='Review the profile and matching molecular inputs, then acknowledge the preparation boundary.';return;}
  submit('membrane_preparation',[['preparation_profile',$('wbProfile').files[0]],['molecular_topology',$('wbTop').files[0]],['molecular_coordinates',$('wbGro').files[0]]],{threads:Number($('wbThreads').value)});};
$('wbImport').onclick=()=>submit('membrane_result_import',[['membrane_bundle',$('wbBundle').files[0]]]);
function selectedResult(){const value=$('wbRuns').value;if(!value)throw Error('Select a recorded run first.');return value.split('/');}
function recordedControls(result,doc){
  $('viewMode').value='recorded';$('shape').value=doc.membrane_geometry?.kind||'plane';
  const box=result.box_angstrom||(result.box_nm||[]).map(x=>x*10);
  $('sizeLabel').textContent='Recorded periodic box width, height (Å)';$('size').value=box.slice(0,2).join(', ');
  $('thickness').value='';
}
async function openResult(base,doc){
  if(doc.artifacts.configuration&&doc.artifacts.source){
    const [pdb,spec]=await Promise.all([fetch(base+'/source').then(r=>{if(!r.ok)throw Error('Source unavailable');return r.text();}),wbApi(base+'/configuration')]);
    await restoreConfiguration(pdb,spec);status(`${doc.stage} · reopened exact recorded placement. Native sidedness unknown.`);return;
  }
  const result=await wbApi(base+'/scene');
  if(!Array.isArray(result.parts))throw Error('No supported molecular display; full system remains available in the result ZIP.');
  recordedControls(result,doc);
  viewer.removeAllModels();viewer.removeAllShapes();boundsModel=null;explicitModels=null;scene=null;pins=[];
  portableModels=[];$('runLayers').replaceChildren();
  for(const part of result.parts){const model=viewer.addModel(part.pdb,'pdb');
    const protein=['ompf','fhua'].includes(part.id);
    model.setStyle({},protein?{cartoon:{color:part.color},line:{color:part.color}}:{line:{color:part.color||'#82d9c5'}});
    portableModels.push({model,part,protein});
    const label=document.createElement('label'),toggle=document.createElement('input');toggle.type='checkbox';toggle.checked=true;
    toggle.onchange=()=>{if(toggle.checked)model.show();else model.hide();viewer.render();};
    label.append(toggle,document.createTextNode(part.label||({ompf:'OmpF',fhua:'FhuA',outer_lps:'Outer LPS',inner_pl:'Inner phospholipids',divalent_ions:'Source Ca / Mg'}[part.id])||part.id));$('runLayers').append(label);
  }
  viewer.zoomTo();viewer.render();molecularLocked=true;setPlacementLocked(true);$('resetView').disabled=false;
  $('viewCaption').textContent=`${doc.stage} · exact saved molecular coordinates; water and ions remain in the full system. Shape editing is locked.`;
  setLegend([[doc.stage,'#82d9c5'],['No equilibration','#f1d59c']]);
  status(`${result.full_atoms?.toLocaleString()||'Recorded'} atoms · ${doc.profile_id||'imported result'} · no biological validation`);
  $('proteinEvidence').replaceChildren();for(const [index,item] of (doc.sources||[]).entries()){const url=typeof item==='string'?item:item.url;if(typeof url==='string'&&url.startsWith('https://'))$('proteinEvidence').append(sourceLink('Source '+(index+1),url));}for(const assumption of doc.assumptions||[])paragraph($('proteinEvidence'),assumption);
}
$('wbOpen').onclick=async()=>{try{const [id,run]=selectedResult(),base=`/api/analyses/${id}/membrane-result/${run}`;await openResult(base,await wbApi(base+'/manifest'));}catch(e){$('wbStatus').textContent=e.message;}};
$('wbExport').onclick=()=>{try{const [id,run]=selectedResult(),a=document.createElement('a');a.href=`/api/analyses/${id}/membrane-result/${run}/bundle`;a.download='MEMBRANE_RESULT.zip';a.click();}catch(e){$('wbStatus').textContent=e.message;}};
const geometryStyleChange=$('proteinStyle').onchange;
$('proteinStyle').onchange=()=>{if(!portableModels)return geometryStyleChange();
  for(const {model,part,protein} of portableModels){if(protein)model.setStyle({},$('proteinStyle').value==='atoms'?{stick:{radius:.12,color:part.color}}:{cartoon:{color:part.color},line:{color:part.color}});}viewer.render();};
async function refreshWorkbench(){
  try{const [cases,jobs]=await Promise.all([wbApi('/api/membrane-cases'),wbApi('/api/jobs')]);
    const selected=$('wbRuns').value,options=[];
    for(const c of cases){const snapshot=await wbApi('/api/analyses/'+c.analysis_id);for(const r of snapshot.analysis.runs||[])options.push(new Option(`${c.question} · ${r.run_id} · ${r.status}`,c.analysis_id+'/'+r.run_id));}
    $('wbRuns').replaceChildren(new Option('Choose a recorded membrane run',''),...options);$('wbRuns').value=selected;
    const ids=new Set(cases.map(c=>c.analysis_id));$('wbJobs').replaceChildren();
    for(const job of jobs.filter(j=>ids.has(j.analysis_id)).slice(0,12)){const line=document.createElement('p');line.className='hint';line.textContent=`${job.analysis_id} · ${job.state} · ${job.progress?.message||job.error||''}`;
      if(['running','queued'].includes(job.state)){const button=document.createElement('button');button.textContent='Cancel job';button.onclick=async()=>{try{await wbApi('/api/jobs/'+job.job_id+'/cancel',{});await refreshWorkbench();}catch(e){$('wbStatus').textContent=e.message;}};line.append(button);}
      $('wbJobs').append(line);}
    if(!jobs.some(j=>ids.has(j.analysis_id)&&['queued','running'].includes(j.state))&&$('wbStatus').textContent.startsWith('Queued locally'))$('wbStatus').textContent='No membrane jobs are active. Inspect the recorded outcomes below.';
  }catch(e){$('wbStatus').textContent=e.message;}
}
$('wbRefresh').onclick=refreshWorkbench;
(async()=>{const session=await localSession;if(session){$('wbControls').hidden=false;$('wbHint').textContent='Inputs and results stay in this Workbench workspace. Preparation stops before equilibration.';await refreshWorkbench();setInterval(()=>{if(!document.hidden&&!wbBusy)refreshWorkbench();},5000);}else{
  const r=await localFetch('/connection').then(r=>r.json());if(r.has_result){const doc=await localFetch('/result').then(r=>r.json());const result=await localFetch('/result/scene').then(r=>r.json());if(result.pdb){scene=result;render(result,true);status(doc.stage+' · recorded result');}else if(result.parts){recordedControls(result,doc);viewer.removeAllModels();for(const part of result.parts){const m=viewer.addModel(part.pdb,'pdb');m.setStyle({},{line:{color:part.color}});}viewer.zoomTo();viewer.render();molecularLocked=true;setPlacementLocked(true);$('resetView').disabled=false;status(doc.stage+' · exact saved coordinates; no equilibration');}}}
})().catch(e=>status(e.message,true));
