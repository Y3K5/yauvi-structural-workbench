'use strict';
// Shared by normal analysis views and registered biological cases.
window.RegionExplorer=(()=>{
 let data=null,region=null,indices=[],first=null,measurement=null,source=null,getViewer=null,redraw=null,componentOnly=false;
 const node=id=>document.getElementById(id);
 const text=(tag,value)=>{const e=document.createElement(tag);e.textContent=value;return e;};
 const identity=a=>`${a.chain} / ${a.resn} ${a.resi}${a.icode||''} / ${a.atom} / alt ${a.label_alt_id||'blank'} / source atom ${a.source_atom_id}`;
 const stable=x=>JSON.stringify(x,(key,value)=>value&&typeof value==='object'&&!Array.isArray(value)?Object.fromEntries(Object.keys(value).sort().map(k=>[k,value[k]])):value);
 function clear(){data=null;region=null;indices=[];first=null;measurement=null;componentOnly=false;
  node('regions-panel').hidden=true;for(const id of ['region-cards','region-detail','region-states','region-atom'])node(id).replaceChildren();node('distance-result').textContent='';}
 function atoms(){node('region-atom').replaceChildren(new Option('Choose a represented atom',''),...indices.map(i=>new Option(identity(data.atoms[i]),String(i))));}
 function choose(r){if(!data||!r||stable(r.binding)!==stable(data.regions[0]?.binding))throw Error('Region binding differs from this view');
  region=r;componentOnly=false;first=null;measurement=null;node('distance-result').textContent='';
  indices=[...new Set([...r.atom_indices,...r.ligand_atom_indices])];atoms();const box=node('region-detail');box.replaceChildren(text('h3',r.title));
  box.append(text('p',`${r.basis.replaceAll('_',' ')} · ${r.mapping_state} · observed ${r.coverage.observed} / ${r.coverage.expected??'unknown'} expected`),source(r.source_id));
  if(r.feature_pointer)box.append(text('p','Reference feature '+r.feature_pointer));
  if(r.missing_positions?.length)box.append(text('p','Missing or unmapped positions: '+r.missing_positions.join(', ')));
  if(r.missing_heavy_atom_names===null)box.append(text('p','Chemical atom coverage is unevaluated: deposited reference atom list unavailable.'));
  else if(r.missing_heavy_atom_names?.length)box.append(text('p','Missing represented component atoms: '+r.missing_heavy_atom_names.join(', ')));
  if(r.measurement){box.append(text('p',`${r.measurement.cutoff_A} Å heavy-atom convention · ${r.measurement.contact_residue_count} contacting residues · minimum ${r.measurement.minimum_distance_A?.toFixed(3)??'unevaluated'} Å`));
   for(const p of r.measurement.closest_pairs){const b=text('button',`${p.residue_id}: ${p.distance_A.toFixed(3)} Å`);b.onclick=()=>distance(p.protein_atom_index,p.ligand_atom_index);box.append(b);}}
  box.append(text('p','Residues: '+(r.residue_ids.join(', ')||'No mapped coordinate residues')),text('p',r.limitations.join(' ')));
  const details=document.createElement('details');details.append(text('summary','Exact view binding'),text('pre',JSON.stringify(r.binding,null,2)));box.append(details);
  for(const b of node('region-cards').children)b.setAttribute('aria-pressed',String(b.dataset.region===r.id));redraw();}
 function distance(i,j){const a=data.atoms[i],b=data.atoms[j];if(!a||!b)throw Error('Atom is unavailable in this view');
  const d=Math.hypot(a.x-b.x,a.y-b.y,a.z-b.z);measurement={schema_version:'1.0',binding:region?.binding||data.regions[0]?.binding,atom_indices:[i,j],atoms:[a,b],distance_A:d,interpretation_limit:'Geometric distance only; no bond or affinity assignment.'};
  node('distance-result').textContent=`${identity(a)} → ${identity(b)}: ${d.toFixed(3)} Å. Geometric distance only.`;redraw();}
 function draw(v){if(!data)return;
  if(region){v.addStyle({serial:region.atom_indices},{stick:{radius:.14,color:'#eea22a'},sphere:{scale:.2,color:'#eea22a'}});
   v.addStyle({serial:region.ligand_atom_indices},{stick:{radius:.2,color:'#ab4e99'},sphere:{scale:.3,color:'#ab4e99'}});
   if(componentOnly){v.setStyle({},{sphere:{scale:.03,color:'#dce3df'}});v.setStyle({serial:region.ligand_atom_indices},{stick:{radius:.2},sphere:{scale:.3}});}}
  if(measurement){const [a,b]=measurement.atom_indices.map(i=>data.atoms[i]);v.addLine({start:a,end:b,color:'#224ea0',dashed:true});
   v.addLabel(measurement.distance_A.toFixed(3)+' Å',{position:{x:(a.x+b.x)/2,y:(a.y+b.y)/2,z:(a.z+b.z)/2},backgroundColor:'#ffffff',fontColor:'#224ea0',backgroundOpacity:.9});}}
 function load(view,link,viewer,render){clear();data=view;source=link;getViewer=viewer;redraw=render;node('regions-panel').hidden=false;
  for(const r of view.regions){const b=text('button',`${r.title} · ${r.basis.replaceAll('_',' ')} · ${r.mapping_state}`);b.dataset.region=r.id;b.setAttribute('aria-pressed','false');b.onclick=()=>choose(r);node('region-cards').append(b);}
  if(!view.regions.length)node('region-cards').append(text('p','No exact annotations or measured component contacts are available. Unmapped evidence remains accessible under Sources.'));
  indices=view.atoms.map(a=>a.index);atoms();
  for(const state of view.recorded_states){const card=text('p',`${state.property}: ${JSON.stringify(state.value)} · ${state.basis} · ${state.state}. ${state.interpretation_limit}`);card.append(source(state.source_id));node('region-states').append(card);}
  if(!view.recorded_states.length)node('region-states').append(text('p','No source-supported static-state label is bound to this exact view. Occupancy alone does not establish a catalytic state.'));
  node('region-whole').onclick=()=>{region=null;componentOnly=false;first=null;measurement=null;indices=view.atoms.map(a=>a.index);atoms();node('region-detail').replaceChildren();node('distance-result').textContent='';redraw();const v=getViewer();if(v){v.zoomTo();v.render();}};
  node('region-focus').onclick=()=>{const v=getViewer();if(v&&indices.length){v.zoomTo({serial:indices});v.render();}};
  node('region-ligand').onclick=()=>{if(!region?.ligand_atom_indices.length){node('distance-result').textContent='This annotation has no represented component selection.';return;}componentOnly=true;indices=region.ligand_atom_indices;atoms();redraw();const v=getViewer();if(v){v.zoomTo({serial:indices});v.render();}};
  node('region-export').onclick=()=>{const blob=new Blob([JSON.stringify({schema_version:'1.0',view_binding:region?.binding||{...view,atoms:undefined,chains:undefined,regions:undefined},region,measurement,limitations:view.limitations},null,2)],{type:'application/json'});const url=URL.createObjectURL(blob),a=document.createElement('a');a.href=url;a.download='REGION_SELECTION.json';a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);};
  node('distance-first').onclick=()=>{if(node('region-atom').value==='')return;first=Number(node('region-atom').value);node('distance-result').textContent='First atom: '+identity(view.atoms[first]);};
  node('distance-second').onclick=()=>{if(first===null||node('region-atom').value===''){node('distance-result').textContent='Select a first and second represented atom.';return;}distance(first,Number(node('region-atom').value));};
 }
 return {clear,load,draw,selectAtom(i){if(!data?.atoms[i])return;if(!indices.includes(i)){indices.push(i);atoms();}node('region-atom').value=String(i);}};
})();
