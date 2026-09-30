'use strict';
// A shared, text-safe renderer. It only displays an already identity-bound graph.
window.renderBioOrient = function(graph, root, sourceUrl) {
  root.replaceChildren();
  const make=(tag,text)=>{const n=document.createElement(tag);if(text!==undefined)n.textContent=text;return n;};
  const heading=make('h2',graph.layer_kind==='component_chemistry'?'Component chemistry · local chirality checks':'Bio-Orient · sides and surroundings');root.append(heading);
  root.append(make('p','Structural availability, orientational availability, and functional availability are separate evidence dimensions.'));
  const binding=make('details'),bindingTitle=make('summary','Exact model, assembly, conformer and coordinate identity');
  binding.append(bindingTitle,make('pre',JSON.stringify(graph.binding,null,2)));root.append(binding);
  if(sourceUrl){const link=make('a',graph.layer_kind==='component_chemistry'?'Open exact chemistry source':'Open exact Bio-Orient source');link.href=sourceUrl;root.append(link);}
  if(graph.membrane){const m=graph.membrane,a=m.side_assignment||{state:'unknown'};root.append(make('p',`Modeled membrane: ${m.context||'declared context'} · sidedness ${a.state} · positive face ${a.positive_compartment||'unsigned'} · negative face ${a.negative_compartment||'unsigned'}`));const evidence=make('details');evidence.append(make('summary','Mapped compartment markers and source'),make('pre',JSON.stringify(a,null,2)));root.append(evidence);}
  const selectors=make('div');selectors.className='bio-controls';
  const sideLabel=make('label','Surface patch '),side=make('select');side.setAttribute('aria-label','Bio-Orient surface patch');
  side.append(new Option('Choose a patch',''));for(const s of graph.sides||[])side.append(new Option(s.side_id,s.side_id));sideLabel.append(side);
  const neighborLabel=make('label','Neighbor '),neighbor=make('select');neighbor.setAttribute('aria-label','Bio-Orient neighbor');neighbor.append(new Option('Select a patch first',''));neighbor.disabled=true;neighborLabel.append(neighbor);
  const componentLabel=make('label','Chemical component '),component=make('select');component.setAttribute('aria-label','Bio-Orient chemical component');
  component.append(new Option('Choose a component',''));for(const c of graph.components||[])component.append(new Option(c.chemical_component_id+' · '+c.component_id,c.component_id));componentLabel.append(component);
  const atomLabel=make('label','Atom '),atom=make('select');atom.setAttribute('aria-label','Bio-Orient atom');atom.disabled=true;atomLabel.append(atom);
  const normalLabel=make('label'),normals=make('input');normals.type='checkbox';normals.checked=true;normalLabel.append(normals,document.createTextNode(' Surface directions'));
  selectors.append(sideLabel,neighborLabel,componentLabel,atomLabel,normalLabel);root.append(selectors);
  const tools=make('div');tools.className='bio-view-tools';
  const reset=make('button','Reset molecular view'),zoomIn=make('button','Zoom in'),zoomOut=make('button','Zoom out'),focus=make('button','Focus selection');
  for(const button of [reset,zoomIn,zoomOut,focus]){button.type='button';tools.append(button);}root.append(tools);
  const backgroundLabel=make('label','Viewer background '),background=make('select');background.setAttribute('aria-label','Bio-Orient viewer background');background.append(new Option('Light','light'),new Option('Dark','dark'));backgroundLabel.append(background);tools.append(backgroundLabel);
  root.append(make('p','Drag to rotate; scroll to zoom. Colored traces show each protein copy; small gray spheres show all other represented atoms.'));
  const canvas=make('div');canvas.className='bio-orient-viewer';canvas.style.height='420px';canvas.setAttribute('aria-label','Bio-Orient molecular view');root.append(canvas);
  const details=make('div');details.className='bio-details';details.setAttribute('aria-live','polite');root.append(details);
  let viewer=null,initialView=null;
  if(typeof $3Dmol!=='undefined'){
    viewer=$3Dmol.createViewer(canvas,{backgroundColor:'#e8efeb'});viewer.addModel().addAtoms(graph.atoms||[]);
    viewer.rotate(65,'y');viewer.zoomTo();
    initialView=viewer.getView();
    viewer.setClickable({},true,a=>{const id=a.properties?.component_id;if(!id)return;component.value=id;component.onchange();atom.value=String(a.index);atom.onchange();});
  }else canvas.textContent='Molecular rendering unavailable; all evidence remains readable below.';
  const words=value=>typeof value==='string'?value.replaceAll('_',' '):value==null?'Unknown':String(value);
  const number=value=>Number.isFinite(value)?value.toFixed(2):'Unevaluated';
  const show=(title,value)=>{
    details.replaceChildren(make('h3',title));
    const rows=make('dl');
    const row=(label,text)=>rows.append(make('dt',label),make('dd',words(text)));
    if(value.side_id)row('Patch',value.side_id);
    if(value.biological_roles?.length)row('Recorded biological roles',value.biological_roles.map(words).join('; '));
    if(value.equivalence_group)row('Declared equivalence group',value.equivalence_group);
    const availability=value.availability||value.relationship?.availability;
    if(availability){row('Structural availability',availability.structural);row('Orientational availability',availability.orientational);row('Functional availability',availability.functional);}
    if(value.coverage?.declared_residues!==undefined)row('Residue coverage',`${value.coverage.observed_residues} of ${value.coverage.declared_residues} represented; ${value.coverage.missing_residue_keys.length} missing`);
    if(value.compartment_orientation){row('Compartment',value.compartment_orientation.compartment||'Unresolved for this whole patch');row('Patch location',value.compartment_orientation.geometry);row('Membrane evidence',value.compartment_orientation.state);}
    if(value.normal_state)row('Outward direction',value.normal_state==='unresolved'?'Unresolved: no coherent sampled direction':value.normal_state);
    if(value.directional_coherence!==undefined)row('Directional coherence',Number.isFinite(value.directional_coherence)?value.directional_coherence.toFixed(4):'Unevaluated');
    if(value.solvent_exposure){row('Sampled exposure',`${number(value.solvent_exposure.sasa_A2)} Å² in the represented assembly`);row('Isolated chain exposure',`${number(value.solvent_exposure.isolated_chain_sasa_A2)} Å²`);}
    if(value.dynamic_state)row('Recorded state',value.dynamic_state.value?.label||value.dynamic_state.state);
    if(value.neighbors){row('Represented neighbors',`${value.neighbors.length}; ${value.neighbors.filter(e=>e.contact_atom_pairs>0).length} with contacts at 5 Å`);}
    if(value.relationship){const e=value.relationship;row('Neighbor',e.adjacent_object);row('Minimum heavy-atom distance',`${number(e.minimum_distance_A)} Å`);row('Contacts at 5 Å',e.contact_atom_pairs??'Unevaluated');row('Measured occlusion',`${number(e.buried_sasa_A2)} Å² of sampled SASA; physical contact area remains unevaluated`);row('Relative angle',`${number(e.relative_angle_degrees)} degrees`);}
    if(value.component_id)row('Component',value.chemical_component_id+' · '+value.component_id);
    if(value.atom)row('Selected atom',value.atom.atom+' · '+value.atom.elem+' · '+value.atom.properties?.component_id);
    if(Array.isArray(value.chirality)){
      const findings=value.chirality,counts={};for(const f of findings)counts[f.state]=(counts[f.state]||0)+1;
      row('Local chirality',findings.length?Object.entries(counts).map(([state,count])=>`${count} ${words(state)}`).join('; '):'No evaluated center at this selection; inspect reference coverage');
      for(const finding of findings)row('Center '+finding.atom_id,`${words(finding.state)}${finding.reason?' · '+words(finding.reason):''}`);
    }else if(value.chirality)row('Chemistry coverage',value.chirality);
    if(value.coverage?.state)row('Chemical reference coverage',value.coverage.state);
    if(value.sidedness){row('Declared sides',value.sidedness.side_count);row('Scoped relationships',`Multisided: ${words(value.sidedness.multisided)}; homosided: ${words(value.sidedness.homosided)}; heterosided: ${words(value.sidedness.heterosided)}`);}
    const source=value.source||value.patch_source;if(source){row('Source',source.id);row('Source citation',source.citation);}
    details.append(rows);
    if(value.relationship)details.append(make('p','Geometric contact or occlusion does not determine biological interaction or functional availability.'));
    const exact=make('details');exact.append(make('summary','Inspect exact recorded fields'),make('pre',JSON.stringify(value,null,2)));details.append(exact);
  };
  function draw(){
    if(!viewer)return;viewer.setBackgroundColor(background.value==='dark'?'#14231f':'#e8efeb');viewer.removeAllShapes();viewer.removeAllLabels();viewer.setStyle({},{sphere:{scale:.13,color:background.value==='dark'?'#c7d9ce':'#a4b4aa'}});
    // Coordinate-only graphs carry no inferred chemical bonds. Render a mapped
    // C-alpha trace and visible atom spheres without inventing ligand connectivity.
    const palette=['#368375','#537cb1','#a87742','#995d93','#6c8c39','#b96858'];
    const chains=[...new Set((graph.components||[]).filter(c=>c.is_polymer&&c.component_kind==='protein').map(c=>c.chain_id))];
    chains.forEach((chain,index)=>{
      const color=palette[index%palette.length],cas=(graph.atoms||[]).filter(a=>a.chain===chain&&a.atom==='CA').sort((a,b)=>a.resi-b.resi||a.icode.localeCompare(b.icode));
      viewer.addStyle({chain,atom:'CA'},{sphere:{radius:.45,color}});
      for(let i=1;i<cas.length;i++){const a=cas[i-1],b=cas[i],distance=Math.hypot(a.x-b.x,a.y-b.y,a.z-b.z);if(b.resi===a.resi+1&&distance<4.6)viewer.addCylinder({start:a,end:b,radius:.23,color,fromCap:1,toCap:1});}
    });
    const selected=(graph.sides||[]).find(s=>s.side_id===side.value);
    const selectedComponent=(graph.components||[]).find(c=>c.component_id===component.value);
    focus.disabled=!selected&&!selectedComponent&&!atom.value;
    if(selected)viewer.addStyle({serial:selected.atom_indices},{sphere:{scale:.48,color:'#0072B2'}});
    if(selectedComponent)viewer.addStyle({serial:selectedComponent.atom_indices},{sphere:{scale:.36,color:'#D55E00'}});
    const selectedNeighbor=(graph.objects||[]).find(o=>o.object_id===neighbor.value);if(selectedNeighbor?.atom_indices)viewer.addStyle({serial:selectedNeighbor.atom_indices},{sphere:{scale:.24,color:'#009E73'}});
    if(atom.value)viewer.addStyle({serial:[Number(atom.value)]},{sphere:{scale:.65,color:'#D55E00'}});
    for(const f of graph.chirality?.findings||[]){if(f.state!=='reference_mismatch')continue;const a=(graph.atoms||[]).find(a=>a.properties?.component_id===f.component_id&&a.atom===f.atom_id);if(a)viewer.addStyle({serial:[a.serial]},{sphere:{scale:.65,color:'#CC3311'}});}
    if(normals.checked)for(const s of graph.sides||[]){if(!s.centroid||!s.surface_normal)continue;const p=s.centroid,n=s.surface_normal;viewer.addArrow({start:{x:p[0],y:p[1],z:p[2]},end:{x:p[0]+8*n[0],y:p[1]+8*n[1],z:p[2]+8*n[2]},radius:.16,color:'#0072B2'});}
    const m=graph.membrane;
    if(m){const n=m.normal,assignment=m.side_assignment||{state:'unknown'};
      // Shift the drawn slab within its own plane to the assembly centroid;
      // this changes presentation only, never the measured membrane geometry.
      const mean=[0,0,0];for(const a of graph.atoms||[])for(let i=0;i<3;i++)mean[i]+=a[['x','y','z'][i]]/graph.atoms.length;
      const displacement=mean.map((x,i)=>x-m.center[i]),depth=displacement.reduce((sum,x,i)=>sum+x*n[i],0),c=m.center.map((x,i)=>x+displacement[i]-depth*n[i]);
      const cross=(a,b)=>[a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]];
      let u=cross(n,Math.abs(n[0])<.8?[1,0,0]:[0,1,0]);const length=Math.hypot(...u);u=u.map(x=>x/length);const v=cross(n,u);
      for(const sign of [-1,1]){const p=c.map((value,i)=>value+sign*m.half_thickness*n[i]),text=sign===1?assignment.positive_compartment:assignment.negative_compartment;
        const color=text?(sign===1?'#d7a13b':'#598dbf'):'#998aad';
        const vertices=[[-1,-1],[1,-1],[1,1],[-1,1]].map(([a,b])=>({x:p[0]+42*(a*u[0]+b*v[0]),y:p[1]+42*(a*u[1]+b*v[1]),z:p[2]+42*(a*u[2]+b*v[2])}));
        viewer.addCustom({vertexArr:vertices,faceArr:[0,1,2,0,2,3],normalArr:vertices.map(()=>({x:n[0],y:n[1],z:n[2]})),color,alpha:.27,doubleSided:true});
        for(let i=0;i<4;i++)viewer.addCylinder({start:vertices[i],end:vertices[(i+1)%4],radius:.12,color,fromCap:1,toCap:1});
        viewer.addLabel('Membrane plane · '+(text||'unsigned face '+(sign===1?'A':'B'))+' · '+assignment.state,{position:{x:p[0]+sign*24*u[0],y:p[1]+sign*24*u[1],z:p[2]+sign*24*u[2]},fontColor:'#18332f',fontSize:12,backgroundColor:'#ffffff',backgroundOpacity:.95,borderColor:color,borderThickness:1,inFront:true});}}
    viewer.resize();viewer.render();
  }
  side.onchange=()=>{component.value='';atom.value='';atom.disabled=true;neighbor.replaceChildren(new Option('Choose a neighbor',''));const s=(graph.sides||[]).find(s=>s.side_id===side.value);neighbor.disabled=!s;if(s){const edges=(graph.edges||[]).filter(e=>e.source_side===s.side_id);for(const e of edges)neighbor.append(new Option(e.adjacent_object,e.adjacent_object));show('Selected side', {...s,neighbors:edges});}draw();};
  neighbor.onchange=()=>{const edge=(graph.edges||[]).find(e=>e.source_side===side.value&&e.adjacent_object===neighbor.value);if(edge)show('Selected side–neighbor relationship',{side_id:side.value,neighbor:(graph.objects||[]).find(o=>o.object_id===neighbor.value),relationship:edge,patch_source:(graph.sides||[]).find(s=>s.side_id===side.value)?.source});draw();};
  component.onchange=()=>{side.value='';neighbor.value='';neighbor.disabled=true;atom.replaceChildren(new Option('Choose an atom',''));const c=(graph.components||[]).find(c=>c.component_id===component.value);atom.disabled=!c;if(c){for(const i of c.atom_indices){const a=graph.atoms[i];atom.append(new Option(a.atom+' · '+a.elem,String(i)));}show('Selected chemical component',{...c,chirality:(graph.chirality?.findings||[]).filter(f=>f.component_id===c.component_id),coverage:(graph.chirality?.component_coverage||[]).find(f=>f.component_id===c.component_id)});}draw();};
  atom.onchange=()=>{const a=graph.atoms?.[Number(atom.value)];if(a&&atom.value)show('Selected atom',{atom:a,chirality:(graph.chirality?.findings||[]).filter(f=>f.component_id===component.value&&f.atom_id===a.atom)});draw();};
  normals.onchange=draw;
  background.onchange=draw;
  reset.onclick=()=>{if(viewer){viewer.setView(initialView);viewer.render();}};
  zoomIn.onclick=()=>{if(viewer){viewer.zoom(1.35);viewer.render();}};
  zoomOut.onclick=()=>{if(viewer){viewer.zoom(1/1.35);viewer.render();}};
  focus.onclick=()=>{if(!viewer)return;const s=(graph.sides||[]).find(s=>s.side_id===side.value),c=(graph.components||[]).find(c=>c.component_id===component.value);const indices=atom.value?[Number(atom.value)]:c?.atom_indices||s?.atom_indices;if(indices?.length){viewer.zoomTo({serial:indices});viewer.render();}};
  if(!viewer)for(const button of [reset,zoomIn,zoomOut,focus])button.disabled=true;
  show('Sidedness and chemistry coverage',{sidedness:graph.sidedness_descriptors,chirality:graph.chirality?.state,coverage:graph.chirality?.component_coverage,recorded_state_evidence:graph.recorded_state_evidence,limitations:graph.limitations});draw();
};
