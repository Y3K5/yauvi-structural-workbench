// Endpoint interpolation is illustrative; it is never a physical trajectory.
let stopStateMotion = () => {};
function buildStateMotion(p) {
  const m=p.motion, el=document.getElementById('state-motion-view');
  const viewer=$3Dmol.createViewer(el,{backgroundColor:bg(),antialias:true});
  viewer.addModel(m.active_pdb,'pdb'); viewer.addModel(m.inactive_pdb,'pdb');
  const slider=document.getElementById('state-motion-slider');
  const region=document.getElementById('state-motion-region');
  const residue=document.getElementById('state-motion-residue');
  const play=document.getElementById('state-motion-play');
  const status=document.getElementById('state-motion-status');
  let timer=null,direction=1;
  const vec=a=>({x:a[0],y:a[1],z:a[2]});
  const fmt=n=>n===null?'unavailable':n.toFixed(1)+'°';
  const difference=(a,b)=>a===null||b===null?'unavailable':(((b-a+540)%360)-180).toFixed(1)+'°';
  function selection(){const r=m.regions.find(x=>x.id===region.value);return m.residues.filter(x=>x.position>=r.start&&x.position<=r.end);}
  function detail(){
    const r=m.residues.find(x=>x.position===Number(residue.value));
    document.getElementById('state-motion-detail').textContent=
      r.resname+' '+r.position+': Cα displacement '+r.displacement.toFixed(2)+' Å after the declared fit.';
    const body=document.getElementById('state-angle-body');body.replaceChildren();
    for(const [key,label] of [['phi','Backbone φ (torsion)'],['psi','Backbone ψ (torsion)'],['chi1','Side-chain χ1 (Asp/Phe only)'],['n_ca_c','N–Cα–C bond angle']]){
      const a=r.active_angles[key],b=r.inactive_angles[key],tr=document.createElement('tr');
      for(const value of [label,fmt(a),fmt(b),key==='n_ca_c'?(a===null||b===null?'unavailable':(b-a).toFixed(1)+'°'):difference(a,b)]){
        const td=document.createElement('td');td.textContent=value;tr.appendChild(td);
      }body.appendChild(tr);
    }
  }
  function draw(){
    const t=Number(slider.value)/100, selected=selection(), ids=selected.map(r=>r.position);
    viewer.removeAllShapes(); viewer.removeAllLabels();
    viewer.setStyle({},{cartoon:{color:'#9CA3AF',opacity:0.16}});
    viewer.setStyle({model:0,resi:ids},{cartoon:{color:'#2A9D8F',opacity:0.4}});
    viewer.setStyle({model:1,resi:ids},{cartoon:{color:'#7B5EA7',opacity:0.4}});
    const focus=Number(residue.value);
    viewer.addStyle({model:0,resi:focus},{stick:{color:'#2A9D8F',radius:0.18}});
    viewer.addStyle({model:1,resi:focus},{stick:{color:'#7B5EA7',radius:0.18}});
    let previous=null;
    for(const r of selected){
      const xyz=r.active.map((x,i)=>x*(1-t)+r.inactive[i]*t),color=r.position===focus?'#E63946':'#E9A23B';
      viewer.addSphere({center:vec(xyz),radius:r.position===focus?0.65:0.28,color,clickable:true,
        callback:()=>{residue.value=String(r.position);detail();draw();}});
      if(previous&&r.position===previous.r.position+1&&
        Math.hypot(...r.active.map((x,i)=>x-previous.r.active[i]))<4.5&&
        Math.hypot(...r.inactive.map((x,i)=>x-previous.r.inactive[i]))<4.5){
        viewer.addCylinder({start:vec(previous.xyz),end:vec(xyz),radius:0.12,color:'#E9A23B',fromCap:1,toCap:1});
      }
      if(r.position===focus){
        if(r.displacement>0.1)viewer.addArrow({start:vec(r.active),end:vec(r.inactive),radius:0.13,color:'#E63946'});
        viewer.addLabel(r.resname+' '+r.position,{position:vec(xyz),fontSize:12,fontColor:'#111',backgroundColor:'#fff',backgroundOpacity:0.9,inFront:true});
      }
      previous={r,xyz};
    }
    status.textContent=t===0?'Active-like endpoint · 2GQG':t===1?'Inactive-like endpoint · 2HYY':
      'Illustrative interpolation · '+slider.value+'% · not an observed structure';
    viewer.render();
  }
  function reset(){
    const selected=selection();
    if(!selected.some(r=>r.position===Number(residue.value)))residue.value=String(selected.reduce((a,b)=>a.displacement>b.displacement?a:b).position);
    draw();viewer.zoomTo({resi:selected.map(r=>r.position)});viewer.zoom(region.value==='dfg'?1.15:0.72);viewer.render();detail();
    const r=m.regions.find(x=>x.id===region.value);
    document.getElementById('state-region-detail').textContent=r.mapped_count+' mapped residues; mean endpoint displacement '+r.mean_displacement.toFixed(2)+' Å. Gaps are left open.';
  }
  stopStateMotion=()=>{if(timer)clearInterval(timer);timer=null;play.textContent='Play illustrative motion';play.setAttribute('aria-pressed','false');};
  play.addEventListener('click',()=>{if(timer){stopStateMotion();return;}play.textContent='Pause';play.setAttribute('aria-pressed','true');
    timer=setInterval(()=>{let n=Number(slider.value)+direction*2;if(n>=100){n=100;direction=-1;}if(n<=0){n=0;direction=1;}slider.value=String(n);draw();},90);});
  slider.addEventListener('input',()=>{stopStateMotion();draw();});
  region.addEventListener('change',()=>{stopStateMotion();reset();});
  residue.addEventListener('change',()=>{const n=Number(residue.value),r=m.regions.find(x=>x.id===region.value);
    if(n<r.start||n>r.end)region.value='all';reset();});
  document.getElementById('state-motion-reset').addEventListener('click',reset);
  document.getElementById('state-motion-export').addEventListener('click',()=>{
    const url=URL.createObjectURL(new Blob([JSON.stringify(m,null,2)],{type:'application/json'}));
    const a=document.createElement('a');a.href=url;a.download='ABL-endpoint-geometry.json';a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);
  });
  document.addEventListener('visibilitychange',()=>{if(document.hidden)stopStateMotion();});
  reset();made[p.id].push(viewer);
}
