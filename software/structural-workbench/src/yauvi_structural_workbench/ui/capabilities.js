/* Native host for the same read-only capability contract used by Research. */
(function(){
  'use strict';
  const e=(tag,text)=>{const n=document.createElement(tag);if(text!==undefined)n.textContent=text;return n;};
  const names={configured:'Named source files present',not_configured:'Provider not configured',source_missing:'Source missing',source_unavailable:'Source unavailable locally',runtime_missing:'Configured interpreter missing',changed_since_view:'Source or configuration changed',not_registered:'Not registered in this distribution'};
  const button=(text,action)=>{const b=e('button',text);b.type='button';b.className='secondary';b.onclick=action;return b;};
  async function get(path){
    const response=await fetch(path,{cache:'no-store',credentials:'same-origin',redirect:'error'});
    const value=await response.json();if(!response.ok)throw Error(value.error||'Local metadata unavailable');return value;
  }
  const field=(root,label,value)=>{const p=e('p');p.append(e('strong',label+' · '),document.createTextNode(value));root.append(p);};
  async function open(){
    const previous=document.activeElement,dialog=e('dialog');dialog.className='capability-dialog';dialog.setAttribute('aria-label','Method capabilities and readiness');
    const heading=e('h2','Method capabilities and readiness'),close=button('Close capability view',()=>dialog.close());
    const intro=e('p','Inspect existing method definitions, local configuration and saved input checks. Source presence is not a runtime or scientific pass.');
    const status=e('p');status.setAttribute('role','status');
    const cards=e('div');cards.className='capability-cards';const detail=e('section');detail.setAttribute('aria-live','polite');
    const refresh=button('Refresh capability status',load);let generation=0;
    dialog.append(close,heading,intro,refresh,status,cards,detail);document.body.append(dialog);
    dialog.addEventListener('close',()=>{generation++;dialog.remove();previous?.focus();},{once:true});dialog.showModal();
    async function inspect(card,analysisId){
      const serial=++generation;detail.replaceChildren();status.textContent='Reading method metadata…';
      try{
        const query=new URLSearchParams({expected_sha256:card.identity_sha256});if(analysisId)query.set('analysis_id',analysisId);
        const value=await get('/api/capabilities/'+encodeURIComponent(card.id)+'?'+query);
        if(!dialog.isConnected||serial!==generation)return;
        const title=e('h3',value.title);title.tabIndex=-1;detail.append(title);
        field(detail,'Local status',names[value.source_state]||'Unknown');
        field(detail,'Owner',value.owner||'Current registry');
        field(detail,'Runtime',value.source?.runtime==='missing_interpreter'?'Interpreter unavailable':'Not probed by inspection');
        detail.append(e('h4','What goes in'));const inputs=e('ul');
        for(const item of value.inputs||[])inputs.append(e('li',item.label+' · '+(item.required?'required':'optional')+' · '+item.formats.join(', ')));
        detail.append(inputs,e('h4','What comes out'));const outputs=e('ul');for(const text of value.outputs||[])outputs.append(e('li',text));detail.append(outputs);
        field(detail,'Scientific limit',value.claim_ceiling||'Not assessed');
        if(value.synthetic_only)field(detail,'Scope','Shipped invented software example only; no real-protein evaluation');
        const selected=value.input_state||{};field(detail,'Selected inputs',selected.label||'Not selected');
        if(selected.missing?.length)field(detail,'Missing inputs',selected.missing.join('; '));
        field(detail,'Last saved check',selected.recorded_check?.label||'No case selected or no saved check');
        if(selected.recorded_check?.limit)detail.append(e('p',selected.recorded_check.limit));
        field(detail,'Next action',value.next_action);
        if(selected.open_inputs_url){const a=e('a','Open existing input check');a.href=selected.open_inputs_url;a.target='_blank';a.rel='noopener';detail.append(a);}
        const choices=e('div');const choose=button('Choose saved inputs',async()=>{
          const selectedGeneration=++generation;choose.disabled=true;
          try{
            const rows=await get('/api/analyses');if(!dialog.isConnected||selectedGeneration!==generation)return;
            const matching=rows.filter(r=>r.analysis_type===card.id);
            if(!matching.length){choices.replaceChildren(e('p','No saved inputs for this method. No example was created.'));return;}
            const label=e('label','Saved input case'),select=e('select');select.setAttribute('aria-label','Saved input case for '+card.title);select.append(new Option('Choose a case…',''));
            for(const row of matching)select.append(new Option((row.subject_id||row.analysis_id)+' · '+(row.question||row.analysis_id),row.analysis_id));
            select.onchange=()=>{if(select.value)inspect(card,select.value);};label.append(select);choices.replaceChildren(label);
          }catch(error){if(dialog.isConnected&&selectedGeneration===generation)choices.replaceChildren(e('p',error.message));}
          finally{choose.disabled=false;}
        });detail.append(choose,choices);
        const exact=e('details');exact.append(e('summary','Inspection identity and source coverage'),e('p',value.source?.scope||'Not checked'),e('code',value.identity_sha256));detail.append(exact);
        status.textContent='Local inspection complete. No tool probe, validation or job was started.';title.focus();
      }catch(error){if(dialog.isConnected&&serial===generation)status.textContent=error.message;}
    }
    async function load(){
      const serial=++generation;refresh.disabled=true;status.textContent='Reading local method definitions…';cards.replaceChildren();detail.replaceChildren();
      try{
        const data=await get('/api/capabilities');if(!dialog.isConnected||serial!==generation)return;
        for(const card of data.capabilities){
          const article=e('article');article.append(e('h3',card.title));
          field(article,'Connection',card.registration==='registered'?'Registered method':'Unavailable in this distribution');
          field(article,'Local status',names[card.source_state]||'Unknown');
          if(card.identity_sha256){field(article,'Runtime',card.source?.runtime==='missing_interpreter'?'Interpreter unavailable':'Not checked by inspection');if(card.synthetic_only)field(article,'Scope','Invented software check only');article.append(button('Inspect '+card.title,()=>inspect(card)));}
          else article.append(e('p','This distribution does not provide this method. No private adapter or substitute example is loaded.'));
          cards.append(article);
        }
        status.textContent='Local metadata only. Public reference retrieval is '+(data.privacy.public_reference_fetch_enabled?'enabled separately':'disabled')+'. Tool network behavior is not checked here.';
      }catch(error){if(dialog.isConnected&&serial===generation)status.textContent='Capability view unavailable: '+error.message;}
      finally{refresh.disabled=false;}
    }
    await load();
  }
  document.getElementById('capability-view')?.addEventListener('click',open);
})();
