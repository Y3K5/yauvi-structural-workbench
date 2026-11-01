"""Source-derived protein records and native unknown-function investigation plans.
No function classifier, scientific executor, ontology assignment or second store.
"""
from __future__ import annotations
from copy import deepcopy
from functools import lru_cache
import hashlib
import html
import json
import os
from .function_catalogue import catalogue, bounded, SHA
from .store import AnalysisError, _canonical

SCHEMA = 'protein-function-record-v1'
INVESTIGATION = 'protein-function-investigation-v1'
SOURCE_SCHEMAS = {'refinement':'protein-function-refinement-v1', 'duf':'local-duf-function-review-v1', 'bam':'bam-domain-function-review-v1'}
GAPS = {'all':'All protein records','function':'Original function unresolved','specificity':'Family / domain specificity needs review','location':'Location unplaced / conflicting','partners':'No partner evidence joined','structure':'No qualified saved form'}
ACTIONS = [
 {'id':'identity','title':'Resolve exact reference and material','question':'Do accession, full sequence, source strain, processed form and numbering refer to the same material?','needed':'Pinned source records, exact sequence hashes and explicit strain/region correspondence.','limit':'Identity correspondence alone does not establish function.'},
 {'id':'domains','title':'Qualify family and domain clues','question':'Which versioned family/domain evidence covers this exact sequence and region?','needed':'Original domain/family calls, boundaries, coverage, database version and method parameters.','limit':'A domain or family label does not establish a specific native substrate or activity.'},
 {'id':'homology','title':'Compare functional explanations','question':'Do homologous references share the relevant domain architecture and the evidence needed for the proposed role?','needed':'Exact reference/material identities, alignment scope and original functional observations.','limit':'Homology and multiple copied annotations are not independent target-function experiments.'},
 {'id':'location','title':'Resolve cellular context','question':'Is there evidence for the relevant compartment, exposure and biological condition?','needed':'Separate localization predictions, source observations, material identity and condition.','limit':'A predicted compartment or schematic placement does not establish native accessibility.'},
 {'id':'partners','title':'Qualify partners and complex membership','question':'Is the relationship direct binding, complex membership, genetic dependency or phenotype association?','needed':'Typed source relationship, exact materials, controls and competing explanations.','limit':'Complex membership, coexpression and cell phenotype do not establish a reciprocal atomic contact.'},
 {'id':'function','title':'Define a distinguishing functional observation','question':'What observation would distinguish the proposed role from its strongest alternative?','needed':'A precise question, endpoint, comparison, material identity and interpretation limits.','limit':'This is a study-planning question; feasibility and method qualification remain open.'},
 {'id':'condition','title':'Identify condition and expression gaps','question':'Under what cellular/environmental condition is this protein present and relevant?','needed':'Source-qualified expression, abundance or phenotype context and experimental unit.','limit':'Reference proteome membership does not establish expression, activity or essentiality.'}
]

def digest(value): return hashlib.sha256(_canonical(value)).hexdigest()
def pointer(value): return str(value).replace('~','~0').replace('/','~1')

@lru_cache(maxsize=8)
def parsed(raw): return json.loads(raw)


def _supplements(data, meta):
    name=os.environ.get('YAUVI_PROTEIN_FUNCTION_EVIDENCE')
    if not name:return {}, {'state':'not_configured','sources':{},'limit':'No optional local review sources joined; base catalogue records remain available.'}
    try:
        raw=bounded(name,16384);config=parsed(raw)
        if config.get('schema')!='protein-function-evidence-config-v1' or config.get('dataset_id')!=data['dataset_id'] or config.get('catalogue_sha256')!=meta['catalogue_sha256']:raise AnalysisError('Optional function evidence is bound to a different catalogue')
        sources=config.get('sources')
        if not isinstance(sources,dict) or not sources or set(sources)-set(SOURCE_SCHEMAS):raise AnalysisError('Unknown optional function-evidence source')
        by={r['key']:r for r in data['records']};docs={};pins={}
        for kind,spec in sources.items():
            if not isinstance(spec,dict) or set(spec)!={'path','sha256'} or not SHA.fullmatch(str(spec['sha256'])):raise AnalysisError('Invalid optional source identity')
            contents=bounded(spec['path'])
            if hashlib.sha256(contents).hexdigest()!=spec['sha256']:raise AnalysisError('Optional function evidence changed: '+kind)
            doc=parsed(contents)
            if doc.get('schema')!=SOURCE_SCHEMAS[kind] or doc.get('base_dataset_id')!=data['dataset_id'] or not isinstance(doc.get('records'),dict):raise AnalysisError('Optional source contract differs: '+kind)
            for key,row in doc['records'].items():
                if key not in by or not isinstance(row,dict) or row.get('key')!=key or row.get('sequence_sha256')!=by[key]['sequence_sha256'] or ('length' in row and row['length']!=by[key]['length']) or any(field in row and row[field]!=by[key][field] for field in ['id','species']):raise AnalysisError('Optional source has an unbound protein identity: '+kind)
            docs[kind]=doc;pins[kind]={'sha256':spec['sha256'],'schema':doc['schema'],'records':len(doc['records'])}
        return docs,{'state':'current','configuration_sha256':hashlib.sha256(raw).hexdigest(),'sources':pins,'limit':'Retained reviews only; original inferred/unknown dispositions are unchanged.'}
    except (OSError,ValueError,KeyError,TypeError,AnalysisError) as exc:
        return {},{'state':'unavailable','sources':{},'limit':'Optional evidence could not be qualified. Base catalogue remains readable.','reason':str(exc)}


def bundle():
    data,meta=catalogue();docs,supp=_supplements(data,meta)
    return data,meta,docs,supp


def relationships(data,profile):
    """Join the complete retained graph, enriching names without assigning identities."""
    atlas=profile.get('atlas_id');joined=[];seen=set()
    if not atlas and not profile.get('connections'):return []
    names={p.get('atlas_id'):p.get('alias') for p in data.get('profiles',{}).values() if p and p.get('atlas_id')}
    for index,edge in enumerate(data.get('source_relationships',[])):
        if atlas and atlas in {edge['a'],edge['b']}:
            other=edge['b'] if edge['a']==atlas else edge['a']
            local=next((x for x in profile.get('connections',[]) if x.get('id')==edge.get('id')),{})
            name=local.get('partner_name') or names.get(other) or other
            joined.append((edge,'/source_relationships/'+str(index),name,other))
            seen.add((edge.get('id'),edge.get('claim'),edge.get('limit'),edge.get('level')))
    for index,edge in enumerate(profile.get('connections',[])):
        signature=(edge.get('id'),edge.get('claim'),edge.get('limit'),edge.get('level'))
        if signature in seen:continue
        joined.append((edge,None,edge.get('partner_name','Source counterpart unresolved'),None));seen.add(signature)
    return joined


def derive(data,meta,docs,supp,row,index):
    key=row['key'];profile=data.get('profiles',{}).get(key) or {};evidence=data['evidence'][key]
    reviews={kind:deepcopy(doc['records'][key]) for kind,doc in docs.items() if key in doc['records']}
    refinement=reviews.get('refinement',{});duf=reviews.get('duf',{});bam=reviews.get('bam',{})
    sources={'catalogue':{'sha256':meta['catalogue_sha256'],'schema':data['schema']}}
    sources.update({kind:spec for kind,spec in supp['sources'].items() if kind in reviews})
    assertions=[]
    def add(predicate,value,kind,locator,*,source='catalogue',limit='',scope=None):
        a={'subject_key':key,'predicate':predicate,'value':deepcopy(value),'source_kind':kind,'source':{'id':source,'sha256':sources[source]['sha256'],'json_pointer':locator},'review_state':'retained_source_unmodified','polarity':'source_assertion','scope':scope or 'Retained reference annotation; exact experimental material/condition not established by this assertion.','limit':limit or 'Source statement only; not newly verified native function.'}
        a['id']='claim-'+digest(a)[:24];assertions.append(a)
    add('annotated_as',row['name'],'retained_annotation',f'/records/{index}/name')
    if profile.get('role'):add('reported_role',profile['role'],'retained_source_context','/profiles/'+pointer(key)+'/role',limit=profile.get('caveat',''),scope=profile.get('join'))
    if profile.get('reviewed_role'):add('reviewed_role',profile['reviewed_role'],'homology_inference','/profiles/'+pointer(key)+'/reviewed_role',limit=profile.get('caveat',''))
    if evidence['location'].get('prediction'):add('location_prediction',evidence['location']['prediction'],'computational_prediction','/evidence/'+pointer(key)+'/location/prediction',limit='Retained predictions/votes, not probabilities or measured native placement.')
    if evidence['location'].get('context'):add('location_context',evidence['location']['context'],'retained_source_context','/evidence/'+pointer(key)+'/location/context',limit=evidence['location'].get('reason',''),scope=profile.get('join'))
    if refinement:add('annotation_refinement',{'detail_group':refinement['detail_group'],'detail_level':refinement['detail_level'],'assignments':refinement.get('assignments',[])},'annotation_organization','/records/'+pointer(key),source='refinement',limit=refinement.get('limit',''))
    if duf:add('duf_review',{'assessment':duf['assessment'],'proposed_role':duf['proposed_role'],'channels':duf.get('channels',{})},'retained_computational_review','/records/'+pointer(key),source='duf',limit=duf.get('disposition','No measured target function established.'))
    if bam:add('bam_reference_role',{'group':bam.get('group_label'),'mechanism':bam.get('mechanism'),'assignment_state':bam.get('assignment_state'),'measured_w83_function':bam.get('measured_w83_function')},'homology_inference','/records/'+pointer(key),source='bam',limit=bam.get('hold',''))
    partners=[]
    for n,(edge,locator,partner_name,partner_entity) in enumerate(relationships(data,profile)):
        if locator is None:locator='/profiles/'+pointer(key)+'/connections/'+str(profile.get('connections',[]).index(edge))
        add('reported_relationship',edge,'retained_source_relationship',locator,limit=edge.get('limit',''),scope=profile.get('study_context') or profile.get('join'))
        partners.append({'name':partner_name,'source_entity_id':partner_entity,'species':edge.get('partner_species'),'relation':edge.get('relation'),'original_evidence_level':edge.get('level'),'claim':edge.get('claim'),'limit':edge.get('limit'),'source_ids':edge.get('source_ids',[]),'assertion_id':assertions[-1]['id'],'selected_reference_transfer':'Not established by this relationship alone; read the studied strain/material context.'})
    original_unresolved=row['primary_group']=='unresolved'
    has_hypothesis=bool(bam or profile.get('reviewed_role') or duf.get('assessment') in {'cached_function_inference','sequence_role_hypothesis','structural_hypothesis'})
    function_state='reference_supported_hypothesis' if has_hypothesis else 'unresolved' if original_unresolved else 'annotation_only'
    level=refinement.get('detail_level')
    specificity=level in {'family_hint','architecture_only','location_only','unresolved'} or original_unresolved
    unknowns=[{'dimension':'function','state':'not_qualified','question':'Which exact activity or causal cellular role is supported for this reference and condition?','action_id':'function'},
              {'dimension':'condition','state':'not_joined','question':'Expression, activity and essentiality in a defined biological condition are not joined by this adapter.','action_id':'condition'}]
    if specificity:unknowns.append({'dimension':'specificity','state':'unresolved','question':'Distinguish the recorded family/domain/architecture from a specific functional role.','action_id':'domains'})
    if evidence['location']['status'] in {'unresolved','conflict'} or evidence['location'].get('zone')=='unplaced':unknowns.append({'dimension':'location','state':evidence['location']['status'],'question':'Qualify cellular location and reconcile any conflicting predictions/source observations.','action_id':'location'})
    else:unknowns.append({'dimension':'location','state':'scope_limited','question':'The recorded location context does not establish native coordinates, accessibility or orientation.','action_id':'location'})
    if not partners:unknowns.append({'dimension':'partners','state':'not_joined','question':'No typed partner relationship joined here; this is missing coverage, not absence of interaction.','action_id':'partners'})
    if not evidence['structure'].get('available'):unknowns.append({'dimension':'structure','state':'not_joined','question':'No qualified saved molecular form is joined; do not infer a fold or native pose.','action_id':'domains'})
    for text in [*duf.get('remaining_unknowns',[]),*(profile.get('study_context') or {}).get('connections',[])]:
        if isinstance(text,str):unknowns.append({'dimension':'source_review','state':'unresolved','question':text,'action_id':'homology'})
        elif isinstance(text,dict):
            for question in text.get('unknowns',[]):unknowns.append({'dimension':'source_review','state':'unresolved','question':question,'action_id':'identity'})
    if row.get('holds') or row.get('quality_flags'):unknowns.append({'dimension':'recorded_holds','state':'held','question':'Original review holds and quality flags remain active; inspect them before transferring claims.','action_id':'identity'})
    if supp['state']!='current':unknowns.append({'dimension':'review_coverage','state':supp['state'],'question':supp.get('reason',supp['limit']),'action_id':'identity'})
    profile_sources=profile.get('sources',[]);ids={x.get('id') if isinstance(x,dict) else x for x in profile_sources}
    for p in partners:ids.update(p['source_ids'])
    method_context={}
    if refinement:method_context['refinement']={k:deepcopy(docs['refinement'][k]) for k in ['boundary','levels','method_sources'] if k in docs['refinement']}
    if duf:method_context['duf']={k:deepcopy(docs['duf'][k]) for k in ['policy','tools','execution'] if k in docs['duf']}
    if bam:method_context['bam']={k:deepcopy(docs['bam'][k]) for k in ['methods','evidence'] if k in docs['bam']}
    result={'schema':SCHEMA,'derivation_version':1,'identity':{k:row[k] for k in ['key','id','species','strain','length','sequence_sha256']},'name':row.get('alias') or row['name'],'binding':{'dataset_id':data['dataset_id'],'catalogue_sha256':meta['catalogue_sha256'],'supplement_state':deepcopy(supp)},'original_record':deepcopy(row),'function':{'state':function_state,'annotation':row['name'],'reported_role':profile.get('role'),'detail_level':level,'hypothesis':bam.get('mechanism') or duf.get('proposed_role'),'original_unresolved':original_unresolved,'terms_state':'No new GO/EC/pathway terms assigned by this adapter.'},'location':deepcopy(evidence['location']),'process':{'state':'recorded_groups_only','groups':row.get('groups',[row['primary_group']]),'specific_pathway':'not_joined'},'condition':{'state':'not_joined','limit':'Reference membership does not establish condition-specific expression or activity.'},'partners':partners,'complex_context':profile.get('same_organism_summary'),'structure':deepcopy(evidence['structure']),'study_context':deepcopy(profile.get('study_context')),'reference_transfer':profile.get('join'),'caveat':profile.get('caveat'),'assertions':assertions,'sources':sources,'bibliography':[deepcopy(x) for x in data.get('source_references',[]) if x['id'] in ids],'profile_bibliography':deepcopy(profile_sources),'retained_reviews':reviews,'review_method_context':method_context,'evidence_inventory':deepcopy(refinement.get('inventory',[])),'unknowns':unknowns,'investigation_actions':deepcopy(ACTIONS),'evidence_base_url':meta['evidence_base_url'],'scope':'Versioned presentation of retained claims and explicit gaps. No new function assignment, ontology inference, measured target function or scientific execution.'}
    result['record_sha256']=digest(result)
    return result


def record(key):
    data,meta,docs,supp=bundle()
    for index,row in enumerate(data['records']):
        if row['key']==key:return derive(data,meta,docs,supp,row,index)
    raise AnalysisError('Protein is not in this exact function catalogue')


def browse(*,species='',group='',q='',gap='all',page='0'):
    if gap not in GAPS or not isinstance(q,str) or len(q)>200 or not isinstance(page,str) or not page.isdecimal() or len(page)>6:raise AnalysisError('Invalid function-record filter')
    data,meta,docs,supp=bundle()
    if species and species not in {x['id'] for x in data['species']}:raise AnalysisError('Unknown reference microbe')
    if group and group not in {x['id'] for x in data['groups']}:raise AnalysisError('Unknown original group')
    summaries=[];counts={k:0 for k in GAPS};query=q.strip().casefold()
    for r in data['records']:
        p=data.get('profiles',{}).get(r['key']) or {};e=data['evidence'][r['key']];ref=docs.get('refinement',{}).get('records',{}).get(r['key'],{});flags={'all':True,'function':r['primary_group']=='unresolved','specificity':ref.get('detail_level') in {'family_hint','architecture_only','location_only','unresolved'} or r['primary_group']=='unresolved','location':e['location'].get('zone')=='unplaced' or e['location']['status']=='conflict','partners':not bool(relationships(data,p)),'structure':not bool(e['structure'].get('available'))}
        if species and r['species']!=species or group and r['primary_group']!=group:continue
        if query and query not in ' '.join(str(r.get(k,'')) for k in ['id','name','alias','strain']).casefold():continue
        for k,yes in flags.items():counts[k]+=int(yes)
        if flags[gap]:summaries.append({'key':r['key'],'id':r['id'],'name':r.get('alias') or r['name'],'strain':r['strain'],'species':r['species'],'primary_group':r['primary_group'],'annotation':r['name'],'original_function_unresolved':flags['function'],'location_label':e['location'].get('label',e['location']['status']),'partner_records':len(relationships(data,p)),'saved_form':bool(e['structure'].get('available'))})
    number=int(page)
    if number and number*25>=len(summaries):raise AnalysisError('Page outside the filtered function records')
    return {'schema':'protein-function-catalogue-v1','total':len(data['records']),'matching':len(summaries),'page':number,'page_size':25,'species':data['species'],'groups':[{'id':x['id'],'label':x['label']} for x in data['groups']],'gap_filters':[{'id':k,'label':v,'count':counts[k]} for k,v in GAPS.items()],'records':summaries[number*25:number*25+25],'binding':{'dataset_id':data['dataset_id'],'catalogue_sha256':meta['catalogue_sha256']},'supplement_state':supp,'limit':'Gap counts describe joined coverage and original review states; they do not measure biological absence.'}


def report_section(study,packet):
    esc=lambda x:html.escape(str(x));r=packet['function_record'];plan=packet['plan']
    assertions=''.join('<li>'+esc(a['predicate'])+': '+esc(json.dumps(a['value'],ensure_ascii=False))+'<br>Source '+esc(a['source']['id'])+' / '+esc(a['source']['sha256'])+' / '+esc(a['source']['json_pointer'])+'<br>'+esc(a['limit'])+'</li>' for a in r['assertions'])
    hypotheses=''.join('<li>'+esc(h['statement'])+'<br>Alternative: '+esc(h['alternative'])+'<br>Distinguishing evidence: '+esc(h['distinguishing_evidence'])+'<br>Support: '+esc(', '.join(h['assertion_ids']) or 'No support selected; unbacked proposal')+'</li>' for h in plan['hypotheses'])
    gaps=''.join('<li>'+esc(x['dimension'])+': '+esc(x['question'])+'</li>' for x in r['unknowns'])
    return '<section><h2>Protein function record and unknown-function investigation</h2><p>'+esc(r['name'])+' · '+esc(r['identity']['id'])+' · '+esc(r['identity']['strain'])+'</p><p>Exact sequence '+esc(r['identity']['sequence_sha256'])+'; record '+esc(r['record_sha256'])+'</p><p>Retained source assertions and proposed hypotheses only. Original annotations/holds remain unchanged.</p><h3>Source-qualified assertions</h3><ul>'+assertions+'</ul><h3>Unresolved questions</h3><ul>'+gaps+'</ul><h3>Proposed investigation</h3><p>'+esc(plan['question'])+'</p><ul>'+hypotheses+'</ul><p>'+esc(plan['notes'])+'</p><pre>'+esc(json.dumps(plan.get('origin_study'),indent=2))+'</pre></section>'
