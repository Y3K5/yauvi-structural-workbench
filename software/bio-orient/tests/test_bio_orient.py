"""Synthetic geometry controls, not biological validation cases."""
import hashlib
import json
from pathlib import Path
import numpy as np
import pytest
from bio_orient.core import analyze, binding, canonical, write_outputs
from structqc.coordinate_scope import load_scope, proper_rotation
from structqc.chirality import check_scope
from memorient.sidedness import resolve_sidedness


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def dump(path, doc):
    path.write_text(json.dumps(doc))
    return path


def pdb_line(serial, name, cid, chain, resid, xyz, element, alt="", het=True):
    record = "HETATM" if het else "ATOM  "
    return f"{record}{serial:5d} {name:^4}{alt or ' '}{cid:>3} {chain}{resid:4d}    {xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}{1:6.2f}{20:6.2f}          {element:>2}\n"


def fixture(tmp_path, cid="ALA", reflect=False, alternate=False):
    atoms = [("CA",(0,0,0),"C"),("N",(-1,0,0),"N"),("C",(0,1,0),"C"),("CB",(0,0,1),"C")]
    text = "".join(pdb_line(i+1,name,cid,"A",1,(x,y,-z if reflect else z),el,"A" if alternate else "",cid!="ALA") for i,(name,(x,y,z),el) in enumerate(atoms))
    if alternate:
        text += pdb_line(8,"CA",cid,"A",1,(0,0,0),"C","B",cid!="ALA")
    text += pdb_line(9,"C1","LIG","B",2,(0,0,3),"C") + "END\n"
    structure = tmp_path / "structure.pdb"
    structure.write_text(text)
    scope = load_scope(structure)
    first = scope["components"][0]
    patches = dump(tmp_path/"patches.json", {"schema_version":"1.0","binding":binding(scope),"sides":[
        {"side_id":"test-side","residue_set":[{"chain_id":first["chain_id"],"auth_seq_id":1,"insertion_code":""}],
         "source":{"id":"synthetic-declaration","citation":"Synthetic software fixture; no biological claims"},"biological_roles":["declared_test_role"]}]})
    cif = f'''data_{cid}
_chem_comp.id {cid}
loop_
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.pdbx_stereo_config
_chem_comp_atom.pdbx_model_Cartn_x_ideal
_chem_comp_atom.pdbx_model_Cartn_y_ideal
_chem_comp_atom.pdbx_model_Cartn_z_ideal
CA C S 0 0 0
N N N -1 0 0
C C N 0 1 0
CB C N 0 0 1
H H N 1 -1 -1
loop_
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
CA N
CA C
CA CB
CA H
'''
    reference = tmp_path / "component.cif"
    reference.write_text(cif)
    refs = dump(tmp_path / "references.json",{"schema_version":"1.0","components":[{"component_id":cid,"path":"component.cif","sha256":digest(reference)}]})
    return structure, patches, refs


@pytest.mark.parametrize("cid",["ALA","DAL","SUG","LIP","LIG"])
@pytest.mark.parametrize("reflect",[False,True])
def test_reference_relative_chirality_for_component_categories(tmp_path,cid,reflect):
    structure,_,refs=fixture(tmp_path,cid,reflect)
    result=check_scope(load_scope(structure),refs)
    assert result["findings"][0]["state"] == ("reference_mismatch" if reflect else "matches_reference")
    assert result["findings"][0]["observed_configuration"] is None
    assert result["component_coverage"][1]["state"] == ("partial_or_unevaluated" if cid == "LIG" else "reference_missing")


def test_alternative_conformers_are_not_a_clean_chirality_result(tmp_path):
    structure,_,refs=fixture(tmp_path,alternate=True)
    f=check_scope(load_scope(structure),refs)["findings"][0]
    assert f["state"]=="unevaluated" and f["reason"]=="alternative_conformers"


def test_missing_atoms_and_references_remain_unknown(tmp_path):
    structure,_,refs=fixture(tmp_path)
    structure.write_text("\n".join(line for line in structure.read_text().splitlines() if " CB " not in line))
    assert check_scope(load_scope(structure),refs)["findings"][0]["reason"]=="missing_required_atoms"
    assert check_scope(load_scope(structure))["state"]=="partial"


@pytest.mark.parametrize("matrix",[np.diag([-1,1,1]),np.diag([2,1,1]),np.full((3,3),np.nan)])
def test_invalid_rigid_transform_is_rejected(matrix):
    with pytest.raises(ValueError):proper_rotation(matrix)


def test_graph_measures_neighbor_burial_without_function_claim(tmp_path):
    s,p,r=fixture(tmp_path)
    doc=analyze(s,p,chemical_reference=r)
    assert doc["edges"] and any(e.get("buried_sasa_A2",0)>0 for e in doc["edges"])
    assert all(e["functional_availability"]["state"]=="unknown" for e in doc["edges"])
    assert all(e["dynamic_probability"] is None for e in doc["edges"])
    assert doc["sides"][0]["availability"]["functional"]=="unknown"
    write_outputs(doc,tmp_path/"out")
    first={p.name:p.read_bytes() for p in (tmp_path/"out").iterdir()}
    write_outputs(analyze(s,p,chemical_reference=r),tmp_path/"out")
    assert first=={p.name:p.read_bytes() for p in (tmp_path/"out").iterdir()}


def test_neighbor_removal_keeps_sample_frame_and_distant_objects_cannot_bury(tmp_path):
    from memorient.sasa import atom_sasa
    coords=np.array([[0.,0.,0.],[0.,0.,3.],[100.,50.,20.]])
    radii=np.array([1.7,1.7,1.7])
    full=atom_sasa(coords,radii)
    retained=atom_sasa(coords,radii,target_indices=[0],excluded_indices=[1])
    direct=atom_sasa(coords[[0,2]],radii[[0,2]])
    assert retained[0]==direct[0] and retained[0]>full[0]
    assert atom_sasa(coords,radii,target_indices=[0],excluded_indices=[2])[0]==full[0]
    s,p,_=fixture(tmp_path)
    text=s.read_text().replace('END\n',pdb_line(20,'C1','FAR','C',3,(100,50,20),'C')+'END\n')
    s.write_text(text)
    declaration=json.loads(p.read_text());declaration['binding']=binding(load_scope(s));dump(p,declaration)
    doc=analyze(s,p)
    edge=next(e for e in doc['edges'] if e['adjacent_object'].endswith('/3/FAR'))
    assert edge['buried_sasa_A2']==0 and edge['contact_atom_pairs']==0


@pytest.mark.parametrize("field,value",[("coordinate_sha256","0"*64),("model_id","2"),("assembly_id","1"),("frame","oriented"),("conformer_policy","mixed")])
def test_wrong_patch_scope_is_rejected(tmp_path,field,value):
    s,p,_=fixture(tmp_path)
    declaration=json.loads(p.read_text());declaration["binding"][field]=value;dump(p,declaration)
    with pytest.raises(ValueError,match="bind"):analyze(s,p)


def test_wrong_chain_copy_is_rejected(tmp_path):
    s,p,_=fixture(tmp_path)
    declaration=json.loads(p.read_text());declaration["sides"][0]["residue_set"][0]["chain_id"]="another-copy";dump(p,declaration)
    with pytest.raises(ValueError,match="chain copy"):analyze(s,p)


def test_missing_patch_coordinates_retains_denominator(tmp_path):
    s,p,_=fixture(tmp_path)
    d=json.loads(p.read_text());row=dict(d["sides"][0]["residue_set"][0]);row["auth_seq_id"]=99;d["sides"][0]["residue_set"].append(row);dump(p,d)
    side=analyze(s,p)["sides"][0]
    assert side["coverage"]["declared_residues"]==2 and side["coverage"]["observed_residues"]==1
    assert side["availability"]["structural"]=="partial_or_missing"


def marker_atoms():
    return [{"atom":"CA","chain_id":"A","auth_seq_id":i,"insertion_code":"","xyz":[0,0,z]} for i,z in [(1,20),(2,-20)]]


@pytest.mark.parametrize("context,compartments",[("gram_negative_om",["extracellular","periplasm"]),("gram_negative_im",["periplasm","cytoplasm"]),("er_membrane",["er_lumen","cytosol"]),("mitochondrial_outer",["cytosol","intermembrane_space"]),("mitochondrial_inner",["intermembrane_space","matrix"]),("organelle_membrane",["custom_lumen","custom_cytosol"])])
def test_context_evidence_assignments(context,compartments):
    t={"compartments":compartments,"source":{"id":"synthetic","citation":"test"},"side_markers":[{"chain_id":"A","auth_seq_id":1,"compartment":compartments[0]}]}
    result=resolve_sidedness(marker_atoms(),np.zeros(3),np.array([0,0,1]),10,context,t)
    assert result["state"]=="supported" and result["positive_compartment"]==compartments[0]


def test_unknown_predicted_and_conflicting_sides():
    atoms=marker_atoms();args=(atoms,np.zeros(3),np.array([0,0,1]),10,"er_membrane")
    assert resolve_sidedness(*args)["state"]=="unknown"
    t={"source":{"id":"synthetic","citation":"test"},"side_markers":[{"chain_id":"A","auth_seq_id":1,"compartment":"er_lumen","basis":"predicted"}]}
    assert resolve_sidedness(*args,t)["state"]=="predicted"
    t["side_markers"].append({"chain_id":"A","auth_seq_id":2,"compartment":"er_lumen"})
    assert resolve_sidedness(*args,t)["state"]=="conflicting"
    t["side_markers"][1]["auth_seq_id"]=999
    with pytest.raises(ValueError,match="exact coordinate"):resolve_sidedness(*args,t)


def test_supporting_state_evidence_keeps_scope_and_no_probability(tmp_path):
    s,p,_=fixture(tmp_path)
    output={"side_records":[{"side_id":"test-side","residue_set":json.loads(p.read_text())["sides"][0]["residue_set"],"dynamic_state":{"label":"open_like","limitations":["Recorded synthetic static comparison"]}}]}
    record={"binding":binding(load_scope(s)),"engine":"state_atlas","method":{"id":"synthetic-state","source_sha256":"a"*64},"output":output,"output_sha256":hashlib.sha256(canonical(output)).hexdigest()}
    e=dump(tmp_path/"state.json",record)
    doc=analyze(s,p,evidence_paths=[e])
    assert doc["sides"][0]["dynamic_state"]["state"]=="recorded_method_output"
    assert all(e["dynamic_probability"] is None for e in doc["edges"])
    record["binding"]["assembly_id"]="wrong";dump(e,record)
    with pytest.raises(ValueError,match="different coordinate"):analyze(s,p,evidence_paths=[e])


def test_supported_membrane_does_not_resolve_a_patch_in_the_core(tmp_path):
    s,p,_=fixture(tmp_path)
    marker={'chain_id':'Axp:identity','auth_seq_id':1,'compartment':'cytosol','basis':'curated'}
    atoms=[{'atom':'CA','chain_id':'Axp:identity','auth_seq_id':1,'xyz':[0,0,0]}]
    assignment=resolve_sidedness(atoms,np.array([0,0,-20]),np.array([0,0,1]),10,'er_membrane',{'source':{'id':'fixture','citation':'Synthetic scope control'},'side_markers':[marker]})
    mem={'frame':'deposited_coordinates','context':'er_membrane','center':[0,0,-20],'normal':[0,0,1],'half_thickness':10,'sidedness':'supported','side_assignment':assignment}
    # The ligand patch is in a separately declared thick slab; its mapped marker
    # remains outside, while the surface itself spans the modeled interface.
    scope=load_scope(s)
    declaration=json.loads(p.read_text())
    declaration['sides'][0]['residue_set']=[{'chain_id':scope['components'][1]['chain_id'],'auth_seq_id':2,'insertion_code':''}]
    dump(p,declaration)
    mem['center']=[0,0,20];mem['normal']=[0,0,-1];mem['half_thickness']=18
    mem['side_assignment']=resolve_sidedness(atoms,np.array(mem['center']),np.array(mem['normal']),18,'er_membrane',{'source':{'id':'fixture','citation':'Synthetic scope control'},'side_markers':[marker]})
    output={'input_coordinate_membrane':mem}
    evidence=dump(tmp_path/'membrane.json',{'binding':binding(scope),'engine':'memorient','method':{'id':'fixture','source_sha256':'a'*64},'output':output,'output_sha256':hashlib.sha256(canonical(output)).hexdigest()})
    result=analyze(s,p,evidence_paths=[evidence])
    assert result['membrane']['sidedness']=='supported'
    assert result['sides'][0]['compartment_orientation']['compartment'] is None
    assert result['sides'][0]['availability']['orientational']=='unresolved'


def test_rigid_rotation_preserves_chirality_contacts_and_exposure(tmp_path):
    s,p,r=fixture(tmp_path)
    first=analyze(s,p,chemical_reference=r)
    scope=load_scope(s)
    rotation=np.array([[0,-1,0],[1,0,0],[0,0,1]])
    text="".join(pdb_line(i+1,a["atom"],a["resn"],"A" if a["resi"]==1 else "B",a["resi"],rotation@np.array([a[k] for k in "xyz"])+[9,4,2],a["elem"],het=a["hetflag"]) for i,a in enumerate(scope["atoms"]))+"END\n"
    s.write_text(text);decl=json.loads(p.read_text());decl["binding"]=binding(load_scope(s));dump(p,decl)
    second=analyze(s,p,chemical_reference=r)
    assert first["chirality"]["findings"][0]["state"]==second["chirality"]["findings"][0]["state"]
    assert first["edges"][0]["contact_atom_pairs"]==second["edges"][0]["contact_atom_pairs"]
    assert first["sides"][0]["solvent_exposure"]["sasa_A2"]==pytest.approx(second["sides"][0]["solvent_exposure"]["sasa_A2"])


def test_wrong_ccd_checksum_is_rejected(tmp_path):
    s,p,r=fixture(tmp_path)
    (tmp_path/"component.cif").write_text((tmp_path/"component.cif").read_text()+"\n# changed\n")
    with pytest.raises(ValueError,match="checksum"):analyze(s,p,chemical_reference=r)


def test_missing_stereo_column_cannot_be_clean(tmp_path):
    s,p,r=fixture(tmp_path)
    c=tmp_path/"component.cif"
    text=c.read_text().replace("_chem_comp_atom.pdbx_stereo_config\n", "")
    for element in ["C", "N", "H"]:
        text=text.replace(f" {element} S ", f" {element} ").replace(f" {element} N ", f" {element} ")
    c.write_text(text)
    refs=json.loads(r.read_text());refs["components"][0]["sha256"]=digest(c);dump(r,refs)
    assert analyze(s,p,chemical_reference=r)["chirality"]["state"]=="partial"


def test_unsupported_atom_element_cannot_be_clean(tmp_path):
    s,p,r=fixture(tmp_path)
    lines=s.read_text().splitlines()
    lines[0]=lines[0][:76]+'Fe'
    s.write_text('\n'.join(lines)+'\n')
    assert check_scope(load_scope(s),r)["findings"][0]["reason"]=="element_identity_mismatch"


def test_declared_d_residue_uses_its_own_reference(tmp_path):
    s,p,r=fixture(tmp_path,"DAL",reflect=True)
    c=tmp_path/"component.cif";c.write_text(c.read_text().replace("CB C N 0 0 1","CB C N 0 0 -1"))
    refs=json.loads(r.read_text());refs["components"][0]["sha256"]=digest(c);dump(r,refs)
    finding=analyze(s,p,chemical_reference=r)["chirality"]["findings"][0]
    assert finding["chemical_component_id"]=="DAL" and finding["state"]=="matches_reference"


def envelope(path,structure,output,engine="site_context"):
    return dump(path,{"binding":binding(load_scope(structure)),"engine":engine,"method":{"id":"synthetic-method","source_sha256":"a"*64},"output":output,"output_sha256":hashlib.sha256(canonical(output)).hexdigest()})


def test_same_side_name_with_wrong_residues_cannot_receive_state(tmp_path):
    s,p,_=fixture(tmp_path)
    e=envelope(tmp_path/"evidence.json",s,{"side_records":[{"side_id":"test-side","residue_set":[],"dynamic_state":{"label":"open"}}]},"state_atlas")
    with pytest.raises(ValueError,match="exact declared residue set"):analyze(s,p,evidence_paths=[e])


def test_only_annotated_edge_receives_functional_evidence(tmp_path):
    s,p,_=fixture(tmp_path)
    original=analyze(s,p)
    edge=original["edges"][0]
    item={"source_side":edge["source_side"],"adjacent_object":edge["adjacent_object"],"residue_set":original["sides"][0]["residue_set"],
          "functional_availability":{"state":"predicted","value":"occluded"},"source":{"id":"synthetic","citation":"Synthetic annotation only"},"interpretation_limit":"No measured affinity or native-accessibility claim"}
    e=envelope(tmp_path/"evidence.json",s,{"edge_records":[item]})
    graph=analyze(s,p,evidence_paths=[e])
    assert graph["edges"][0]["functional_availability"]["state"]=="recorded_evidence"
    assert all(x["functional_availability"]["state"]=="unknown" for x in graph["edges"][1:])
    item["residue_set"]=[];envelope(e,s,{"edge_records":[item]})
    with pytest.raises(ValueError,match="residue scope"):analyze(s,p,evidence_paths=[e])


def test_sidedness_terms_have_independent_scopes(tmp_path):
    s,p,_=fixture(tmp_path)
    declaration=json.loads(p.read_text());a=declaration["sides"][0];a["equivalence_group"]="same_family"
    b={**a,"side_id":"second-side","biological_roles":["different_role"]}
    declaration["sides"].append(b);dump(p,declaration)
    descriptors=analyze(s,p)["sidedness_descriptors"]
    assert descriptors["multisided"] and descriptors["homosided"] and descriptors["heterosided"]


def test_workbench_reports_export_and_review_flags_are_nonblocking(tmp_path):
    from yauvi_platform.structural_workbench import StructuralAnalysisStore
    s,p,r=fixture(tmp_path,reflect=True)
    provenance=dump(tmp_path/"provenance.json",{"class":"predicted","method":"Synthetic software control"})
    validation=dump(tmp_path/"validation.json",{"method":"Synthetic validator fixture, no biological validation","clashscore":0})
    store=StructuralAnalysisStore(tmp_path)
    store.create("bio-case",analysis_type="bio_orient",question="Inspect synthetic adjacency",subject_id="Synthetic control")
    for role,path in [("structure",s),("patch_declaration",p),("chemical_reference",r),("chemical_component",tmp_path/"component.cif"),("provenance",provenance),("validation_report",validation)]:
        store.add_file("bio-case",role=role,path=path)
    assert store.preflight("bio-case")["valid"]
    run=store.run("bio-case")
    assert run["status"]=="completed",run
    snapshot=store.snapshot("bio-case")
    graph=next(d["document"] for d in snapshot["report"]["documents"] if d["path"].endswith("/BIO_ORIENT.json"))
    assert graph["chirality"]["state"]=="review_required"
    import zipfile
    with zipfile.ZipFile(store.artifact_path("bio-case",run["run_id"],"RAW_EVIDENCE.zip")) as archive:
        assert "outputs/bio_orient/BIO_ORIENT.json" in archive.namelist()
        assert "outputs/bio_orient/BIO_ORIENT_EDGES.tsv" in archive.namelist()
    assert store.run("bio-case")["run_id"]==run["run_id"]
