from pathlib import Path
from structqc.core import analyze

def test_linked_modified_residue_preserves_chemical_identity(tmp_path):
    rows=[]
    # Explicit peptide C-N connections; final PTR is isolated and must stay out.
    for i,(record,name,num,offset) in enumerate([('ATOM','ALA',1,0),('HETATM','PTR',2,3.9),('ATOM','GLY',3,7.8),('HETATM','PTR',99,40)]):
        for j,(atom,x,element) in enumerate([('N',0,'N'),('CA',1.3,'C'),('C',2.6,'C'),('O',3.0,'O')]):
            rows.append(f'{record:<6}{i*4+j+1:5d} {atom:^4} {name:3} A{num:4d}    {offset+x:8.3f}{0:8.3f}{0:8.3f}{1:6.2f}{20:6.2f}          {element:>2}')
    p=tmp_path/'modified.pdb';p.write_text('\n'.join(rows)+'\nEND\n')
    document=analyze(p,reference_sequence='AYG')
    records=document['residues']
    assert len(records)==3
    modified=next(r for r in records if r['resname']=='PTR')
    assert modified['chemical_component_id']=='PTR'
    assert modified['normalized_parent_component_id']=='TYR'
    assert modified['one_letter']=='Y'
    assert modified['sequence_normalization']=='explicit_modified_parent_v1'
    assert modified['auth_seq_id']==2 and modified['model_index']==0
