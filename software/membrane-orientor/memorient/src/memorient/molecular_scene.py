"""Bounded heavy-solute display derived from exact final GRO coordinates."""
from pathlib import Path
import json
from .md_preparation import topology_summary


def build_display(run,profile):
    root=Path(run);summary=topology_summary((root/'system.top').read_text())
    rows=(root/'em.gro').read_text().splitlines();offset=0;parts=[];mapping=[];lines=[];serial=0
    # Cartesian coordinates are displayed as saved. No curvature remapping,
    # coordinate centering or guessed molecule unwrapping is applied.
    def flush():
        nonlocal lines,serial
        if not lines:return
        parts.append({'id':'solute-'+str(len(parts)+1),'pdb':'\n'.join(lines+['END'])+'\n',
                      'color':'#82d9c5','view_atoms':serial})
        lines=[];serial=0
    for name,count in summary['sequence']:
        aa=summary['types'][name]
        for instance in range(count):
            for a in aa:
                row=rows[2+offset];offset+=1
                if name in {'SOL','SOD','CLA'} or a['mass']<=2:continue
                if len(mapping)>=100000:raise ValueError('Display exceeds 100000 heavy solute atoms; full coordinates are retained')
                if serial==50000:flush()
                serial+=1;xyz=[float(row[l:r])*10 for l,r in ((20,28),(28,36),(36,44))]
                fields=[f'{v:8.3f}' for v in xyz]
                if any(len(f)!=8 for f in fields):raise ValueError('Display exceeds PDB field limits')
                element='Ca' if name=='CAL' else 'Mg' if name=='MG' else a['name'][0]
                lines.append(f'HETATM{serial:5d} {a["name"][:4]:>4} {a["residue"][:3]:>3} A{a["resid"]%10000:4d}    '+''.join(fields)+f'  1.00  0.00          {element:>2}')
                mapping.append(f'{offset}\tsolute-{len(parts)+1}\t{serial}\t{name}\t{instance+1}\t{a["name"]}')
    flush()
    doc={'schema_version':'memorient.molecular_scene.v1','units':'angstrom','stage':'minimized',
         'organism':profile['organism'],'parts':parts,'full_atoms':summary['atoms'],'view_atoms':len(mapping),
         'box_nm':[float(x) for x in rows[-1].split()],
         'limits':['Heavy solute only; bulk water and NaCl retained in full system.',
                   'Coordinates displayed as saved; periodic molecules can cross box boundaries.',
                   'Restrained minimum only; not equilibrated or biologically validated.']}
    (root/'scene.json').write_text(json.dumps(doc,allow_nan=False)+'\n')
    (root/'ATOM_MAP.tsv').write_text('system_atom_index_1based\tpart\tdisplay_serial\tmolecule_type\tinstance\tatom\n'+'\n'.join(mapping)+'\n')
    return doc
