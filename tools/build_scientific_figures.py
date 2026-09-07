from pathlib import Path
import os
import tempfile
os.environ.setdefault('MPLCONFIGDIR', str(Path(tempfile.gettempdir()) / 'yauvi-figure-cache'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from pathlib import Path
import json,hashlib,csv,shutil
R=Path(__file__).resolve().parents[1];P=R/'yauvi-structural-workbench/paper';F=P/'figures';F.mkdir(exist_ok=True)
S=R/'yauvi-structural-workbench/implementation-evidence/2026-09-06'
qp=S
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'svg.fonttype':'none'})
navy='#193641';teal='#26786b';gold='#a16929'
def save(fig,name):
 fig.savefig(F/(name+'.png'),dpi=200,bbox_inches='tight',facecolor='white');fig.savefig(F/(name+'.svg'),bbox_inches='tight',facecolor='white');plt.close(fig)
 svg=F/(name+'.svg')
 svg.write_text('\n'.join(line.rstrip() for line in svg.read_text().splitlines())+'\n')
fig,ax=plt.subplots(figsize=(11,3.6));ax.axis('off');ax.set_xlim(0,11);ax.set_ylim(0,3.6)
ax.text(0,3.3,'A shared path from input to interpretation',fontsize=18,fontweight='bold',color=navy)
labels=[('Input','Structure + source\ncontent hashes'),('Mapping','Model, assembly,\nchain, residue'),('Analysis','Shared engine\nCLI or browser'),('Evidence','Measurements, gaps,\nimmutable run'),('Interpretation','Bounded claims\nportable report')]
for i,(title,body) in enumerate(labels):
 x=i*2.2;ax.add_patch(FancyBboxPatch((x,.9),1.9,1.6,boxstyle='round,pad=.08',facecolor='#edf4f2',edgecolor=teal));ax.text(x+.95,2.05,title,ha='center',fontweight='bold',color=navy);ax.text(x+.95,1.45,body,ha='center',va='center',fontsize=10)
 if i<4:ax.annotate('',xy=(x+2.08,1.7),xytext=(x+1.95,1.7),arrowprops={'arrowstyle':'->','color':navy})
ax.text(0,.3,'Execution outcome • Scientific qualification • Independent reproduction • Publication approval',fontsize=11,color=navy)
ax.text(0,-.1,'Four separate states; success in one does not establish the others.',fontsize=10,color=gold)
save(fig,'architecture')
e=json.loads((qp/'installed-example-STRUCTURE_EVIDENCE.json').read_text());rows=list(csv.DictReader((qp/'installed-example-RESIDUE_QUALITY.tsv').open(),delimiter='\t'))
fig,(ax,side)=plt.subplots(1,2,figsize=(10,3.6),gridspec_kw={'width_ratios':[1,1.5]});xyz=[[float(c) for c in r['ca_xyz'].split(';')] for r in rows];ax.plot([c[0] for c in xyz],[c[1] for c in xyz],color=teal,lw=2)
for r,c in zip(rows,xyz):ax.scatter(*c[:2],s=120,color=teal);ax.annotate(f"{r['chain_id']}:{r['auth_seq_id']} {r['chemical_component_id']}",c[:2],xytext=(8,-18),textcoords='offset points')
ax.set_xlabel('x (Å)');ax.set_ylabel('y (Å)');ax.set_title('Cα coordinates · synthetic example');ax.set_xlim(-.6,4);ax.set_ylim(-.8,2.4);ax.spines[['top','right']].set_visible(False)
side.axis('off');side.text(0,.94,'What was actually measured?',fontsize=16,fontweight='bold',color=navy)
lines=[f"Mapped residues: {e['completeness']['mapped_residues']} of {e['completeness']['reference_length']}",f"Sequence identity: {e['completeness']['identity_fraction']:.0%} within this toy mapping",f"Imported clashscore: {e['external_validation']['metrics']['clashscore']} (synthetic fixture)",'PAE: not supplied; not inferred','Input and validation records identified by SHA-256','Limit: these two residues demonstrate software behavior.','They do not establish structure validity or function.']
for i,t in enumerate(lines):side.text(0,.78-i*.105,t,color=navy if i<5 else gold,fontsize=10.5)
save(fig,'worked-example')
summary_path=R/'yauvi-structural-workbench/benchmarks/qualification-v2/results/EXECUTION_SUMMARY.json';summary=json.loads(summary_path.read_text());by={p['workflow']:p for p in summary['panels']}
fig,ax=plt.subplots(figsize=(10.8,3.8));ax.axis('off');ax.set_title('Historical collection 2.9 · changed code requires fresh qualification',loc='left',fontsize=15,color=navy,pad=20)
pretty=[('structure_qc','StructQC'),('functional_site_state','SiteContext'),('assembly_interface','AssemblyContext'),('conformational_state','ABL StateAtlas'),('sf_csa','SF-CSA'),('membrane_orientation','MembraneOrient')]
rs=[]
for key,label in pretty:
 p=by.get(key);rs.append([label,'14: 4 refs + 10 held-out' if key=='conformational_state' else '16', 'Unexecuted' if not p else f"{p['cases']['passed']}/{p['cases']['total']}", 'Pending source recovery' if key=='sf_csa' else 'Experimental' if key=='membrane_orientation' else 'Curated case execution'])
t=ax.table(cellText=rs,colLabels=['Workflow','Cases in compared panel','Historical result','Interpretation'],loc='center',cellLoc='left',colLoc='left',colWidths=[.2,.26,.18,.36]);t.auto_set_font_size(False);t.set_fontsize(10.5);t.scale(1,1.7)
for (r,c),cell in t.get_celld().items():cell.set_edgecolor('white');cell.set_facecolor(navy if r==0 else '#edf4f2' if r%2 else '#f7f9f8');cell.set_text_props(color='white' if r==0 else navy)
ax.text(0,-.07,'Six archived runners: four required panels agree; membrane is 5/16 in five and 4/16 in one.\nControls are separate (four per runner). Unexecuted alpha-helical requirements are not shown as tested cases.\nSmall curated panels do not establish broad biological validity; no combined accuracy score.',transform=ax.transAxes,fontsize=9,color=gold)
save(fig,'qualification-history')
record={'schema_version':'1.0','figure_source_sha256':{'historical_summary':hashlib.sha256(summary_path.read_bytes()).hexdigest(),'installed_example':hashlib.sha256((qp/'installed-example-STRUCTURE_EVIDENCE.json').read_bytes()).hexdigest(),'installed_residue_table':hashlib.sha256((qp/'installed-example-RESIDUE_QUALITY.tsv').read_bytes()).hexdigest()},'figures':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(F.iterdir()) if p.suffix in {'.png','.svg'}},'limits':'Worked example synthetic; qualification figure historical collection 2.9; new code not scientifically qualified.'}
(F/'FIGURE_PROVENANCE.json').write_text(json.dumps(record,indent=2)+'\n')
print('Regenerated three figures from versioned evidence.')
