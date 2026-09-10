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
fig,ax=plt.subplots(figsize=(11,4.6));ax.axis('off');ax.set_xlim(0,11);ax.set_ylim(-.95,4.6)
ax.text(0,4.28,'Six structural evidence questions',fontsize=18,fontweight='bold',color=navy)
ax.text(0,3.98,'One gated manifest, five independent analyses',fontsize=10.5,color=navy)
ax.add_patch(FancyBboxPatch((0,.58),2.4,3.27,boxstyle='round,pad=.06',facecolor='#edf4f2',edgecolor=teal))
for y,t,kw in [(3.42,'Coordinate trust',{'fontweight':'bold','fontsize':12}),(3.14,'can I trust this file?',{'fontsize':10}),(2.86,'structqc',{'fontsize':10,'color':teal,'family':'DejaVu Sans Mono'}),(2.36,'identity · provenance',{'fontsize':9.5}),(2.14,'SHA-256 · residue map',{'fontsize':9.5}),(1.72,'provenance is recorded,',{'fontsize':9,'style':'italic','color':gold}),(1.52,'never inferred',{'fontsize':9,'style':'italic','color':gold}),(1.08,'no downstream module may',{'fontsize':8.5}),(.90,'interpret unbound coordinates',{'fontsize':8.5})]:
 ax.text(1.2,y,t,ha='center',color=kw.pop('color',navy),**kw)
ax.plot([.35,2.05],[2.62,2.62],color='#b9c8c4',lw=.8)
scopes=[('Membrane orientation','memorient','which way is out?','geometry, not intact-cell exposure'),('Conformational state','state-atlas','active-like or inactive-like?','resemblance, not activity'),('Functional site','site-context · actstate','is the site intact, and competent?','annotation, not observed catalysis'),('Assembly interface','assembly-context','what does it touch, and how much is buried?','one coordinate state, not affinity'),('Comparative analysis','sf-csa','what else looks like this — by fold, and by sequence?','a shared fold is not a function')]
ax.plot([2.4,3.05],[2.215,2.215],color=navy,lw=.9);ax.plot([3.05,3.05],[.855,3.575],color=navy,lw=.9)
rend=fig.canvas.get_renderer()
for i,(title,cli,q,limit) in enumerate(scopes):
 b=3.30-i*.68;c=b+.275
 ax.add_patch(FancyBboxPatch((3.38,b),7.62,.55,boxstyle='round,pad=.04',facecolor='#f7f9f8',edgecolor='#d3d9de'))
 ax.annotate('',xy=(3.36,c),xytext=(3.05,c),arrowprops={'arrowstyle':'->','color':navy,'lw':.9})
 tt=ax.text(3.58,c+.10,title,fontsize=11.5,fontweight='bold',color=navy,va='center')
 ax.text(tt.get_window_extent(renderer=rend).transformed(ax.transData.inverted()).x1+.18,c+.10,cli,fontsize=8.5,color='#5f7570',family='DejaVu Sans Mono',va='center')
 ax.text(3.58,c-.15,q,fontsize=9.5,color=navy,va='center')
 ax.text(10.85,c,limit,fontsize=8.5,style='italic',color=gold,ha='right',va='center')
ax.add_patch(FancyBboxPatch((0,-.52),11,.44,boxstyle='round,pad=.02',facecolor='#edf4f2',edgecolor=teal))
ax.text(5.5,-.30,'Evidence is reported separately — no combined score',fontsize=12,fontweight='bold',color=navy,ha='center',va='center')
ax.text(0,-.78,'The five analyses are independently installable and consume the same checksum-bound manifest.\nMissing references, mappings or runtimes stay visible as missing, and no analysis feeds another a verdict.',fontsize=9,color=gold,va='top')
save(fig,'six-evidence-questions')
record={'schema_version':'1.0','figure_source_sha256':{'historical_summary':hashlib.sha256(summary_path.read_bytes()).hexdigest(),'installed_example':hashlib.sha256((qp/'installed-example-STRUCTURE_EVIDENCE.json').read_bytes()).hexdigest(),'installed_residue_table':hashlib.sha256((qp/'installed-example-RESIDUE_QUALITY.tsv').read_bytes()).hexdigest()},'figures':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(F.iterdir()) if p.suffix in {'.png','.svg'}},'limits':'Worked example synthetic; qualification figure historical collection 2.9; new code not scientifically qualified.'}
(F/'FIGURE_PROVENANCE.json').write_text(json.dumps(record,indent=2)+'\n')
print('Regenerated three figures from versioned evidence.')
