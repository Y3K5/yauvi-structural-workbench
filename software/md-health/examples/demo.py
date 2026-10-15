"""Invented diagnostics only: no GROMACS simulation or biological input."""
import argparse
import hashlib
import json
from pathlib import Path
import sys

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from md_health.core import assess, sha, write_report


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    if a.out.exists():raise SystemExit('Preserve existing demo evidence; choose a new output directory.')
    a.out.mkdir(parents=True)
    policy={'expected_end_ps':1,'required_metrics':['Temperature'],'bands':{'Temperature':[290,310]},
            'basis':'Invented unit-test bands, not a scientific recommendation'}
    cases=[('bounded-clean-fixture','GROMACS synthetic fixture\nLINCS: order 4\nFinished mdrun on fixture\n','0 300\n1 300\n',True),
           ('early-instability-fixture','GROMACS synthetic fixture\nStep 10, LINCS WARNING\nFatal error: constraint failure\n','0 300\n1 nan\n',False),
           ('incomplete-fixture','GROMACS synthetic fixture\nStep 0\n','0 300\n',False)]
    results={}
    for name,log,values,receipt in cases:
        source=a.out/name;source.mkdir();lp=source/'run.log';lp.write_text(log)
        x=source/'metrics.xvg';x.write_text('@ xaxis label "Time (ps)"\n@ s0 legend "Temperature"\n'+values)
        rp=None
        if receipt:
            rp=source/'MOCK_RECEIPT.json';rp.write_text(json.dumps({'kind':'gromacs_mdrun_receipt_v1','exit_code':0,'log_sha256':sha(lp),'synthetic_only':True}))
        result=assess(lp,x,policy,rp);result['synthetic_only']=True;result['real_engine_executed']=False
        write_report(result,source/'report');results[name]=result['status']
    links=''.join(f'<li><a href="{name}/report/REPORT.html">{name.replace("-"," ")}</a> → {status.replace("_"," ")}</li>'for name,status in results.items())
    (a.out/'index.html').write_text('''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>MD Health · software demonstration</title><style>body{max-width:860px;margin:60px auto;padding:0 24px;background:#f5f7f5;color:#20322e;font:18px/1.7 system-ui}h1{line-height:1.2}li{margin:20px 0}a{color:#226c68}.notice{padding:20px;background:white;border-left:4px solid #226c68}</style><p>YAUVI WORKBENCH · MD HEALTH PROTOTYPE</p><h1>Catch early problems.<br>Keep the evidence.</h1><p>Existing GROMACS files → a recorded test copy → an explicit short pilot → diagnostics from logs and named metric series.</p><p class="notice">These three reports use invented test data. No GROMACS simulation ran. A clean short-run report never establishes equilibrium, convergence or biological validity.</p><ul>'''+links+'''</ul><p>The preparation module preserves original coordinates, topology dependencies and run settings. A changed test length is recorded separately. GROMACS performs the authoritative topology/coordinate preprocessing check. Missing inputs, warnings and interrupted runs stay visible.</p><p>Next real-engine check: supply the existing prepared input bundle and an available GROMACS runtime. File conversion, automatic physical-parameter repair and a browser upload/run form are outside this first prototype.</p></html>''')
    record={'synthetic_only':True,'real_engine_executed':False,'results':results,
            'files':{str(f.relative_to(a.out)):sha(f)for f in a.out.rglob('*')if f.is_file()}}
    (a.out/'DEMO_RECEIPT.json').write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps(results))


if __name__=='__main__':main()
