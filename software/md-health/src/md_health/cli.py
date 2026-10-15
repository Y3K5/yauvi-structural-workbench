import argparse
import json
from pathlib import Path
import sys

from .core import InputError, LIMIT, assess, inspect_bundle, load, prepare, write_report
from .runner import pilot, check_trajectory


def describe():
    return {'module_id':'md_preparation_health','component_id':'md_health','version':'0.1.0','readiness':'prototype',
            'display_name':'GROMACS MD preparation and early health','commands':['describe','validate','prepare','pilot','check-trajectory','assess'],
            'inputs':['GRO coordinates','topology with relative include dependencies','supplied MDP','declared input manifest','optional log/metric XVG and log-bound exit receipt'],
            'outputs':['PREPARATION.json','pilot.mdp','immutable input copies','engine execution receipts','REPORT.json','REPORT.html','CHECKSUMS.json'],
            'limitations':[LIMIT,'No force-field, protonation, solvation, equilibration or automatic repair is performed.',
                           'GROMACS runtime required only for explicit pilot; text/diagnostic checks run offline.',
                           'v1 supports single-process gmx, self-contained topology bundles and GRO; GUI analysis registration is a later integration.']}


def main(argv=None):
    p=argparse.ArgumentParser(prog='md-health',description='Prepare immutable GROMACS pilot inputs and inspect early failure evidence.')
    sub=p.add_subparsers(dest='action',required=True);sub.add_parser('describe')
    for name in ['validate','prepare']:
        s=sub.add_parser(name);s.add_argument('--bundle',type=Path,required=True);s.add_argument('--manifest',type=Path,required=True)
        if name=='prepare':s.add_argument('--out',type=Path,required=True);s.add_argument('--pilot-steps',type=int,required=True)
    s=sub.add_parser('pilot');s.add_argument('--prepared',type=Path,required=True);s.add_argument('--out',type=Path,required=True);s.add_argument('--gmx',default='gmx');s.add_argument('--timeout-seconds',type=int,default=60);s.add_argument('--threads',type=int,default=1)
    s=sub.add_parser('assess');s.add_argument('--log',type=Path,required=True);s.add_argument('--series',type=Path);s.add_argument('--policy',type=Path);s.add_argument('--run-receipt',type=Path);s.add_argument('--out',type=Path,required=True)
    s=sub.add_parser('check-trajectory');s.add_argument('--trajectory',type=Path,required=True);s.add_argument('--tpr',type=Path,required=True);s.add_argument('--out',type=Path,required=True);s.add_argument('--gmx',default='gmx');s.add_argument('--timeout-seconds',type=int,default=60)
    a=p.parse_args(argv)
    try:
        if a.action=='describe':value=describe();code=0
        elif a.action=='validate':value=inspect_bundle(a.bundle,load(a.manifest));code=0
        elif a.action=='prepare':value=prepare(a.bundle,load(a.manifest),a.out,a.pilot_steps);code=0
        elif a.action=='pilot':
            value=pilot(a.prepared,a.out,gmx=a.gmx,timeout_seconds=a.timeout_seconds,threads=a.threads)
            code=0 if value['status']=='engine_exited_without_error'else 2
        elif a.action=='check-trajectory':
            value=check_trajectory(a.trajectory,a.tpr,a.out,gmx=a.gmx,timeout_seconds=a.timeout_seconds)
            code=1 if value['exit_code']==0 else 2
        else:
            value=assess(a.log,a.series,load(a.policy)if a.policy else None,a.run_receipt)
            write_report(value,a.out);code=2 if value['status']=='failure_signal_observed'else 0 if value['status']=='no_failure_signal_observed'else 1
        print(json.dumps(value,indent=2,allow_nan=False));return code
    except (InputError,OSError)as exc:
        print('MD Health blocked: '+str(exc),file=sys.stderr);return 2


if __name__=='__main__':raise SystemExit(main())
