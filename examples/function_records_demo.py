"""Write an explicitly invented, offline function-record example; no biological data.

Run from a source checkout: python examples/function_records_demo.py OUTPUT_DIR
The output directory must not already exist. Nothing is sent or executed.
"""
import hashlib
import json
from pathlib import Path
import sys


def create(destination):
    root = Path(destination).absolute()
    root.mkdir(parents=True, exist_ok=False)
    data = {'schema': 'functional-protein-catalog-v1', 'dataset_id': 'a' * 64,
            'source_ledger_sha256': 'b' * 64,
            'species': [{'id': 'example', 'label': 'Invented organism'}],
            'groups': [{'id': 'unresolved', 'label': 'Unresolved function'},
                       {'id': 'envelope', 'label': 'Envelope annotation'}],
            'records': [], 'evidence': {}, 'profiles': {},
            'source_relationships': [], 'source_references': []}
    for index, group in enumerate(['unresolved', 'envelope']):
        key = 'example|EXAMPLE' + str(index)
        data['records'].append({'key': key, 'id': 'EXAMPLE' + str(index), 'name': 'Invented protein ' + str(index),
                                'species': 'example', 'strain': 'invented reference', 'primary_group': group,
                                'length': 20, 'sequence_sha256': str(index) * 64, 'holds': ['Invented software example'],
                                'quality_flags': [], 'review_category': 'unevaluated'})
        data['evidence'][key] = {'location': {'status': 'unresolved', 'zone': 'unplaced', 'label': 'Unresolved',
                                             'reason': 'Invented missing location evidence'},
                                'structure': {'status': 'not_joined', 'available': False,
                                              'label': 'No joined structure', 'reason': 'Invented missing form'}}
    raw = json.dumps(data, indent=2).encode()
    (root / 'catalogue.json').write_bytes(raw)
    config = {'schema': 'function-research-config-v1', 'label': 'INVENTED SOFTWARE EXAMPLE',
              'catalogue_path': str(root / 'catalogue.json'), 'catalogue_sha256': hashlib.sha256(raw).hexdigest()}
    (root / 'config.json').write_text(json.dumps(config, indent=2) + '\n')
    print('Invented local example written. Set YAUVI_FUNCTION_RESEARCH_CONFIG to config.json before starting Workbench.')


if __name__ == '__main__':
    if len(sys.argv) != 2:
        raise SystemExit('usage: python examples/function_records_demo.py OUTPUT_DIR')
    create(sys.argv[1])
