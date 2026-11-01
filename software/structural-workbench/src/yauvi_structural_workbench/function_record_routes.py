"""Public-profile routes over shared record/investigation cores.

Called only after the native Host/origin/session guards. No provider, model,
analysis or experiment is dispatched by these routes.
"""
from urllib.parse import parse_qs, urlsplit
from yauvi_platform.structural_workbench.function_records import browse, record
from yauvi_platform.structural_workbench.function_investigations import FunctionInvestigations


def get(handler, parts):
    relevant = parts[:2] == ['api', 'protein-functions']
    bundle = len(parts) == 4 and parts[:2] == ['api', 'studies'] and parts[3] == 'bundle'
    if not relevant and not bundle:
        return False
    values = parse_qs(urlsplit(handler.path).query, keep_blank_values=True)
    if any(len(v) != 1 for v in values.values()):
        raise ValueError('Choose one filter or exact revision')
    if parts == ['api', 'protein-functions', 'catalogue'] and not set(values) - {'species', 'group', 'q', 'gap', 'page'}:
        handler._send(200, browse(**{k: v[0] for k, v in values.items()}))
        return True
    if parts in (['api', 'protein-functions', 'record'], ['api', 'protein-functions', 'export']) and set(values) == {'key'}:
        handler._send(200, record(values['key'][0]), attachment='protein-function-record.json' if parts[-1] == 'export' else None)
        return True
    service = FunctionInvestigations(handler.server.store.workspace)
    if parts == ['api', 'protein-functions', 'investigations'] and set(values) == {'key'}:
        handler._send(200, service.list(values['key'][0]))
        return True
    if len(parts) == 4 and parts[2] == 'investigations' and not set(values) - {'revision'}:
        handler._send(200, service.load(parts[3], values.get('revision', [None])[0]))
        return True
    if bundle and not set(values) - {'revision'}:
        handler._send(200, service.studies.bundle(parts[2], revision_sha256=values.get('revision', [None])[0]),
                      'application/zip', attachment='protein-function-investigation.zip')
        return True
    raise ValueError('Unknown function-record request')


def mutate(handler, parts, data):
    if parts[:2] != ['api', 'protein-functions']:
        return False
    if parts != ['api', 'protein-functions', 'investigations'] or handler.command != 'POST':
        raise ValueError('Function investigations require the declared POST route')
    handler._send(201, FunctionInvestigations(handler.server.store.workspace).save(data))
    return True
