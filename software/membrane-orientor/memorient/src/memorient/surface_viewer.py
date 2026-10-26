"""Loopback-only surface viewer; fixed assets, no filesystem browsing or writes."""
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
import json
from .pdb_surface import build_scene
from .organism_context import read_marse_context
from .starter_catalogue import catalogue_index, protein_record, organism_overview
from .explicit_scene import load_explicit_scene, model_file
from .md_scene import preparation_status, minimized_scene, minimized_coordinates, minimized_topology, minimized_atom_map

ASSETS = Path(__file__).with_name('viewer_assets')


def get_response(path,manifest=None,workbench=False):
    """Shared read-only routes for standalone and native Workbench viewers."""
    assets={'/':('index.html','text/html; charset=utf-8'),'/index.html':('index.html','text/html; charset=utf-8'),
            '/app.js':('app.js','text/javascript'),'/workbench.js':('workbench.js','text/javascript'),
            '/3Dmol-min.js':('3Dmol-min.js','text/javascript')}
    if path in assets:
        name,mime=assets[path];return 200,(ASSETS/name).read_bytes(),mime
    result=None
    if path=='/connection':result={'workbench':workbench,'has_result':manifest is not None}
    elif path=='/catalogue':result=catalogue_index()
    elif path.startswith('/catalogue/organism/'):result=organism_overview(path.split('/')[-1])
    elif path.startswith('/catalogue/protein/'):result=protein_record(path.split('/')[-1])
    elif path in {'/result','/result/scene'}:
        if manifest is None:result={'stage':'held','detail':'Select a recorded result or start viewer with --manifest'}
        else:
            from .run_bundle import load_manifest,scene
            result=scene(manifest) if path.endswith('/scene') else load_manifest(manifest)
    if result is None:return None
    return 200,json.dumps(result,allow_nan=False).encode(),'application/json'

class ViewerHandler(BaseHTTPRequestHandler):
    def log_message(self, *args):
        pass  # Receptor content and request details never enter a request log.

    def allowed(self):
        host = f'127.0.0.1:{self.server.server_port}'
        return self.headers.get('Host') == host and self.headers.get('Origin', f'http://{host}') == f'http://{host}'

    def respond(self, code, body, mime='application/json'):
        self.send_response(code)
        self.send_header('Content-Type', mime)
        self.send_header('Content-Length', str(len(body)))
        self.send_header('Cache-Control', 'no-store')
        self.send_header('X-Content-Type-Options', 'nosniff')
        self.send_header('Content-Security-Policy', "default-src 'self'; script-src 'self'; style-src 'self' 'unsafe-inline'; connect-src 'self'; img-src 'self' data: blob:; worker-src 'self' blob:; object-src 'none'; frame-ancestors 'none'; base-uri 'none'; form-action 'none'")
        self.end_headers(); self.wfile.write(body)

    def do_GET(self):
        if not self.allowed():
            return self.respond(403, b'{}')
        try:
            response=get_response(self.path,getattr(self.server,'manifest',None))
            if response is not None:return self.respond(*response)
        except KeyError:
            return self.respond(404,b'{"error":"Unknown catalogue record"}')
        except (ValueError,OSError,TypeError):
            return self.respond(503,b'{"error":"Local source or result unavailable or checksum changed"}')
        if self.path in ('/md/ecoli-k12','/md/ecoli-k12/scene','/md/ecoli-k12.gro','/md/ecoli-k12.top','/md/ecoli-k12-map.tsv'):
            try:
                if self.path == '/md/ecoli-k12.gro':
                    return self.respond(200,minimized_coordinates(),'chemical/x-gromos87')
                if self.path == '/md/ecoli-k12.top':
                    return self.respond(200,minimized_topology(),'text/plain; charset=utf-8')
                if self.path == '/md/ecoli-k12-map.tsv':
                    return self.respond(200,minimized_atom_map(),'text/tab-separated-values')
                result = preparation_status() if self.path == '/md/ecoli-k12' else minimized_scene()
                return self.respond(200,json.dumps(result,allow_nan=False).encode())
            except (ValueError,OSError,KeyError):
                return self.respond(503,b'{"error":"Accepted minimized coordinates are not yet available"}')
        if self.path in ('/explicit/ecoli-k12','/explicit/ecoli-k12.cif','/explicit/ecoli-k12-bonds.tsv','/explicit/ecoli-k12-map.tsv'):
            try:
                if self.path == '/explicit/ecoli-k12':
                    return self.respond(200,json.dumps(load_explicit_scene(),separators=(',',':')).encode())
                name={'/explicit/ecoli-k12.cif':'ECOLI_K12_OM.cif','/explicit/ecoli-k12-bonds.tsv':'LIPID_BONDS.tsv','/explicit/ecoli-k12-map.tsv':'ATOM_MAP.tsv'}[self.path]
                return self.respond(200,model_file(name),'chemical/x-mmcif' if name.endswith('.cif') else 'text/tab-separated-values')
            except (ValueError,OSError):
                return self.respond(503,b'{"error":"Explicit model identity check failed or model is absent"}')
        if self.path == '/catalogue' or self.path.startswith('/catalogue/protein/') or self.path.startswith('/catalogue/organism/'):
            try:
                if self.path == '/catalogue': result=catalogue_index()
                elif self.path.startswith('/catalogue/protein/'): result=protein_record(self.path.removeprefix('/catalogue/protein/'))
                else: result=organism_overview(self.path.removeprefix('/catalogue/organism/'))
                return self.respond(200,json.dumps(result,allow_nan=False).encode())
            except KeyError:
                return self.respond(404,b'{"error":"Unknown catalogue record"}')
            except (ValueError,OSError):
                return self.respond(503,b'{"error":"Catalogue source integrity check failed"}')
        assets = {'/':('index.html','text/html; charset=utf-8'), '/app.js':('app.js','text/javascript'),
                  '/3Dmol-min.js':('3Dmol-min.js','text/javascript')}
        if self.path not in assets:
            return self.respond(404, b'{}')
        filename,mime = assets[self.path]
        self.respond(200, (ASSETS/filename).read_bytes(), mime)

    def do_POST(self):
        if not self.allowed():
            return self.respond(403, b'{}')
        if self.path not in ('/place', '/marse-context'):
            return self.respond(404, b'{}')
        try:
            length = int(self.headers.get('Content-Length', '0'))
            if not 0 < length <= 10_000_000:
                raise ValueError('Request limit is 10 MB')
            if self.headers.get('Content-Type') != 'application/json':
                raise ValueError('Expected application/json')
            spec = json.loads(self.rfile.read(length))
            result = read_marse_context(spec['text']) if self.path == '/marse-context' else build_scene(spec)
            payload = json.dumps(result, allow_nan=False).encode()
        except (ValueError, KeyError, TypeError, OverflowError) as exc:
            return self.respond(400, json.dumps({'error':str(exc)}).encode())
        self.respond(200, payload)


def serve(port=8956,manifest=None):
    if manifest is not None:
        from .run_bundle import load_manifest
        load_manifest(manifest)
    server = ThreadingHTTPServer(('127.0.0.1', port), ViewerHandler)
    server.manifest=manifest
    print(f'MemOrient local viewer: http://127.0.0.1:{server.server_port}', flush=True)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        server.server_close()
