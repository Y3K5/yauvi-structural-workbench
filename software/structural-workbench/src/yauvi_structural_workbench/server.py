"""Standalone loopback browser, using the same store and engines as the CLI."""
from __future__ import annotations
import json
import mimetypes
import secrets
import socket
import threading
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from importlib.resources import files
from urllib.parse import urlsplit, unquote
from yauvi_platform.structural_workbench import (AnalysisError, StructuralAnalysisStore, StructuralSourceStore,
    StructuralSourceError, analysis_definitions, structural_source_descriptors, template_artifact)
from .jobs import JobManager

MAX_JSON=1024*1024
MAX_CHUNK=8*1024*1024
MAX_FILE=250*1024*1024

class WorkbenchServer(ThreadingHTTPServer):
    daemon_threads=True
    def __init__(self,address,workspace,allow_reference_fetch=False):
        if address[0] not in {'127.0.0.1','localhost','::1'}:raise AnalysisError('loopback host required')
        if address[0]=='::1':self.address_family=socket.AF_INET6
        self.store=StructuralAnalysisStore(workspace)
        self.jobs=JobManager(self.store)
        self.allow_reference_fetch=allow_reference_fetch
        self.protein_import_lock=threading.Lock()
        self.token=secrets.token_urlsafe(32)
        try:super().__init__(address,Handler)
        except Exception:self.jobs.close();raise
    def server_close(self):
        super().server_close();self.jobs.close()

class Handler(BaseHTTPRequestHandler):
    def log_message(self, *_):pass
    def _guard(self,mutating=False):
        port=self.server.server_port
        allowed={f'127.0.0.1:{port}',f'localhost:{port}',f'[::1]:{port}'}
        host=self.headers.get('Host','')
        if host not in allowed:raise PermissionError('invalid local host')
        origin=self.headers.get('Origin')
        if origin is not None and origin != 'http://'+host:raise PermissionError('same-origin requests required')
        if self.headers.get('Sec-Fetch-Site')=='cross-site':raise PermissionError('cross-site request refused')
        if mutating and not secrets.compare_digest(self.headers.get('X-Yauvi-Token',''),self.server.token):
            raise PermissionError('valid local session required')
    def _body(self,raw=False):
        if self.headers.get('Transfer-Encoding'):raise ValueError('transfer encoding is unsupported')
        length=int(self.headers.get('Content-Length','0'))
        if length<1 or length>(MAX_CHUNK if raw else MAX_JSON):raise ValueError('request is empty or too large')
        self.connection.settimeout(20)
        body=self.rfile.read(length)
        if len(body)!=length:raise ValueError('incomplete request body')
        if raw:return body
        if self.headers.get_content_type()!='application/json':raise ValueError('application/json required')
        value=json.loads(body)
        if not isinstance(value,dict):raise ValueError('JSON object required')
        return value
    def _send(self,status,body,kind='application/json',*,attachment=None):
        if kind=='application/json' and not isinstance(body, bytes):body=json.dumps(body,sort_keys=True).encode()
        self.send_response(status);self.send_header('Content-Type',kind);self.send_header('Content-Length',str(len(body)))
        self.send_header('Cache-Control','no-store');self.send_header('X-Content-Type-Options','nosniff')
        self.send_header('Referrer-Policy','no-referrer')
        self.send_header('Content-Security-Policy',"default-src 'self'; script-src 'self'; style-src 'self' 'unsafe-inline'; img-src 'self' data:; connect-src 'self'; object-src 'none'; base-uri 'none'; frame-ancestors 'none'")
        if attachment:self.send_header('Content-Disposition',f'attachment; filename="{attachment}"')
        self.end_headers();self.wfile.write(body)
    def do_GET(self):self._dispatch(False)
    def do_POST(self):self._dispatch(True)
    def do_PUT(self):self._dispatch(True)
    def _dispatch(self,mutating):
        try:
            self._guard(mutating)
            parts=[unquote(p) for p in urlsplit(self.path).path.strip('/').split('/') if p]
            if any(p in {'.','..'} or '/' in p or '\\' in p for p in parts):raise ValueError('invalid path')
            if mutating:return self._mutate(parts)
            return self._get(parts)
        except PermissionError as exc:self._send(403,{'error':str(exc)})
        except (AnalysisError,StructuralSourceError,ValueError,KeyError,TypeError) as exc:self._send(400,{'error':str(exc)})
        except (OSError,TimeoutError):self._send(400,{'error':'request or artifact unavailable'})
        except Exception:self._send(500,{'error':'unexpected local request failure; inspect the case and job records before retrying'})
    def _get(self,p):
        store=self.server.store
        if p==['api','session']:return self._send(200,{'token':self.server.token,'reference_fetch_enabled':self.server.allow_reference_fetch})
        if p==['api','definitions']:return self._send(200,analysis_definitions())
        if p==['api','sources']:return self._send(200,structural_source_descriptors())
        if p==['api','analyses']:return self._send(200,store.list())
        if p==['api','jobs']:return self._send(200,self.server.jobs.list())
        if len(p)==3 and p[:2]==['api','templates']:
            name,kind,content=template_artifact(p[2]);return self._send(200,content,kind,attachment=name)
        if len(p)==3 and p[:2]==['api','analyses']:return self._send(200,store.snapshot(p[2]))
        if len(p)==5 and p[:2]==['api','analyses'] and p[3]=='runs':
            return self._send(200,store.snapshot(p[2],run_id=p[4]))
        if len(p)==5 and p[:2]==['api','analyses'] and p[3]=='inputs':
            path=store.input_path(p[2],p[4]);return self._send(200,path.read_bytes(),'text/plain')
        if len(p)==6 and p[:2]==['api','analyses'] and p[3]=='artifacts':
            if p[5] not in {'REPORT_DATA.json','REPORT.html','RAW_EVIDENCE.zip','CHECKSUMS.json','RUN_MANIFEST.json'}:raise ValueError('unknown report artifact')
            path=store.artifact_path(p[2],p[4],p[5]);return self._send(200,path.read_bytes(),'application/octet-stream',attachment=p[5])
        if len(p)==5 and p[:2]==['api','analyses'] and p[3]=='oriented-structure':
            path=store.artifact_path(p[2],p[4],'outputs/memorient/ORIENTED_STRUCTURE.pdb')
            return self._send(200,path.read_bytes(),'chemical/x-pdb',attachment='ORIENTED_STRUCTURE.pdb')
        if len(p)==5 and p[:2]==['api','analyses'] and p[3]=='membrane-layer':
            path=store.artifact_path(p[2],p[4],'outputs/memorient/MEMBRANE_LAYER.json')
            return self._send(200,path.read_bytes(),'application/json',attachment='MEMBRANE_LAYER.json')
        name='/'.join(p) if p else 'index.html'
        if name not in {'index.html','app.js','style.css','vendor/3Dmol-min.js','vendor/membrane-bilayer.js','vendor/3DMOL-LICENSE.txt'}:return self._send(404,{'error':'not found'})
        resource=files('yauvi_structural_workbench').joinpath('ui',*name.split('/'))
        return self._send(200,resource.read_bytes(),mimetypes.guess_type(name)[0] or 'application/octet-stream')
    def _mutate(self,p):
        store=self.server.store
        if len(p)==4 and p[:2]==['api','ingests'] and self.command=='PUT':
            return self._send(200,store.ingest_chunk(p[2],int(p[3]),self._body(raw=True)))
        data=self._body()
        if p==['api','proteins'] or (len(p)==4 and p[:2]==['api','proteins']):
            from yauvi_platform.structural_workbench.protein_import import ProteinImportStore
            if p==['api','proteins'] or p[3]=='prepare':
                if data.get('allow_public_fetch') is not True:
                    raise PermissionError('Choose Retrieve protein files to allow this public-accession download.')
            if not self.server.protein_import_lock.acquire(blocking=False):
                raise AnalysisError('A protein import is already in progress. Please wait for it to finish.')
            try:
                importer=ProteinImportStore(store)
                if p==['api','proteins']:return self._send(201,importer.lookup(data['link']))
                if p[3]=='prepare':return self._send(200,importer.prepare(p[2]))
                if p[3]=='adopt':
                    if any(j['analysis_id']==data['analysis_id'] and j['state'] in {'queued','running'} for j in self.server.jobs.list()):
                        raise AnalysisError('Wait for this analysis to finish before changing its inputs.')
                    return self._send(200,importer.adopt(p[2],data['analysis_id'],data['revision_sha256']))
                raise ValueError('Unknown protein import action.')
            finally:self.server.protein_import_lock.release()
        if p==['api','examples','structure-qc']:
            from .example import create_example
            return self._send(201,create_example(store.workspace,'example-'+secrets.token_hex(6),without_validation=data.get('without_validation') is True))
        if p==['api','analyses']:
            return self._send(201,store.create(data['analysis_id'],analysis_type=data['analysis_type'],question=data['question'],subject_id=data.get('subject_id','')))
        if len(p)==4 and p[:2]==['api','analyses']:
            aid,action=p[2:]
            if action=='parameters':return self._send(200,store.update_parameters(aid,data['parameters']))
            if action=='validate':return self._send(200,store.preflight(aid))
            if action=='run':return self._send(202,self.server.jobs.submit(aid,experimental_acknowledged=data.get('experimental_acknowledged') is True))
            if action=='ingest':
                size=data['size']
                if type(size)!=int or not 0<size<=MAX_FILE:raise ValueError('file limit is 250 MiB')
                return self._send(201,store.begin_ingest(aid,role=data['role'],file_name=data['file_name'],size=size,expected_sha256=data.get('sha256','')))
            if action=='adopt':
                sources=StructuralSourceStore(store.workspace)
                return self._send(200,store.adopt_source_artifact(aid,role=data['role'],artifact_manifest=sources.load(data['acquisition_id'])['artifact'],path=sources.artifact_path(data['acquisition_id'])))
        if len(p)==4 and p[:2]==['api','ingests'] and p[3]=='finalize':return self._send(200,store.finalize_ingest(p[2]))
        if len(p)==4 and p[:2]==['api','jobs'] and p[3]=='cancel':return self._send(200,self.server.jobs.cancel(p[2]))
        if p==['api','sources']:
            if not self.server.allow_reference_fetch:raise PermissionError('reference retrieval disabled; restart with --allow-reference-fetch to opt in')
            return self._send(201,StructuralSourceStore(store.workspace).acquire(data['artifact_type'],data['identifier']))
        return self._send(404,{'error':'unknown action'})

def serve(workspace,host='127.0.0.1',port=8931,allow_reference_fetch=False):
    server=WorkbenchServer((host,port),workspace,allow_reference_fetch)
    hostname='[::1]' if host=='::1' else host
    print(f'YAUVI local workbench: http://{hostname}:{server.server_port}',flush=True)
    try:server.serve_forever()
    except KeyboardInterrupt:pass
    finally:server.server_close()
    return 0
