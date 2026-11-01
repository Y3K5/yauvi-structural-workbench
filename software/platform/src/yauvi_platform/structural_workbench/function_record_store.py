"""Narrow native study custody for the public function-record profile.
Uses the existing studies/revisions and content-addressed object layout. It does
not expose general result intake, case linking, model execution or legacy report
plugins. The full Research Home injects its existing StudyStore instead.
"""
from copy import deepcopy
from io import BytesIO
import hashlib,html,json,mimetypes,zipfile
from .store import AnalysisError,StructuralAnalysisStore,SAFE_ID,SHA256,FIXED_ZIP_TIME,_canonical,_sha_file,_write_json,_serialized

def _id(value):
    if not isinstance(value,str) or not SAFE_ID.fullmatch(value):raise AnalysisError('Invalid native study identifier')
    return value

def _text(value,field,limit=4000):
    if not isinstance(value,str) or not value.strip() or len(value)>limit:raise AnalysisError(f'{field} must be a nonempty bounded string')
    return value.strip()

def _digest(value):return hashlib.sha256(_canonical(value)).hexdigest()

class FunctionRecordStore:
    def __init__(self, workspace):
        self.analysis = StructuralAnalysisStore(workspace)
        self.root = self.analysis.root / 'studies'
        self.root.mkdir(exist_ok=True)


    def _transaction(self):
        return self.analysis._transaction()


    def _path(self, study_id):
        target=self.root / _id(study_id)
        if target.is_symlink():raise AnalysisError('Study shortcuts are not admitted')
        return target


    def load(self, study_id):
        try:
            record = json.loads((self._path(study_id)/'STUDY.json').read_text())
            expected = record['revision_sha256']
            clean = {key: value for key, value in record.items() if key != 'revision_sha256'}
            if _digest(clean) != expected:
                raise AnalysisError('Study revision identity is corrupted')
            frozen = self._path(study_id)/'revisions'/f"{record['revision']:08d}.json"
            if frozen.read_bytes() != _canonical(record):
                raise AnalysisError('Study revision archive differs from current record')
            return record
        except (OSError, ValueError, KeyError) as exc:
            raise AnalysisError('Study is unavailable or its record is invalid') from exc


    def _save(self, record):
        target = self._path(record['study_id'])
        revisions = target/'revisions'; revisions.mkdir(parents=True, exist_ok=True)
        numbers = [int(p.stem) for p in revisions.glob('*.json') if p.stem.isdecimal()]
        record = deepcopy(record); record.pop('revision_sha256', None)
        record['revision'] = max(numbers, default=0) + 1
        record['revision_sha256'] = _digest(record)
        with (revisions/f"{record['revision']:08d}.json").open('xb') as handle:
            handle.write(_canonical(record))
        _write_json(target/'STUDY.json', record)
        return record


    @_serialized
    def create(self, study_id, *, title, question, scope):
        target = self._path(study_id)
        if target.exists():
            raise AnalysisError('Study already exists')
        return self._save({'schema_version': '1.0', 'contract_id': 'study_record',
                           'study_id': study_id, 'title': _text(title, 'Title', 160),
                           'question': _text(question, 'Question'), 'scope': _text(scope, 'Scope'),
                           'subjects': [], 'artifacts': [], 'evidence': [], 'case_links': [],
                           'scientifically_approved': False})


    def list(self):
        return [self.load(p.name) for p in sorted(self.root.iterdir())
                if p.is_dir() and SAFE_ID.fullmatch(p.name) and not self._unwritten_directory(p)]


    @staticmethod
    def _unwritten_directory(path):
        """A writer can stop after mkdir but before its first revision is saved.

        Only an empty directory or one empty, ordinary revisions directory is
        scaffolding. Any pointer, archive, shortcut or other retained content
        must still go through load() and its existing integrity refusals.
        """
        pointer = path / 'STUDY.json'
        if path.is_symlink() or pointer.exists() or pointer.is_symlink():
            return False
        entries = list(path.iterdir())
        if not entries:
            return True
        return (len(entries) == 1 and entries[0].name == 'revisions'
                and not entries[0].is_symlink() and entries[0].is_dir()
                and not any(entries[0].iterdir()))


    @staticmethod
    def _subject(record, subject_id):
        _id(subject_id)
        subject = next((s for s in record['subjects'] if s['subject_id'] == subject_id), None)
        if subject is None:
            raise AnalysisError('Subject does not belong to this study')
        return subject


    @_serialized
    def add_subject(self, study_id, *, subject_id, label, kind):
        record = self.load(study_id); _id(subject_id)
        if kind not in KINDS or any(s['subject_id'] == subject_id for s in record['subjects']):
            raise AnalysisError('Invalid subject kind or duplicate subject')
        record['subjects'].append({'subject_id': subject_id, 'label': _text(label, 'Label', 160),
                                   'kind': kind, 'artifact_refs': []})
        return self._save(record)


    def _object(self, digest):
        path = self.analysis.object_path(digest)
        if _sha_file(path) != digest:
            raise AnalysisError('Content object checksum mismatch')
        return path


    def artifact(self, study_id, digest):
        """Return only an attached checksum-verified object, never an arbitrary file."""
        record = self.load(study_id)
        meta = next((a for a in record['artifacts'] if a['sha256'] == digest), None)
        if meta is None:
            raise AnalysisError('Artifact does not belong to this study')
        return meta, self._object(digest).read_bytes()


    @staticmethod
    def _attach(record, subject, *, filename, digest, size, origin):
        artifact_id = 'artifact-' + digest
        if not any(a['artifact_id'] == artifact_id for a in record['artifacts']):
            record['artifacts'].append({'artifact_id': artifact_id, 'sha256': digest,
                                        'file_name': filename, 'size_bytes': size, 'source_kind': origin,
                                        'media_type': mimetypes.guess_type(filename)[0] or 'application/octet-stream'})
        if artifact_id not in subject['artifact_refs']:
            subject['artifact_refs'].append(artifact_id)
        return artifact_id


    def snapshot(self,study_id,*,revision_sha256=None):
        record=self.load(study_id)
        if revision_sha256 is not None:
            if not isinstance(revision_sha256,str) or not SHA256.fullmatch(revision_sha256):raise AnalysisError('Invalid study revision identity')
            if record['revision_sha256']!=revision_sha256:
                found=None
                for path in sorted((self._path(study_id)/'revisions').glob('*.json')):
                    candidate=json.loads(path.read_text())
                    if candidate.get('revision_sha256')==revision_sha256:
                        clean={k:v for k,v in candidate.items() if k!='revision_sha256'}
                        if _digest(clean)!=revision_sha256:raise AnalysisError('Study archive identity differs')
                        found=candidate;break
                if found is None:raise AnalysisError('Requested study revision is unavailable')
                record=found
        for artifact in record['artifacts']:self._object(artifact['sha256'])
        if record['case_links']:raise AnalysisError('Case-linked studies require the general native study workspace')
        return {'study':record,'case_links':[],'adapters':[],'status':'review_required','scientifically_approved':False,
                'limitations':['Narrow function-investigation custody only. No result intake or scientific execution.']}

    @_serialized
    def bundle(self,study_id,*,revision_sha256=None):
        snapshot=self.snapshot(study_id,revision_sha256=revision_sha256);record=snapshot['study']
        f=record.get('protein_function',{})
        if f.get('schema')!='protein-function-investigation-v1':raise AnalysisError('Only function investigations can be exported by this profile')
        from .function_records import report_section
        packet=json.loads(self._object(f['packet_sha256']).read_bytes())
        if packet['schema']!=f['schema'] or packet['function_record']['identity']['key']!=f['key']:raise AnalysisError('Investigation packet identity differs')
        report='<!doctype html><html lang="en"><meta charset="utf-8"><title>Protein function investigation</title><body>'+report_section(record,packet)+'</body></html>'
        entries={'STUDY_RECORD.json':_canonical(record),'STUDY_REPORT.json':_canonical(snapshot),'STUDY_REPORT.html':report.encode()}
        total=sum(len(b) for b in entries.values())
        for item in record['artifacts']:
            path=self._object(item['sha256']);total+=path.stat().st_size
            if total>16*1024*1024:raise AnalysisError('Investigation export exceeds its bounded public profile')
            entries['objects/'+item['sha256']]=path.read_bytes()
        entries['CHECKSUMS.json']=_canonical({'schema_version':'1.0','files':{n:hashlib.sha256(b).hexdigest() for n,b in entries.items()}})
        out=BytesIO()
        with zipfile.ZipFile(out,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=9) as archive:
            for name,data in sorted(entries.items()):
                info=zipfile.ZipInfo(name,FIXED_ZIP_TIME);info.compress_type=zipfile.ZIP_DEFLATED;info.external_attr=0o100644<<16;archive.writestr(info,data)
        return out.getvalue()
