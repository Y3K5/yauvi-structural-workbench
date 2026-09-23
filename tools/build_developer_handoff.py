"""Prepare an explicit source-and-preview handoff; never publish."""
import argparse
from pathlib import Path
import hashlib
import json
import re
import zipfile
from build_product_tester import ROOT, PACKAGES, screen, add_denylist_arguments, denylist_from

NAME = 'yauvi-developer-handoff-2026-09-19'
OUT = ROOT/'build/developer-handoff'

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    add_denylist_arguments(parser)
    denylist, denylist_source = denylist_from(parser.parse_args())
    files = {}
    templates = ROOT/'tools/developer-kit'
    for name in ['README.md','dev.py','HANDOFF_VALIDATION.md']:
        files[name] = (templates/name).read_bytes()
    for name in ['UX_UI_ACCEPTANCE.md','ISSUE_TEMPLATE.md']:
        files['docs/'+name] = (templates/name).read_bytes()
    for rel in ['pyproject.toml','yauvi-structural-workbench/README.md','yauvi-structural-workbench/LICENSE']:
        files['source/'+rel] = (ROOT/rel).read_bytes()
    for project, package in PACKAGES.items():
        folder = ROOT/'software'/project/'src'/package
        for path in sorted(folder.rglob('*')):
            if path.is_file() and '__pycache__' not in path.parts and path.suffix != '.pyc' and path.name != '.DS_Store':
                files['source/'+path.relative_to(ROOT).as_posix()] = path.read_bytes()
    for name in ['test_structural_workbench.py','test_workbench_launch.py','test_protein_import.py','test_reference_context.py']:
        files['source/tests/'+name] = (ROOT/'software/platform/tests'/name).read_bytes()
    kit = ROOT/'build/product-tester/yauvi-product-tester-2026-09-19.zip'
    with zipfile.ZipFile(kit) as archive:
        for name in archive.namelist():
            if not name.endswith('/'):
                files['preview/'+name.split('/',1)[1]] = archive.read(name)
    page = files['preview/showcase/index.html'].decode()
    styles = re.findall(r'<style>([\s\S]*?)</style>', page)
    scripts = re.findall(r'<script>([\s\S]*?)</script>', page)
    assert len(styles)==1 and len(scripts)==2
    template = page.replace(styles[0], '__SHOWCASE_STYLE__', 1)
    template = template.replace('<script>'+scripts[0]+'</script>', '__SHOWCASE_VIEWER__', 1)
    template = template.replace('<script>'+scripts[1]+'</script>', '__SHOWCASE_APP__', 1)
    files['source/showcase/template.html'] = template.encode()
    files['source/showcase/style.css'] = styles[0].encode()
    files['source/showcase/app.js'] = scripts[1].removeprefix('\n').encode()
    files['source/showcase/welcome-template.html'] = (ROOT/'tools/tester-kit/START_HERE.html').read_bytes()
    files['source/tools/build_showcase_snapshot.py'] = (templates/'build_showcase_snapshot.py').read_bytes()
    audit = screen(files, denylist, denylist_source)
    # Keep deliberate generic path-prefix assertions in tests; never scrub them.
    reviewed = []
    unresolved = []
    for finding in audit['findings']:
        if finding['file']=='source/tests/test_workbench_launch.py' and finding['category']=='absolute_local_path' and finding['occurrences']==2:
            data=files[finding['file']].decode()
            assert "assert '/Users/' not in json.dumps(original)" in data
            assert "assert '/home/' not in json.dumps(original)" in data
            reviewed.append(dict(finding, disposition='Retained: generic path-prefix assertions, no personal path or account value.'))
        else:
            unresolved.append(finding)
    audit['reviewed_test_literals'] = reviewed
    audit['unresolved_findings'] = unresolved
    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/'CONTENT_SCREEN.json').write_text(json.dumps(audit,indent=2)+'\n')
    if unresolved:
        print(json.dumps(unresolved,indent=2));raise SystemExit('Unreviewed content findings; package held.')
    files['CONTENT_SCREEN.json'] = (json.dumps(audit,indent=2)+'\n').encode()
    manifest={'package':NAME,'files':[{'path':name,'sha256':hashlib.sha256(data).hexdigest(),'bytes':len(data)} for name,data in sorted(files.items())]}
    files['MANIFEST.json']=(json.dumps(manifest,indent=2)+'\n').encode()
    destination=OUT/NAME
    archive_path=OUT/(NAME+'.zip')
    with zipfile.ZipFile(archive_path,'w',compression=zipfile.ZIP_DEFLATED) as z:
        for name,data in sorted(files.items()):
            path=destination/name;path.parent.mkdir(parents=True,exist_ok=True);path.write_bytes(data)
            mode=0o755 if name.endswith(('.sh','.command')) else 0o644
            path.chmod(mode)
            info=zipfile.ZipInfo(NAME+'/'+name,(2026,9,19,0,0,0));info.external_attr=(0o100000|mode)<<16;info.compress_type=zipfile.ZIP_DEFLATED
            z.writestr(info,data)
    with zipfile.ZipFile(archive_path) as z:
        assert z.testzip() is None
        for entry in manifest['files']:
            assert hashlib.sha256(z.read(NAME+'/'+entry['path'])).hexdigest()==entry['sha256']
    digest=hashlib.sha256(archive_path.read_bytes()).hexdigest()
    (OUT/(NAME+'.zip.sha256')).write_text(digest+'  '+archive_path.name+'\n')
    print(json.dumps({'archive':str(archive_path),'bytes':archive_path.stat().st_size,'files':len(files),'sha256':digest,'reviewed_generic_test_literals':len(reviewed),'unresolved_findings':len(unresolved)},indent=2))

if __name__=='__main__':main()
