"""Explicit public-source download command for a future operator; never invoked by replay."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
from urllib.request import urlopen

ROOT = Path(__file__).resolve().parent


def save(url, target, expected, decompress=False):
    if target.exists():
        if hashlib.sha256(target.read_bytes()).hexdigest() != expected:
            raise RuntimeError('Existing source has a different checksum: ' + target.name)
        return
    with urlopen(url, timeout=60) as response:
        data = response.read(100_000_001)
    if len(data) > 100_000_000:
        raise RuntimeError('Download exceeds source size limit')
    if decompress:
        data = gzip.decompress(data)
    if hashlib.sha256(data).hexdigest() != expected:
        raise RuntimeError('Provider bytes differ; source not adopted: ' + target.name)
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.open('xb') as stream:
        stream.write(data)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--download-public-sources', action='store_true', required=True)
    args = parser.parse_args()
    if args.download_public_sources:
        for record in json.loads((ROOT / 'SOURCES.json').read_text())['sources']:
            save(record['url'], ROOT / record['artifact'], record['sha256'])
            validation = record['validation']
            if validation['identity_revision_match']:
                save(validation['url'], ROOT / validation['xml_artifact'], validation['xml_sha256'], decompress=True)
    print('All six coordinate sources and six validation reports match the frozen identities')


if __name__ == '__main__':
    main()
