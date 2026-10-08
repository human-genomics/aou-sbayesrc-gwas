"""AoU union phasing and local ancestry: prepare, dry-run, run/resume, status."""
import argparse
from pathlib import Path
import subprocess

from .common import read_json
from .settings import DEFAULT_WORK, PACKAGE


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage', choices=['prepare', 'dry-run', 'run', 'status'])
    parser.add_argument('--run-id', required=True)
    parser.add_argument('--work', type=Path, default=DEFAULT_WORK)
    parser.add_argument('--config', type=Path, default=PACKAGE/'config.json')
    parser.add_argument('--seeds', type=Path, default=PACKAGE/'seeds.json')
    parser.add_argument('--resources', type=Path, default=PACKAGE/'resources.json')
    for name in ('project', 'workspace-bucket', 'pipeline-uri', 'source-uri', 'output-uri',
                 'data-version', 'service-account', 'network', 'subnetwork'):
        parser.add_argument('--'+name)
    parser.add_argument('--region', default='us-central1')
    parser.add_argument('--once', action='store_true', help='Run one coordinator poll; useful for checking a prepared run')
    args = parser.parse_args()
    # Revalidate run IDs even on read/resume paths; never interpret a path as an ID.
    import re
    if not re.fullmatch(r'[a-z][a-z0-9-]{0,47}', args.run_id):
        parser.error('Run ID must start with a letter and use lowercase letters/digits/hyphens (max48)')
    if args.stage == 'prepare':
        from .prepare import prepare
        prepare(args)
    elif args.stage == 'dry-run':
        from .scheduler import command
        manifest = read_json(args.work.resolve()/'runs'/args.run_id/'manifest.json')
        subprocess.run(command(manifest, 'phase', 22)[0]+['--dry-run'], check=True)
    elif args.stage == 'run':
        from .scheduler import loop
        loop(args)
    else:
        from .scheduler import status
        status(args)


if __name__ == '__main__':
    main()
