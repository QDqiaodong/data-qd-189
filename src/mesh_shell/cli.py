import argparse
import json
import sys

from .io import load, publish
from .service import process


def main(argv=None):
    parser = argparse.ArgumentParser(prog='mesh-shell')
    parser.add_argument('--source', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--material', action='store_true',
                        help='add the material settlement block to the report')
    args = parser.parse_args(argv)
    try:
        result = process(load(args.source), material=args.material)
        publish(args.output, result)
    except (OSError, ValueError, UnicodeError, RecursionError) as error:
        print(f'error: {error}', file=sys.stderr)
        return 2
    print(json.dumps(result['summary'], ensure_ascii=False, sort_keys=True))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
