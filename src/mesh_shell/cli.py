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
    parser.add_argument('--previous',
                        help='previous mesh JSON (used with --inspection-ledger)')
    parser.add_argument('--inspection-ledger',
                        help='inspection ledger JSON (used with --previous)')
    parser.add_argument('--layer-plan',
                        help='laying work order JSON (z_levels, min_support_ratio)')
    args = parser.parse_args(argv)
    if (args.previous is None) != (args.inspection_ledger is None):
        print('error: --previous and --inspection-ledger must be used together',
              file=sys.stderr)
        return 2
    try:
        data = load(args.source)
        handoff = None
        if args.previous is not None:
            handoff = (load(args.previous), load(args.inspection_ledger))
        plan = load(args.layer_plan) if args.layer_plan is not None else None
        result = process(data, material=args.material, handoff=handoff,
                         layer_plan=plan)
        publish(args.output, result)
    except (OSError, ValueError, UnicodeError, RecursionError) as error:
        print(f'error: {error}', file=sys.stderr)
        return 2
    print(json.dumps(result['summary'], ensure_ascii=False, sort_keys=True))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
