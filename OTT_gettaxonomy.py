#!/usr/bin/env python3
import sys
import copy
from opentree import OTCommandLineTool, get_suppressed_taxon_flag_expl_url

cli = OTCommandLineTool(usage='Display taxonomic info (based on the Open Tree '
                              'Taxonomy) for a taxon ID',
                        common_args=("ott-id", ))
cli.parser.add_argument('--input', type=str, required=True,
                        help='path to a 2 coulmn tab-delimited file, the '
                             'first column identifiers and the second '
                             'semicolon-separated taxonomy strings')
cli.parser.add_argument("--context-name", default=None,
                        help='If you know the named Open Tree name searching '
                             'index, you can supply it here to limit your '
                             'search to only a subset of the taxonomny.')
cli.parser.add_argument("--do-approximate-matching", action="store_true",
                        help='Enables fuzzy matching')
cli.parser.add_argument("--include-suppressed", action="store_true",
                        help='Return taxa that are normally suppressed from '
                             'TNRS results. See '
                             f'{get_suppressed_taxon_flag_expl_url()}')
OT, args = cli.parse_cli()

kwargs = {'include_children': False,
          'include_lineage': False,
          'include_terminal_descendants': False,
          }

# Parse the taxstrings

uniqtips = set()
outdata = []
with open(args.input, 'r') as ih:
    for line in ih.readlines():
        idx, taxstring = line.strip().split('\t')
        taxa = taxstring.split(';')
        tip = [t for t in taxa if t][-1]
        outdata.append([idx, taxstring, tip])
        uniqtips.add(tip)

# Search the uniq tax
matchout = OT.tnrs_match(uniqtips,
                         context_name=args.context_name,
                         do_approximate_matching=args.do_approximate_matching,
                         include_suppressed=args.include_suppressed
                         ).response_dict

# Parse the results
defaultout = [0, None, None, None, None, None]

results = dict()
for match in matchout['results']:
    out = copy.deepcopy(defaultout)
    out[0] = len(match['matches'])
    if out[0] == 1:
        dat = match['matches'][0]
        out[1] = dat['unique_name'] if 'unique_name' in dat else match['name']
        out[2] = dat['taxon']['rank']
        out[3] = dat['taxon']['ott_id']
        out[4] = f"ott{out[3]}"
        taxsources = [d.split(':') for d in dat['taxon']['tax_sources']]
        taxsources = {ts[0]: ts[1] for ts in taxsources}
        out[5] = taxsources['ncbi'] if 'ncbi' in taxsources else None
    results[match['name']] = out

# Add in missing results
for missing in [u for u in list(uniqtips) if u not in results]:
    results[missing] = defaultout

# Write out
for line in outdata:
    sys.stdout.write('\t'.join(str(i) for i in line + results[line[2]]) + '\n')