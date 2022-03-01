from itertools import chain

import primer3

from primer4.models import PrimerPair


def design_primers(masked, constraints, params, previous=[]):
    '''
    # 'SEQUENCE_EXCLUDED_REGION'
    # https://primer3.org/manual.html#SEQUENCE_EXCLUDED_REGION
    
    # 'SEQUENCE_INCLUDED_REGION'
    
    # 'PRIMER_PRODUCT_SIZE_RANGE'
    # https://github.com/libnano/primer3-py/issues/18

    No primers? primer3 excludes Ns by default. Relax this and other conditions:

    - https://primer3.org/manual.html#PRIMER_EXPLAIN_FLAG
    - https://primer3.ut.ee/primer3web_help.htm#PRIMER_MAX_NS_ACCEPTED

    Note that the recursive way we implement this here is greedy, i.e. we take
    the best, then mask, then take the next best, ie no "global optimization".
    '''
    # print('Another round!')
    size_range = constraints['size_range']

    # https://primer3.ut.ee/primer3web_help.htm
    # SEQUENCE_PRIMER_PAIR_OK_REGION_LIST=100,50,300,50 ; 900,60,, ; ,,930,100
    only_here = list(chain(*constraints['only_here']))
    # print(only_here)
    # pdb.set_trace()

    if previous:
        x = previous[-1]
        start, end = x.fwd.start, x.fwd.end
        mask_fwd = masked[:start] + ('N' * (end - start)) + masked[end:]

        start, end = x.rev.start, x.rev.end
        mask_rev = masked[:start] + ('N' * (end - start)) + masked[end:]

        masked = ''.join(
            [i if i==j else 'N' for i, j in zip(mask_fwd, mask_rev)])

    spec =  [
        {
            'SEQUENCE_TEMPLATE': masked,
            'SEQUENCE_PRIMER_PAIR_OK_REGION_LIST': only_here,
        },
        {
            'PRIMER_NUM_RETURN': params['n_return'],
            'PRIMER_MIN_SIZE': params['size_min'],
            'PRIMER_OPT_SIZE': params['size_opt'],
            'PRIMER_MAX_SIZE': params['size_max'],

            'PRIMER_MAX_NS_ACCEPTED': params['Ns_max'],
    
            'PRIMER_MIN_TM': params['tm_min'],
            'PRIMER_OPT_TM': params['tm_opt'],
            'PRIMER_MAX_TM': params['tm_max'],
            'PRIMER_MIN_GC': params['GC_min'],
            'PRIMER_MAX_GC': params['GC_max'],
    
            'PRIMER_MAX_POLY_X': params['homopolymer_max_len'],
            'PRIMER_MAX_END_GC': params['3prime_max_GC'],
            'PRIMER_MAX_END_STABILITY': params['3prime_stability'],
            'PRIMER_MAX_NS_ACCEPTED': params['Ns_max'],
            
            'PRIMER_PRODUCT_SIZE_RANGE': [size_range],
    
            # defaults, here to be explicit
            'PRIMER_SALT_MONOVALENT': params['salt_monovalent'],
            'PRIMER_SALT_DIVALENT': params['salt_divalent'],
            'PRIMER_DNTP_CONC': params['conc_dNTP'],
            'PRIMER_DNA_CONC': params['conc_DNA'],
        }
    ]

    # https://libnano.github.io/primer3-py/quickstart.html#workflow
    design = primer3.bindings.designPrimers(*spec)
    
    try:
        best = PrimerPair(parse_design(design, 1)[0])
        previous.append(best)
        yield from design_primers(masked, constraints, params, previous)

    except KeyError:
        # No primers found
        yield previous


def parse_design(design, n):
    
    primers = {}
    for i in range(n):
        
        fwd_start, fwd_len = design[f'PRIMER_LEFT_{i}']
        rev_start, rev_len = design[f'PRIMER_RIGHT_{i}']
        
        # c .. candidate
        primers[i] = {
            'fwd': {
                'start': fwd_start,
                'end': fwd_start + fwd_len,
                # 'sanity': template[fwd_start:fwd_start + fwd_len],
                'sequence': design[f'PRIMER_LEFT_{i}_SEQUENCE'],
                'Tm': round(design[f'PRIMER_LEFT_{i}_TM'], 2),
            },
            'rev': {
                # That Python 0-based, end-exclusive indexing thing ...
                'start': rev_start - rev_len + 1,
                'end': rev_start + 1,
                # 'sanity': rc(template[rev_start - rev_len + 1:rev_start + 1]),
                'sequence': design[f'PRIMER_RIGHT_{i}_SEQUENCE'],
                'Tm': round(design[f'PRIMER_RIGHT_{i}_TM'], 2),
            },
            'insert': design[f'PRIMER_PAIR_{i}_PRODUCT_SIZE'],
            'penalty': round(design[f'PRIMER_PAIR_{i}_PENALTY'], 4),
        }

    return primers


def sort_penalty(primers):
    '''
    It seems primer3 already sorts primers from best to worst, but just to
    make sure.
    '''
    loss = {}
    for k, v in primers.items():
        loss[k] = v['penalty']
    return [k for k, v in sorted(loss.items(), key=lambda x: x[1])]


# def remove_overlapping(primers, length):
#     occupied = [0] * length
#     for k, v in primers.items():
#         for i in ['fwd', 'rev']:
#             start, end = v[i]['start'], v[i]['end']
#             overlap = sum(occupied[start:end])
#             if overlap < 10:
#                 print(k, i, overlap, v['penalty'])
#             # Not +1, primer3 py package has pythonic coords
#             for j in range(start, end + 1):
#                 occupied[j] = 1



