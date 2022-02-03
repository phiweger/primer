import datetime
from difflib import get_close_matches
from math import floor
import pdb
# pdb.set_trace()
import sys

import click
import hgvs.parser
# https://github.com/biocommons/hgvs
import primer3


def infer_coordinates(variant, db):
    '''
    Ah!
    
    https://www.biostars.org/p/104870/
    
    - Exons = gene - introns
    - CDS   = gene - introns - UTRs
    
    and by arithmetic:
    
    - CDS   = Exons - UTRs
    
    So need to look at CDS, not exons.
    
    gffutils assigns IDs to features; for CDS it calls the first CDS in a transcript like "cds-NP_612468.1" (note the prefix) and the next cds-NP_612468.1_1 and so on. Tested for:
    
    - NM_000546.6
    - NM_138459.5
    '''
    hp = hgvs.parser.Parser()
    try:
        v = hp.parse_hgvs_variant(variant)
    except hgvs.exceptions.HGVSParseError:
        click.echo(log('Invalid HGVS syntax, exit.'))
        sys.exit(-1)

    # NM_000350.3:c.4773+3A>G
    # v.posedit.pos.start
    # BaseOffsetPosition(base=4773, offset=3, ...)
    c_pos = v.posedit.pos.start.base
    # coding pos, base (e.g. 4773 from 4773+3)
    offset = v.posedit.pos.start.offset
    if offset:
        click.echo(log('Non-coding variant'))
    name = v.ac

    # Generate a list of valid IDs so we can validate input:
    IDs = set(i.id.replace('rna-', '') for i in db.features_of_type('mRNA'))

    if not name in IDs:
        m = get_close_matches(name, IDs)
        q = name.split('.')[0]
        
        ID_version = {k: v for k, v in [i.split('.') for i in IDs]}
        v = ID_version[q]
        name = f'{q}.{v}'
        click.echo('\n' + log('Warning!\n'))
        click.echo(f'Transcript ID not found, version mismatch? Will use {name} instead.')
        click.echo(f'Not what you want? How about: {m[0]}, {m[1]} or {m[2]}?')

    chromosome = db[f'rna-{name}'].chrom


    coding = []
    for i in db.children(f'rna-{name}', featuretype='CDS', order_by='start'):
        coding.append(i)
    
    strand = db[f'rna-{name}'].strand
    l = 0
    if strand == '-':
        coding = coding[::-1]  # reverse CDS order
        
        for i in coding:
            if l + len(i) >= c_pos:
                break
            else:
                l += len(i)
        # i.end holds genomic coordinate
        g_pos = i.end - (c_pos - l) + 1
    
    else:
        for i in coding:
            if l + len(i) >= c_pos:
                break
            else:
                l += len(i)
        g_pos = i.start + c_pos - 1


    if offset:
        # The offset is in relation to the coding sequence, eg "+3" means
        # three bases downstream of the previous coding sequence.
        if strand == '-':
            g_pos -= offset
        else:
            g_pos += offset

    click.echo(log(f'Variant on chromosome {chromosome}, g. position {g_pos}'))
    return name, chromosome, g_pos


def exon_boundaries(name, chromosome, g_pos, db):
    '''
    !grep 'exon-NM_000546.6-4' GRCh37_latest_genomic.gff
    NC_000017.10    BestRefSeq  exon    7579312 7579590 .   -   .   ID=exon-NM_000546.6-4;Parent=rna-NM_000546.6;Dbxref=GeneID:7157,Genbank:NM_000546.6,HGNC:HGNC:11998,MIM:191170;gbkey=mRNA;gene=TP53;product=tumor protein p53%2C transcript variant 1;tag=RefSeq Select;transcript_id=NM_000546.6

    ex
    <Feature exon (NC_000017.10:7579312-7579590[-]) at 0x7ff39418fb50>
    '''
    g = db.region(region=(chromosome, g_pos, g_pos + 1), featuretype='exon')
    exons = [i for i in g]

    # A variant can be intronic, in which case we don't find any exon
    if not exons:
        return None

    else:
        for ex in exons:
            # An exon can be part of multiple transcripts; choose the exon
            # annotation particular to the transcript of interest. We assume
            # that the exons of any particular transcript do not overlap.
            if name in ex.id:
                break

    # To span an exon, we need Sanger pads and some search space for the 
    # primers on both sides of the exon.
    if ex.strand == '-':
        start = ex.start - 1
        end = ex.end
    else:
        start = ex.start
        end = ex.end + 1

    return ex, start, end


def variant_context(name, genome, chromosome, g_pos, db, params):

    # Primers need 18-30 nt and then we leave another 30 for Sanger "burn-in";
    # so we'll introduce a padding parameter
    # https://eu.idtdna.com/pages/support/faqs/what-is-the-optimal-length-of-a-primer-
    pad = params['burnin_sanger']
    binding = params['binding_site']  # minimum search space for each primer
    # Minimum, maximum amplicon size
    _, mx = params['size_range_PCR']

    try:
        ex, start, end = exon_boundaries(name, chromosome, g_pos, db)
        # Now put an amplicon-sized window across the mutation, such that the
        # mutation has enough space to both sides. Primer3 can later find 
        # primers in this window subject to the (mn, mx) constraints
        required = binding + pad + len(ex) + pad + binding
    
    except TypeError:
        # TypeError: cannot unpack non-iterable NoneType object
        # Don't ask for permission, ask for forgiveness
        required = 1e10
        # Some insanely large value to trigger variant-centered primer design
        # below.

    click.echo(log('Extracting target sequence; depending on chromosome size might take a while'))
    if mx > required:  # as in the sequence
        click.echo(log('Will span exon'))
        # Mask the entire exon and try to find primers spanning the exon
        # Then fill up left and right to maximum amplicon size
        fill = floor((mx - required) / 2)
        
        # See coordinate discussion above
    
        left  = start - pad - binding - fill
        right =   end + pad + binding + fill
    
        s = genome[chromosome][left:right].__str__()
        
        s1 = s[:binding + fill].upper()
        s2 = 'N' * (len(ex) + (2 * pad))
        s3 = s[len(s1 + s2):].upper()
        template = s1 + s2 + s3
        assert len(s) == len(template)

    else:
        click.echo(log('Could not span exon'))
        # We have to do something else
        fill = floor((mx - (2 * pad)) / 2)
    
        left  = g_pos - pad - fill
        right = g_pos + pad + fill
    
        s = genome[chromosome][left:right].__str__()
    
        s1 = s[:fill].upper()
        s2 = 'N' * (2 * pad + 1)  # add + 1 for g_pos
        s3 = s[len(s1 + s2):].upper()
        template = s1 + s2 + s3
        assert len(s) == len(template)

    return template, (left, right)


def common_variants(template, chromosome, boundaries, variants):
    left, right = boundaries
    # Don't put primers in the vicinity of the mutation
    # common, mask = [], []
    mask = set()
    # Don't put them across positions with common variants
    for i in variants.fetch(chromosome, left, right):
         if i.info.get('COMMON'):
            # print(i)
            # common.append(i.ref)
            mask.add(i.pos - left - 1)
            # print(i.ref, s[i.pos - left - 1])
            # print(genome[chromosome][i.pos - 1], i.ref)
            # Not sure why, but we have to subtract 1 to obtain the same letter
    
    masked = ''.join(
        ['N' if ix in mask else i for ix, i in enumerate(template)])
    
    click.echo(log('Template ("x" means no primers here):\n'))
    print(''.join(['x' if i == 'N' else '-' for i in masked]))
    return masked


def design_primers(method, params, masked):
    '''
    # 'SEQUENCE_EXCLUDED_REGION'
    # https://primer3.org/manual.html#SEQUENCE_EXCLUDED_REGION
    
    # 'SEQUENCE_INCLUDED_REGION'
    
    # 'PRIMER_PRODUCT_SIZE_RANGE'
    # https://github.com/libnano/primer3-py/issues/18
    '''
    size_range = params[f'size_range_{method}']
    
    spec =  [
        {
            'SEQUENCE_TEMPLATE': masked,
        },
        {
            'PRIMER_NUM_RETURN': params['n_return'],
            'PRIMER_MIN_SIZE': params['size_min'],
            'PRIMER_OPT_SIZE': params['size_opt'],
            'PRIMER_MAX_SIZE': params['size_max'],
    
            'PRIMER_MIN_TM': params['tm_min'],
            'PRIMER_OPT_TM': params['tm_opt'],
            'PRIMER_MAX_TM': params['tm_max'],
            'PRIMER_MIN_GC': params['GC_min'],
            'PRIMER_MAX_GC': params['GC_max'],
    
            'PRIMER_MAX_POLY_X': params['homopolymer_max_len'],
            'PRIMER_MAX_END_GC': params['3prime_max_GC'],
    
            'PRIMER_MAX_NS_ACCEPTED': params['Ns_max'],
            
            'PRIMER_PRODUCT_SIZE_RANGE': [size_range],
    
            # defaults, here to be explicit
            'PRIMER_SALT_MONOVALENT': params['salt_monovalent'],
            'PRIMER_SALT_DIVALENT': params['salt_divalent'],
            'PRIMER_DNTP_CONC': params['conc_dNTP'],
        }
    ]

    # https://libnano.github.io/primer3-py/quickstart.html#workflow
    design = primer3.bindings.designPrimers(*spec)
    return design


def log(message):
    now = datetime.datetime.now().strftime('%Y-%m-%dT%H:%M:%S')
    return f'[{now}]\t{message}'