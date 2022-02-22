from primer4.space import context


def sanger(template, params): 
    '''
    mask_sanger(v, t, params)
    '''
    mn, mx = params['size_range_PCR']
    burnin = params['burnin_sanger']

    pad = burnin * 2 + 100
    # x2 bc/ burnin on both sides of variant, 100 is space to search primers in 
    ex = context(v, db, 'exon')
    
    if ex and (len(ex) + pad < mx):
        # We can span the exon
        # exon boundary + pad, rest primer
        lb = ex.start - burnin  # left boundary
        rb = ex.end + burnin    # right boundary
    else:
        # Cannot span exon
        # center + pad, rest primer
        
        lb = template.data.g_start - burnin
        rb = template.data.g_end + burnin

    # Primer3 takes (start, length) constraints
    rlb = template.relative_pos(lb)
    rrb = template.relative_pos(rb)
    return {
        'only_here': ((0, rlb), (rrb, len(template) - rrb)),
        'size_range': (mn, mx)
    }


def qpcr():
    pass


def mrna():
    pass