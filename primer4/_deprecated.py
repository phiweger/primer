    def map_to_genomic_manually(self, db):
        # TODO: Not congruent w/ automatic for NM_015015.2:c.2441+1G>A
        name = self.tx
        c_pos = self.base
        offset = self.offset

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

        # click.echo(log(f'Variant on chromosome {chromosome}, g. position {g_pos}'))
        return g_pos


    def load_variation2_(self, databases):
        for name, db in databases.items():
            variants = VariantFile(db)

            # dbSNP names chromosomes like "NC_000007.13", others like "7"
            if name != 'dbSNP':
                chrom = convert_chrom(self.feat.chrom)
            else:
                chrom = self.feat.chrom
            vv = variants.fetch(chrom, self.feat.start, self.feat.end)
            self.variation[name] = vv

            for i in vv:
                # .info.get(...) raises ValueError: Invalid header if not there
                info = dict(i.info)
                
                if name == 'dbSNP':
                    if info.get('COMMON'):
                        self.mask.add(i.pos - self.feat.start - 1)
                
                elif name == '1000Genomes':
                    if info['AF'][0] >= 0.01:
                        self.mask.add(i.pos - self.feat.start - 1)
                
                elif name == 'ESP':
                    if float(info['MAF'][0]) >= 1:
                        self.mask.add(i.pos - self.feat.start - 1)

                else:
                    print(f'"{name}" is not a valid variant database')
        return None










def get_mrna(transcript, feature_db):







# TODO: Have a universal translate fn that translates coords btw/ genomic, transcript, coding and relative -- these need to be object methods of the Template obj.


for i in list(db.children(f'rna-{v.tx}', featuretype='CDS', order_by='start')):
    print(i.frame)




def qpcr(template, feature_db, params):
    '''
    One primer inside the exon, one outside
    '''
    pass


'''
1 inside, 1 outside

should not be too difficult
'''


s = mrna(tmp, db)



mrna = db[f'rna-{v.tx}']

def get_exon(db, name, exon):
    try:
        return db[f'exon-{name}-{exon}']
    except gffutils.exceptions.FeatureNotFoundError:
        return None


# https://pythonhosted.org/gffutils/autodocs/gffutils.FeatureDB.create_introns.html
def get_intron(db, name, exon, relative=-1):



qry = 3
tx = v.tx

ex = retrieve_exon(db, tx, qry)



# https://github.com/seandavi/GFFutils#imputing-introns
exons = {}

for e in db.children('rna-NM_000546.6', featuretype='exon', order_by='start'):
    exons[int(e.id.split('-')[-1])] = e

mrna = ''
for k in sorted(exons.keys()):
    mrna += exons[k].sequence(genome)
    # Validated manually using screed that this considers strand, ie for "-"
    # we get the revcomp sequence.
'''
As a sanity check I can the composed sequence in Blastn and it returned:
"Homo sapiens tumor protein p53 (TP53), transcript variant 1, mRNA" with 100%
identity and 100% query cover.
'''

# qrna
exons = list(db.children('rna-NM_000546.6', featuretype='exon', order_by='start'))
introns = list(db.interfeatures(exons))
# Naturally, we have more exons than introns
assert len(exons) == len(introns) + 1




l = twolists(exons, introns)
ix = [int(i.id.split('-')[-1]) if i.id else None for i in l]

ix.index(1)
# and then l[20+1] and l[20-1]









def qpcr(template, feature_db, params, n_exon):
    '''
    qpcr(tmp, db, params, 5)
    '''
    mn, mx = params['size_range_qPCR']

    exons = list(feature_db.children(
        template.feat.id, featuretype='exon', order_by='start'))
    introns = list(feature_db.interfeatures(exons))
    l = twolists(exons, introns)
    ix = [int(i.id.split('-')[-1]) if i.id else None for i in l]
    mid = ix.index(n_exon)
    left = mid - 1 
    right = mid + 1

    l[left].start
    rlb = template.relative_pos(l[left].start)  # rlb .. relative left bound
    rmb = template.relative_pos(l[mid].start)  # rmb .. mid

    return {
        'only_here': (
            (rlb, len(l[left])),
            (rmb, len(l[mid]))
            ),
        'size_range': (mn, mx)
        }


'''
db = gffutils.create_db(
    fp, 
    dbfn='hg19-p13_annotation.db',
    force=True,
    keep_order=True,
    merge_strategy='merge',
    sort_attribute_values=True)
'''
# https://github.com/daler/gffutils/issues/111
introns = list(db.create_introns())
# 819,463
db.update(
    introns,
    disable_infer_transcripts=True,
    disable_infer_genes=True,
    verbose=True,
    merge_strategy='merge')
# backup=True


# f'exon-{name}-{exon}
# intron = db['exon-NR_026818.1-2', 'exon-NR_026818.1-3']


def mrna(template, feature_db, params):
    '''
    One primer in exon 1, the other in exon 2
    '''
    pass


def get_contraints(x):
    return x.start, x.end - x.start



# TODO: closest features
# https://github.com/seandavi/GFFutils#closest-features


'''
input: just name exon for now, later bold on the syntax to parse automatically

NM_000546.6, exon 6

class mRNA()

has a dict {1: (start, end)}

then select exon + and - (exon +- 1)

place primer there (constrains, see syntax of the primer3 param)

possibility to select only one splice site (eg +1 or -1 if no arg take both)
'''







# -----------------------------------------------------------------------------


'''
Problems encountered/ solved:

- DSL
- SNVs database preparations
- SNVs different chromosome namings and options fields (AF vs. MAF vs. COMMON)
- recursion necessary for multiple pairs (otherwise get 10000 takes forever but bc/ combinatorics still in the same place)
- automatic transcript conversion (no tests yet)
- test suite


also:

visit mibi



'''




# -----------------------------------------------------------------------------

'''
hdp = JSONDataProvider(['cdot-0.2.1.refseq.grch37_grch38.json.gz'])
db = gffutils.FeatureDB('/Users/phi/Dropbox/repos/primer/data/hg19-p13_annotation.db', keep_order=True)
genome = Fasta('/Users/phi/Dropbox/repos/primer/data/GRCh37_latest_genomic.fna')


import json
with open('/Users/phi/Dropbox/repos/primer/config.json', 'r') as file:
    params = json.load(file)


# v = Variant('NM_015015.2:c.2441+1G>A', hdp, db)
v = Variant('NM_000546.6:c.215C>G', hdp, db)
# ex = context(v, db, 'exon')
t = Template(v, db)

# TODO: Template() needs to work w/ ('NM_000546.6', '4') exon coords, too
# We have genomic coords in this case, so good.
# -- Create another class Exon() and do have the same interface, but do the
# manual coord conversion in there.

# Actually, we can parse this rather easily:
# 
# g.(?_234567)_(345678_?)del           -- deleted exon is (234567, 345678)
# c.(4071+1_4072-1)_(5154+1_5155-1)del -- deleted exon is (4072, 5154)
# 
# Genomic we can look up, coding we'd have to translate using existing code.



constraints = mask_sanger(v, t, params)
# TODO: the mask_x fn should not need the variant
# ((0, 7544), (7882, 19070))
design_primers('PCR', params, t.get_sequence(genome), constraints)


# TODO design(search_space, template, params)



method = 'mRNA'

before = neighbor(ex, db, -1)
after = neighbor(ex, db, 1)


design(template, placements, mask, params)


placements = set(
    {'left': (10, 15), 'right': (45, 89)},
    {'left': (45, 89), 'right': (91, 98)},
    )

# then pass each to primer3



Target(v, ex, params, 'mRNA')

'''








##fileformat=VCFv4.2
##fileDate=20210513
##source=dbSNP
##dbSNP_BUILD_ID=155
##reference=GRCh37.p13
##phasing=partial

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