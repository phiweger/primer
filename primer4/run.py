# cd tmp/primer/
# conda activate primer
# py
import json

from cdot.hgvs.dataproviders import JSONDataProvider
import gffutils
from pyfaidx import Fasta

from primer4.models import Variant, ExonDelta, SingleExon, Template
from primer4.design import design_primers
# pick_primers


fp_data = '/Users/phi/Dropbox/repos/primer4/data'
fp_config = '/Users/phi/Dropbox/repos/primer4/config.json'

fp_genome = f'{fp_data}/GRCh37_latest_genomic.fna'
fp_coords = f'{fp_data}/cdot-0.2.1.refseq.grch37_grch38.json.gz'
fp_annotation = f'{fp_data}/hg19-p13_annotation_bak.db'

fp_snvs_1 = f'{fp_data}/GRCh37_latest_dbSNP_all.vcf.gz'
fp_snvs_2 = f'{fp_data}/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz'
fp_snvs_3 = f'{fp_data}/ESP6500SI-V2-SSA137.GRCh38-liftover.snps_indels.vcf.gz'


genome = Fasta(fp_genome)
hdp = JSONDataProvider([fp_coords])
db = gffutils.FeatureDB(fp_annotation, keep_order=True)

with open(fp_config, 'r') as file:
    params = json.load(file)


# Sanger
code = 'NM_000546.6:c.215C>G'      # -- strand -, no offset
# code = 'NM_015015.3:c.2441+1G>A'   # -- strand +, offset 1
# code = 'NM_000546.6:c.672+3C>G'    # -- strand -, offset 3
method = 'sanger'

v = Variant(code, hdp, db)
tmp = Template(v, db)
# TODO: Add this to tests
if tmp.feat.strand == '-':
    assert tmp.c_to_g[v.start] - v.start_offset == v.g_start
else:
    assert tmp.c_to_g[v.start] + v.start_offset == v.g_start
# ((0, 7544), (7882, 11188))
# (250, 600)
tmp.load_variation_({
    'dbSNP': fp_snvs_1,
    '1000Genomes': fp_snvs_2,
    'ESP': fp_snvs_3
    })
masked = tmp.mask_sequence(genome)
# tmp.mask_sequence(genome, unmasked='-')[:1000]
constraints = tmp.apply(method, db, params)
primers = [p for p in next(design_primers(masked, constraints, params, []))]
# [7461-7481:8019-8039, loss: 0.1454,
#  7521-7541:8076-8096, loss: 0.2891,
#  7408-7429:7940-7960, loss: 1.7632,
#  7498-7520:7996-8017, loss: 4.1442]


# qPCR, starting from HGVS code
method = 'qpcr'
code = 'NM_000546.6:c.(?_560-1)_(672+1_?)del'
v = ExonDelta(code, db)

if len(v.data) == 1:
    exon = next(iter(v.data))
    n_exon = int(exon.id.split('-')[-1])

# For qPCR, the fwd primer can be in the 5' intron or exon, the rev in exon or
# 3' intron, so we have two sets of constraints.
primers = []
for constraints in tmp.apply('qpcr', db, params, n_exon):
    x = [p for p in next(design_primers(masked, constraints, params, []))]
    primers.extend(x)


# qPCR, starting from explicit exon annotation
method = 'qpcr'
# code = ('NM_000546.6', 5)
code = ('NM_001145408.2', 6)
v = SingleExon(*code)

tmp = Template(v, db)
tmp.load_variation_({
    'dbSNP': fp_snvs_1,
    '1000Genomes': fp_snvs_2,
    'ESP': fp_snvs_3
    })
masked = tmp.mask_sequence(genome)

primers = []
for constraints in tmp.apply('qpcr', db, params):
    print(constraints)
    x = [p for p in next(design_primers(masked, constraints, params, []))]
    primers.extend(x)







# mRNA

# reconstruct mRNA as template

def get_mrna(tx, feature_db):
    exons = {}
    for e in feature_db.children(
        f'rna-{tx}', featuretype='exon', order_by='start'):
        exons[int(e.id.split('-')[-1])] = e

    reconstruction = ''
    coords = []
    segmentation = []

    for k in sorted(exons.keys()):
        ex = exons[k]
        seq = ex.sequence(genome).upper()  # accounts for strand
        reconstruction += seq
        # Validated manually using screed that this considers strand, ie for "-"
        # we get the revcomp sequence.
    
        pos = list(range(ex.start, ex.end+1))
        if ex.strand == '-':
            pos = list(reversed(pos))
            # +1 bec intervals INCLUDE the last position, eg 7573008 below, but
            # the reversed fn() excludes it:
            # <Feature exon (NC_000017.10:7571739-7573008[-]) at 0x7fb7fe96eca0>
            # list(reversed(range(ex.start, ex.end)))[0]   is 7573007
            # list(reversed(range(ex.start, ex.end+1)))[0] is 7573008 
        assert len(pos) == len(ex) == len(seq)
    
        segmentation.extend([k] * len(seq))
        coords.extend(pos)
    
    return reconstruction, exons, coords, segmentation


tx = 'NM_000546.6'
a, b, c, d = get_mrna(tx, db)







'''
As a sanity check I can the composed sequence in Blastn and it returned:
"Homo sapiens tumor protein p53 (TP53), transcript variant 1, mRNA" with 100%
identity and 100% query cover.
'''


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

