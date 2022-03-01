import re

import click
from gffutils.feature import Feature
import hgvs
from hgvs.assemblymapper import AssemblyMapper
import pyfaidx
from pysam import VariantFile

from primer4.methods import sanger, qpcr, mrna
from primer4.space import pythonic_boundaries
from primer4.utils import (
    convert_chrom,
    log,
    manual_c_to_g,
    sync_tx_with_feature_db
    )


class Variant():
    '''
    hdp = JSONDataProvider(['cdot-0.2.1.refseq.grch37_grch38.json.gz'])
    v = Variant('NM_000546.6:c.215C>G', hdp)

    HGVS syntax examples:

    - NM_000546.6:c.215C>G
    - NM_000546.6:c.215_250del
    - NM_000546.6:c.215_250del

    - g.(?_234567)_(345678_?)del           -- deleted exon is (234567, 345678)
    - c.(4071+1_4072-1)_(5154+1_5155-1)del -- deleted exon is (4072, 5154)
    - (6278_6438+69)_(7310-43_7575)del)

    - NM_015015.3:c.2441+1G>A
    '''
    def __init__(self, code, coords, feature_db):
        self.code = code
        self.data = self.parse_code(code)
        self.g = None
        self.chrom = None
        self.start = self.data.posedit.pos.start.base
        self.start_offset = self.data.posedit.pos.start.offset
        self.end = self.data.posedit.pos.end.base
        self.end_offset = self.data.posedit.pos.end.offset        
        self.tx = self.data.ac

        if self.start_offset:
            self.is_coding = False
        else:
            self.is_coding = True

        self.g = self.map_to_genomic(coords)
        self.chrom = self.g.ac
        self.g_start = self.g.posedit.pos.start.base
        self.g_end = self.g.posedit.pos.end.base

        # Overwrite tx version to sync w/ annotation; v.tx and v.data.ac
        self.tx = sync_tx_with_feature_db(self.tx, feature_db)

    def __repr__(self):
        return self.code

    def parse_code(self, code):
        hp = hgvs.parser.Parser()
        try:
            return hp.parse_hgvs_variant(code)
        except hgvs.exceptions.HGVSParseError:
            if ed := is_exon_deletion(code):
                return ed
            else:
                # click.echo(log('Invalid HGVS syntax, exit.'))
                print('Variant cannot be HGVS-parsed!')
                return None

    def map_to_genomic(self, tx_map, genome_version='GRCh37', method='splign'):
        am = AssemblyMapper(
            tx_map,
            assembly_name=genome_version,
            alt_aln_method=method,
            replace_reference=True)
        return am.c_to_g(self.data)


class ExonDelta():
    def __init__(self, code, feature_db):
        self.code = code
        self.tx, self.is_delta, self.data = self.is_exon_delta(code, feature_db)
        
        if not self.data:
            self.is_unique = None
        else:
            self.is_unique = True if len(self.data) == 1 else False
        
        if self.is_unique:
            x = list(self.data)[0]
            self.chrom = x.chrom
            self.g_start = x.start
            self.g_end = x.end


    def __repr__(self):
        return self.data.__repr__()

    def is_exon_delta(self, code, feature_db):
        '''
        NM_000546.6:g.(?_234567)_(345678_?)del
        NM_000546.6:g.(123_234567)_(345678_?)del
        NM_000546.6:c.(4071+1_4072-1)_(5154+1_5155-1)del
    
        NM_000546.6:c.(?_560-1)_(672+1_?)del
    
        g.(?_234567)_(345678_?)del           -- deleted exon is (234567, 345678)
        c.(4071+1_4072-1)_(5154+1_5155-1)del -- deleted exon is (4072, 5154)
    
        is_exon_deletion(code, db)
        # <Feature exon (NC_000017.10:7578177-7578289[-]) at 0x7f9c32908d60>
        '''
        num = r'[\?]?\d*?[+-]?\d*?'
        cap = r'(\d*)[+-]?\d*?'  # cap .. capture
        pattern = rf'(.*):([gc])\.\({num}_{cap}\)_\({cap}_{num}\)(.*)'
    
        m = re.match(pattern, code)
        
        try:
            # Try to parse code
            tx = m.group(1)
            tx = sync_tx_with_feature_db(tx, feature_db)
        except AttributeError:
            # Not an exon delta;
            # AttributeError: 'NoneType' object has no attribute 'group'
            tx = code.split(':')[0]
            tx = sync_tx_with_feature_db(tx, feature_db)
            return tx, False, set()

        is_coding = True if m.group(2) == 'c' else False
        start     = int(m.group(3))  # coding start coord
        end       = int(m.group(4))  # ... end ...

        if is_coding:
            start = manual_c_to_g(tx, start, feature_db)
            end = manual_c_to_g(tx, end, feature_db)
            if start > end: start, end = end, start
            # Otherwise gffutils does not find a .region()

        A = set(feature_db.children(f'rna-{tx}', featuretype='exon', order_by='start'))
        
        B = set(feature_db.region(region=(list(A)[0].chrom, start, end), featuretype='exon'))
        AB = A.intersection(B)
        
        return tx, True, AB


'''
Template:

- gets tx and retrieves/ masks variants
- gets either a variant or and exon (chrom, start, end)
- applies sanger to var or qpcr/ splice to exon

(decouple input from fn applied, bc/ n:n mapping)

    t = Template(v, db)
    t.relative_pos(t.start)
    # 0

    Needs to interact w/

    - variant
    - manual parse

    Needs to contain:

    - tx coords and relative translation
    - target position in tx
    - padded sequences and other masks (cannot think of any but there might be)
    - SNVs

    Then

    t = Template(...)
    constraints = t.search_for('sanger')
    design(t.masked_sequence, constraints, params['sanger'])

    Template(v, db)
    Template(ed, db)
'''
class Template():
    def __init__(self, mutation, feature_db, featuretype='exon'):
        self.data = mutation
        self.type = type(mutation)  # Template(v, db).type == Variant
        self.feat = feature_db[f'rna-{mutation.tx}']
        self.region = feature_db.region(
            region=(self.feat.chrom, self.feat.start, self.feat.end),
            featuretype=featuretype)
        self.start, self.end = pythonic_boundaries(self.feat)
        self.mask = set()
        self.methods = {
            'sanger': sanger,
            'qpcr':   qpcr,
            'mrna':   mrna,
        }
        self.accepted = {
            'sanger': Variant,
            'qpcr': ExonDelta,
            'mrna': ExonDelta,
        }
        self.variation = {}

    def __repr__(self):
        return self.feat.__repr__()

    def __len__(self):
        return len(self.feat)

    # TODO: apply mask as fn
    def get_sequence(self, genome):
        s = genome[self.feat.chrom][self.start:self.end].__str__()
        # assert len(s) == len(self.feat)
        return s

    def relative_pos(self, n):
        # Primer3 needs positions relative to sequence
        return n - self.start
    
    def apply(self, fn, feature_db, params):
        # Check that we apply the right fn to the right data type
        assert self.type == self.accepted.get(fn)
        # Apply fn
        return self.methods[fn](self, feature_db, params)

    def load_variation_(self, databases):
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

    def mask_sequence(self, genome, mask='N', unmasked=''):
        s = self.get_sequence(genome)
        
        if unmasked:
            masked = ''.join([mask if ix in self.mask else unmasked for ix, i in enumerate(s)])
        else:
            masked = ''.join([mask if ix in self.mask else i for ix, i in enumerate(s)])
        return masked.upper()


class PrimerPair():
    '''
    https://stackoverflow.com/questions/1305532/convert-nested-python-dict-to-object/9413295#9413295
    '''
    def __init__(self, d):
        for a, b in d.items():
            if isinstance(b, (list, tuple)):
               setattr(self, a, [PrimerPair(x) if isinstance(x, dict) else x for x in b])
            else:
               setattr(self, a, PrimerPair(b) if isinstance(b, dict) else b)

    def __repr__(self):
        return f'{self.fwd.start}-{self.fwd.end}:{self.rev.start}-{self.rev.end}, loss: {self.penalty}'




'''
# cd tmp/primer/
# conda activate primer
# py
import json

from cdot.hgvs.dataproviders import JSONDataProvider
import gffutils
from pyfaidx import Fasta

from primer4.models import Variant, Template
from primer4.design import design_primers
# pick_primers


fp_data = '/Users/phi/Dropbox/repos/primer4/data'
fp_config = '/Users/phi/Dropbox/repos/primer4/config.json'

fp_genome = f'{fp_data}/GRCh37_latest_genomic.fna'
fp_coords = f'{fp_data}/cdot-0.2.1.refseq.grch37_grch38.json.gz'
fp_annotation = f'{fp_data}/hg19-p13_annotation.db'

fp_snvs_1 = f'{fp_data}/GRCh37_latest_dbSNP_all.vcf.gz'
fp_snvs_2 = f'{fp_data}/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz'
fp_snvs_3 = f'{fp_data}/ESP6500SI-V2-SSA137.GRCh38-liftover.snps_indels.vcf.gz'


genome = Fasta(fp_genome)
hdp = JSONDataProvider([fp_coords])
db = gffutils.FeatureDB(fp_annotation, keep_order=True)

code = 'NM_000546.6:c.215C>G'
method = 'sanger'

with open(fp_config, 'r') as file:
    params = json.load(file)


v = Variant(code, hdp, db)
tmp = Template(v, db)
constraints = tmp.apply(method, db, params)
# ((0, 7544), (7882, 11188))
# (250, 600)

tmp.load_variation_({
    'dbSNP': fp_snvs_1,
    '1000Genomes': fp_snvs_2,
    'ESP': fp_snvs_3
    })


masked = tmp.mask_sequence(genome)
# tmp.mask_sequence(genome, unmasked='-')[:1000]

primers = [p for p in next(design_primers(masked, constraints, params, []))]
# [7461-7481:8019-8039, loss: 0.1454,
#  7521-7541:8076-8096, loss: 0.2891,
#  7408-7429:7940-7960, loss: 1.7632,
#  7498-7520:7996-8017, loss: 4.1442]



'''







'''
Problems encountered/ solved:

- DSL
- SNVs database preparations
- SNVs different chromosome namings and options fields (AF vs. MAF vs. COMMON)
- recursion necessary for multiple pairs (otherwise get 10000 takes forever but bc/ combinatorics still in the same place)
- automatic transcript conversion (no tests yet)
- test suite
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

