# cd tmp/primer/
# conda activate primer
# py
import json

from cdot.hgvs.dataproviders import JSONDataProvider
import gffutils
from pyfaidx import Fasta

from primer4.models import Variant, ExonDelta, SingleExon, ExonSpread, Template
from primer4.design import design_primers, check_for_multiple_amplicons
from primer4.utils import mask_sequence, reconstruct_mrna
# pick_primers


fp_data = '/path/to/primer4/data'
fp_config = '/path/to/primer4/config.json'

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

vardbs = {
    'dbSNP': fp_snvs_1,
    '1000Genomes': fp_snvs_2,
    'ESP': fp_snvs_3
    }


# Sanger
method = 'sanger'
code = 'NM_000546.6:c.215C>G'      # -- strand -, no offset
# code = 'NM_015015.3:c.2441+1G>A'   # -- strand +, offset 1
# code = 'NM_000546.6:c.672+3C>G'    # -- strand -, offset 3

v = Variant(code, hdp, db)
tmp = Template(v, db)
# TODO: Add this to tests
if tmp.feat.strand == '-':
    assert tmp.c_to_g[v.start] - v.start_offset == v.g_start
else:
    assert tmp.c_to_g[v.start] + v.start_offset == v.g_start
# ((0, 7544), (7882, 11188))
# (250, 600)
tmp.load_variation_(vardbs)
# masked = tmp.mask_sequence(genome)

masked = mask_sequence(tmp.get_sequence(genome), tmp.mask)

# tmp.mask_sequence(genome, unmasked='-')[:1000]
constraints = tmp.apply(method, db, params)
primers = [p for p in next(design_primers(masked, constraints, params, []))]
# [7461-7481:8019-8039, loss: 0.1454,
#  7521-7541:8076-8096, loss: 0.2891,
#  7408-7429:7940-7960, loss: 1.7632,
#  7498-7520:7996-8017, loss: 4.1442]



'''
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
'''


# qPCR, starting from explicit exon annotation
method = 'qpcr'
# code = ('NM_000546.6', 5)
code = ('NM_001145408.2', 6)

v = SingleExon(*code)
tmp = Template(v, db)
tmp.load_variation_(vardbs)
# masked = tmp.mask_sequence(genome)
masked = mask_sequence(tmp.get_sequence(genome), tmp.mask)

primers = []
for constraints in tmp.apply('qpcr', db, params):
    print(constraints)
    x = [p for p in next(design_primers(masked, constraints, params, []))]
    primers.extend(x)
# [10527-10548:10612-10632, loss: 1.1868,
#  10568-10592:10667-10686, loss: 5.3035,
#  10860-10880:10990-11010, loss: 0.4314,
#  10822-10842:10944-10964, loss: 2.2463,
#  10783-10803:10896-10918, loss: 2.619]


# mRNA

# reconstruct mRNA as template
# NM_000546.6:c.215C>G
code = ('NM_000546.6', 6, 7)
method = 'mrna'

v = ExonSpread(*code)
tmp = Template(v, db)

tmp.mrna = reconstruct_mrna(tmp.feat, db, genome, vardbs)
constraints = tmp.apply('mrna', db, params)
masked = tmp.mrna[0]
'''
As a sanity check I can the composed sequence in Blastn and it returned:
"Homo sapiens tumor protein p53 (TP53), transcript variant 1, mRNA" with 100%
identity and 100% query cover.
'''
primers = [p for p in next(design_primers(masked, constraints, params, []))]
# [777-797:888-908, loss: 0.5132,
#  706-727:816-836, loss: 1.474,
#  735-755:861-879, loss: 6.427]



# primers.save('foo')







'''
makeblastdb -in GRCh37_latest_genomic.fna -dbtype nucl

blastn -dust no -word_size 7 -evalue 10 -outfmt 6 -query ../test.fna -db GRCh37_latest_genomic.fna -out result


import pandas as pd
names = 'qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore'.split(' ')
df = pd.read_csv('result', sep='\t', names=names)

df[]

> I have logic that does round-robin pairwise comparisons (i.e. compares every hit to every other hit) to determine if each hit pair could result in an amplicon (on same chromosome, opposite-stranded, perfect match at 3'-ends, plus-strand hit is upstream of minus-strand hit, are separated by <= 1000 bp).


https://bioinformatics.stackexchange.com/questions/7262/hits-in-primer-blast-not-found-with-programmatic-blastn-query
https://www.metagenomics.wiki/tools/blast/blastn-output-format-6
https://www.metagenomics.wiki/tools/blast/megablast

blastn -task megablast -dust no -word_size 15 -evalue 2 -outfmt "6 sam" -query ../test.fna -db GRCh37_latest_genomic.fna -out result -perc_identity 0.9 -strand both && cat result



'''




# TODO: now back translate primer coord from blast


results = check_for_multiple_amplicons(primers, fp_genome)

tmp.g_to_c[7578213]







'''
TODO: 

We now could use pybedtools 

https://daler.github.io/pybedtools/index.html

nearby = genes.closest(intergenic_snps, ...)

if primer coord (start, end) in exon, use g_to_c, else use pybedtools get_closest and then minus the corresponding start coord

nearby = genes.closest(intergenic
'''



'''
{'qseqid': '0605c189-3b95-49ac-a4d3-bc52947a8f0c.fwd',
 'sseqid': 'NC_000002.11',
 'pident': 100.0,
 'length': 16,
 'mismatch': 0,
 'gapopen': 0,
 'qstart': 1,
 'qend': 16,
 'sstart': 71186481,
...

'''



'''
    tmp = tempfile.TemporaryDirectory()
    p = tmp.name
    print(f'Aligning {Path(target).name} to {Path(query).name}')

    steps = [
        f'foldseek createdb {target} {p}/targetDB',
        f'foldseek createdb {query} {p}/queryDB',
        f'foldseek search {p}/queryDB {p}/targetDB {p}/aln {p}/tmp -a --cov-mode {mode} --tmscore-threshold {minscore}',
        f'foldseek aln2tmscore {p}/queryDB {p}/targetDB {p}/aln {p}/aln_tmscore',
        f'foldseek createtsv {p}/queryDB {p}/targetDB {p}/aln_tmscore {p}/aln_tmscore.tsv'
    ]

    command = '; '.join(steps)
    log = subprocess.run(command, capture_output=True, shell=True)
'''







