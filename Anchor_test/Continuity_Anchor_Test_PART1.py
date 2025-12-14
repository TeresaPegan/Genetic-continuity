# High Coverage Anchor Continuity Test. Testing for a deviation from expected freq of derived sites in a more recent sample, from those in a more ancient sample, conditioning upon heterozygosity of site in the ancient sample. Idea being that drift alone will not affect the expected frequency of derived sites in this case.
# Chromosomes and given anchor inds are implemented using an external bash script 'ScriptLoop_PART1.sh'. 
# A position is only called heterozygous if it has genotype 0/1 or 1/0, if one of the REF or ALT alleles is present with full support as ancestral allele.
# This script creates an outfile for a given anchor ind listing het positions and ref/alt alleles.
# Note this script can take a snpAD or GATK called vcf.

import sys
import gzip
import get_file_name
import random
from zipfile import ZipFile
import pandas as pd
import os

################################ Functions ###################################

def make_out_str(a_list):
    b_str=''
    for x in a_list:
        b_str+=str(x)+'\t'
    return b_str[:-1]

def check_if_pass_coverage(a_coverage,LOW_COV_THRESH,HIGH_COV_THRESH):
    if (a_coverage>LOW_COV_THRESH and a_coverage<HIGH_COV_THRESH):
        return 1
    else:
        return 0

def get_genotype(a_list):
    b_geno = ''
    coverage = 0
    
    for x in a_list:
        d = x.split(':')

        # genotype missing
        if d[0] == './.':
            return [0, '']

        b_geno = d[0]

        # SNPAD should still have DP, just not at d[1]
        # TRY DP FIRST IF PRESENT
        if len(d) > 2 and d[2].isdigit():
            coverage += int(d[2])
            continue

        # fallback: if DP is in d[1] and is numeric
        if len(d) > 1 and d[1].isdigit():
            coverage += int(d[1])
            continue

        # otherwise: cannot interpret coverage
        return [0, '']

    return [coverage, b_geno]

def get_allele_support(a_list,ref_nt,alt_nt):
    ref_cov=0
    alt_cov=0
    for x in a_list:
        d=x.split(':')
        if len(d)>5: # a snpAD vcf has "GT:DP:A:C:G:T:PP:GQ" in sample field
            nt_dict = {'A':2,'C':3,'G':4,'T':5}
            ref_cov = sum(int(x) for x in d[nt_dict[ref_nt]].split(',') if x.isdigit())
            alt_cov = sum(int(x) for x in d[nt_dict[alt_nt]].split(',') if x.isdigit())
        else: # it is a GATK called vcf with "GT:DP" (homozygote or missing), or "GT:AD:DP:GQ:PL" (heterozygote) in sample field
            ref_cov = int(d[1].split(',')[0])
            alt_cov = int(d[1].split(',')[1])
    return [ref_cov,alt_cov]


def check_if_ok_and_get_var_form(anc_nt,ref_nt,alt_nt):
    set1=set([anc_nt,ref_nt,alt_nt]).difference('.')
    if set1.issubset(nt_set):
        if len(set1)==1:
            return 'OK_NO_VARIATION'
        if len(set1)==2:
            return 'OK_POLY'
    return ''

#############################################################################

#Throw an error unless you have at least two arguments to command line
if len(sys.argv)<6:
    sys.exit('Input Error: The Anchor continuity test PART1 requires a command line input of format: python script.py the_chr anchor_ind ancPath cov_path vcf_path outpath')

the_chr=sys.argv[1]
anchorind=sys.argv[2]
ancPath=sys.argv[3]
cov_path=sys.argv[4]
vcf_path=sys.argv[5]
outpath=sys.argv[6]
#print(the_chr)
#print(anchorind)
#input()

###########################################################################################
##########################
##########################

NUCL=['A','C','G','T']
nt_set=set(NUCL)
ANCESTRAL_FILTER=['A','C','G','T']

############################## USER CONSIDER QUAL THRESHOLDS################################# 
# Get coverage distribution for given individuals and define site coverage thresholds based on that
# Upper DP threshold is 95% of coverage dist
# Lower threshold is higher of 5% of coverage dist, or 8X (to ensure reliable diploid calls)

header_list=['freq','cov']
df=pd.read_csv(cov_path, sep="\t", names=header_list)

# filter out the top and bottom 5% of read lengths (as outliers)
df['cum_percent'] = 100*(df.freq.cumsum() / df.freq.sum())
df_filt = df[(df['cum_percent'] > 5) & ~(df['cum_percent'] > 95)]

# set site coverage thresholds
min_depth = 6.0
LOW_COV_THRESH = max(min_depth, min(df_filt['cov']))
HIGH_COV_THRESH = max(df_filt['cov'])
#print(LOW_COV_THRESH, HIGH_COV_THRESH)
############################################################################################

out_dict={the_chr:{}}

with open(ancPath, 'r') as anc_file:
    with gzip.open(vcf_path, 'rt', encoding='utf-8') as AnchorInd:

        # ---- Parse VCF header; find anchor individual column ----
        line = AnchorInd.readline()
        while line.startswith("##"):
            line = AnchorInd.readline()

        # Now at #CHROM line
        header_fields = line.strip().split('\t')
        sample_names = header_fields[9:]

        if anchorind not in sample_names:
            sys.exit(f"ERROR: Anchor individual '{anchorind}' not found in VCF header.")

        anchor_idx = 9 + sample_names.index(anchorind)

        # Start reading data lines
        l = AnchorInd.readline()
        anc_l = anc_file.readline()

        while l and anc_l:
            vcf_data = l.strip().split()
            anc_d = anc_l.strip().split()

            vcf_pos = vcf_data[1]
            anc_pos = anc_d[1]

            # Sync positions between VCF and ancestral states file
            while vcf_pos != anc_pos:
                if int(vcf_pos) < int(anc_pos):
                    l = AnchorInd.readline()
                else:
                    anc_l = anc_file.readline()

                if not l or not anc_l:
                    break

                vcf_data = l.strip().split()
                anc_d = anc_l.strip().split()
                vcf_pos = vcf_data[1]
                anc_pos = anc_d[1]

            if not (l and anc_l and vcf_pos == anc_pos):
                break

            at_pos = int(anc_pos)
            anc_support = anc_d[3]

            # Only consider ancestral states with full support
            if anc_support == '3':
                qual = vcf_data[5]
                if qual != '.' and int(float(qual)) >= 30:
                    flag = vcf_data[6]
                    if flag not in ['FAIL','FAIL1','FAIL2','FAIL3']:

                        anc_nt = anc_d[2]
                        ref_nt = vcf_data[3]
                        alt_nt = vcf_data[4]

                        # --- ONLY USE THE ANCHOR INDIVIDUAL’S GENOTYPE FIELD ---
                        sample_field = [vcf_data[anchor_idx]]

                        coverage, genotype = get_genotype(sample_field)

                        if check_if_pass_coverage(coverage, LOW_COV_THRESH, HIGH_COV_THRESH):
                            if anc_nt in ANCESTRAL_FILTER:
                                var_form = check_if_ok_and_get_var_form(anc_nt, ref_nt, alt_nt)

                                if var_form == 'OK_POLY':
                                    if genotype in ['0/1', '1/0']:
                                        ref_cov, alt_cov = get_allele_support(sample_field, ref_nt, alt_nt)

                                        anchor_frac_ref = ref_cov / coverage
                                        anchor_frac_alt = alt_cov / coverage

                                        if min(anchor_frac_ref, anchor_frac_alt) >= 0.3333:
                                            if anc_nt in [ref_nt, alt_nt]:
                                                if anc_nt == ref_nt:
                                                    out_dict[the_chr][vcf_pos] = [ref_nt, alt_nt]
                                                else:
                                                    out_dict[the_chr][vcf_pos] = [alt_nt, ref_nt]

            # Read next lines
            l = AnchorInd.readline()
            anc_l = anc_file.readline()
            

#Write chr out_dict to text file
outfile = os.path.join(outpath, f"{the_chr}_{anchorind}_HetPos.txt")
with open(outfile, 'w') as outf:
    for pos, alleles in out_dict[the_chr].items():
        outf.write(f"{pos}\t{alleles[0]}\t{alleles[1]}\n")

print(f"Done. Wrote {len(out_dict[the_chr])} heterozygous anchor positions to:")
print(f"  {outfile}")

##############################################################################################
##############################################################################################
##############################################################################################

