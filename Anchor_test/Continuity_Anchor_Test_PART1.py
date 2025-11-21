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
    b_geno=''
    coverage=0
    for x in a_list:
        d=x.split(':')
        if d[0]=='./.':
            return [0,'']
        b_geno=d[0]
        if len(d)>5: # a snpAD vcf has "GT:DP:A:C:G:T:PP:GQ" in sample field
            coverage+=int(d[1])
        else: # it is a GATK called vcf with "GT:DP" (homozygote or missing), or "GT:AD:DP:GQ:PL" (heterozygote) in sample field
            if len(d)>2:
                if d[2]=='.':
                    return [0,'']
                coverage+=int(d[2])
            else:
                coverage+=int(d[1])
    return [coverage,b_geno]

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
if len(sys.argv)<5:
    sys.exit('Input Error: The Anchor continuity test PART1 requires a command line input of format: python script.py the_chr anchor_ind ancPath vcf_path outpath')

the_chr=sys.argv[1]
anchorind=sys.argv[2]
ancPath=sys.argv[3]
vcf_path=sys.argv[4]
outpath=sys.argv[5]
#print(the_chr)
#print(anchorind)
#input()

###########################################################################################

#Get vcf file name
file_dict = get_file_name.get_name_file_dict()
anchorind_vcf = vcf_path+'/'+file_dict[anchorind]
in_vcf = anchorind_vcf.split('.vc')
anchor_vcf = in_vcf[0] + the_chr + '.vc' + in_vcf[1]

##########################
##########################

NUCL=['A','C','G','T']
nt_set=set(NUCL)
ANCESTRAL_FILTER=['A','C','G','T']

############################## USER CONSIDER QUAL THRESHOLDS################################# 
# Get coverage distribution for given individuals and define site coverage thresholds based on that
# Upper DP threshold is 95% of coverage dist
# Lower threshold is higher of 5% of coverage dist, or 8X (to ensure reliable diploid calls)

header_list=['cov','freq']
df=pd.read_csv('DIR_cov_distr/'+anchorind, sep="\t", names=header_list)

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

#Begin reading in VCF and ANC files, and ensure lines in-sync
with ZipFile(ancPath+'/Ancestral_states.zip','r') as z:
    with z.open('Ancestral_states/chr'+the_chr+'.txt','r') as anc_file:
        with gzip.open(anchor_vcf,'rt',encoding='utf-8') as AnchorInd:
            l='##'
            while l[0]=='#':
                l=AnchorInd.readline()
            anc_l=anc_file.readline().decode('utf-8')
            while l and anc_l:
                vcf_data=l.strip().split()
                anc_d=anc_l.strip().split()
                vcf_pos=vcf_data[1]
                anc_pos=anc_d[0]
        
                while not vcf_pos==anc_pos:
                    if int(vcf_pos) == min(int(vcf_pos), int(anc_pos)):
                        l=AnchorInd.readline()
                    elif int(anc_pos) == min(int(vcf_pos), int(anc_pos)):
                        anc_l = anc_file.readline().decode('utf-8')
                    if l and anc_l:
                        vcf_data = l.strip().split()
                        vcf_pos = vcf_data[1]
                        anc_d = anc_l.strip().split()
                        anc_pos = anc_d[0]
                    else:
                        break
                if not anc_pos == vcf_pos:   # ugly way to catch mismatched positions at end of files and allow to finish without error
                    break

                at_pos=int(float(anc_pos))
                anc_support=anc_d[2]
                if anc_support=='3':    #ANC has full support in outgroup
                    qual=vcf_data[5]
                    if not qual=='.' and int(float(qual)) >= 30:
                        flag=vcf_data[6]
##############################CONSIDER IF YOU WANT OTHER FILTERS#############################
                        if not (flag in ['FAIL','FAIL1','FAIL2','FAIL3']):
                            anc_nt=anc_d[1]
                            ref_nt=vcf_data[3]
                            alt_nt=vcf_data[4]
                            [coverage,genotype]=get_genotype(vcf_data[9:])
                            if check_if_pass_coverage(coverage,LOW_COV_THRESH,HIGH_COV_THRESH):
#############################################################################################
                                if anc_nt in ANCESTRAL_FILTER:
                                    var_form=check_if_ok_and_get_var_form(anc_nt,ref_nt,alt_nt)
                                    if var_form=='OK_POLY':
                                        if (genotype in ['1/0','0/1']): # site is heterozygote
                                            [ref_coverage,alt_coverage]=get_allele_support(vcf_data[9:],ref_nt,alt_nt) # find coverages of alleles
                                            anchor_frac_ref = ref_coverage/coverage
                                            anchor_frac_alt = alt_coverage/coverage 
                                            if not min(anchor_frac_ref,anchor_frac_alt)<0.3333: # ensure proportion of alleles at least a third (to reduce risk of sequencing errors)         
                                                if (anc_nt in [ref_nt,alt_nt]):
                                                    if anc_nt==ref_nt: # if the vcf site is heterozygous with the ref allele matching the ancestral state
                                                        out_dict[the_chr].update({vcf_pos: [ref_nt,alt_nt]})
                                                    else: # heterozygous with the alt allele matching ancestral state
                                                        out_dict[the_chr].update({vcf_pos: [alt_nt,ref_nt]})
                l=AnchorInd.readline()
                anc_l=anc_file.readline().decode('utf-8')



#Write chr out_dict to text file
with open(os.path.join(outpath, f'chr{the_chr}_{anchorind}_HetPos.txt'), 'w') as outf:    #write chr pos to chr file
    for k,v in out_dict[the_chr].items():
        out_str=k+'\t'+'\t'.join(v)
        outf.write(out_str+'\n')


##############################################################################################
##############################################################################################
##############################################################################################

