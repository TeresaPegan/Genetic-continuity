# High Coverage Anchor Continuity Test PART2. Testing for a deviation from expected freq of derived sites in a more recent sample, from those in a more ancient sample, conditioning upon heterozygosity of site in the ancient sample. Idea being that drift alone will not affect the expected frequency of derived sites in this case.
# This version takes an anchor individual and a single vcf (more recent individual) as argument to command line.
# The iteration across chromosomes is implemented within an external bash script. 
# A position is only called heterozygous if it has genotype 0/1 or 1/0, if one of the REF or ALT alleles is present with full support as ancestral allele in outgroups.
# This version is the same as 'Continuity_Anchor_Test_PART2.py' except it differentiates between how the transition and transversion het sites in the anchor are found as ancestral, derived and other. 
# Arguments to this python script are supplied from arrays in 'ScriptLoop_PART2.sh'

import sys
import gzip
import get_file_name
import random

################################ Functions ###################################
def check_if_missingness(format_str, sample_str):
    d = sample_str.split(':')
    gt = d[0]
    if gt in ['./.', '.|.']:
        return 1
    return 0

def get_haploid_genotype(ref_nt, alt_nt, format_str, sample_str):
    """
    Extract haploid genotype call and coverage for a single sample in a multi-individual VCF.

    - ref_nt, alt_nt: reference and (single) alternate allele from VCF
    - format_str: FORMAT column (e.g., "GT:AD:DP:GQ:PL" or "GT:DP:A:C:G:T:PP:GQ")
    - sample_str: sample's genotype field corresponding to format_str

    Returns [coverage, haploid_nt] where haploid_nt is one of ref_nt / alt_nt / ''.
    """

    fmt_fields = format_str.split(':')
    d = sample_str.split(':')

    # Missing genotype
    gt = d[0]
    if gt in ['./.', '.|.']:
        return [0, '']

    # Build map FORMAT field -> index
    fmt_idx = {name: i for i, name in enumerate(fmt_fields)}

    # Get coverage (DP if present, else from AD, else 0)
    coverage = 0
    if 'DP' in fmt_idx and fmt_idx['DP'] < len(d) and d[fmt_idx['DP']] not in ['.', '']:
        coverage = int(d[fmt_idx['DP']])
    elif 'AD' in fmt_idx and fmt_idx['AD'] < len(d) and d[fmt_idx['AD']] not in ['.', '']:
        # Sum AD if DP missing
        ad_vals = d[fmt_idx['AD']].split(',')
        coverage = sum(int(x) for x in ad_vals if x not in ['.', ''])
    else:
        coverage = 0

    # Case 1: SNPAD-style VCF with per-nucleotide fields "A","C","G","T"
    if all(x in fmt_idx for x in ['A', 'C', 'G', 'T']):
        # Extract per-nt coverages by letter
        nt_pos = {nt: fmt_idx[nt] for nt in ['A', 'C', 'G', 'T']}
        ref_cov = int(d[nt_pos[ref_nt]]) if ref_nt in nt_pos and d[nt_pos[ref_nt]].isdigit() else 0

        if len(alt_nt) == 1 and alt_nt in nt_pos and d[nt_pos[alt_nt]].isdigit():
            alt_cov = int(d[nt_pos[alt_nt]])
        else:
            # multiallelic or unexpected; skip
            return [0, '']

        if ref_cov == 0 and alt_cov == 0:
            return [coverage, '']
        if ref_cov == 0:
            return [coverage, alt_nt]
        if alt_cov == 0:
            return [coverage, ref_nt]
        # Both present: randomly choose one
        return [coverage, random.choice([ref_nt, alt_nt])]

    # Case 2: GATK-like with AD field ("GT:AD:DP:GQ:PL" or similar)
    if 'AD' in fmt_idx and fmt_idx['AD'] < len(d):
        temp = d[fmt_idx['AD']].split(',')
        if len(temp) >= 2:
            ref_cov = int(temp[0]) if temp[0] not in ['.', ''] else 0
            alt_cov = int(temp[1]) if temp[1] not in ['.', ''] else 0

            if coverage == 0 and (ref_cov + alt_cov) > 0:
                coverage = ref_cov + alt_cov

            if ref_cov == 0 and alt_cov == 0:
                return [coverage, '']
            if ref_cov == 0:
                return [coverage, alt_nt]
            if alt_cov == 0:
                return [coverage, ref_nt]
            return [coverage, random.choice([ref_nt, alt_nt])]

    # Fallback: interpret genotype directly when AD info is not usable
    if gt in ['0/0', '0|0']:
        return [coverage, ref_nt]
    elif gt in ['1/1', '1|1']:
        return [coverage, alt_nt]
    elif gt in ['0/1', '1/0', '0|1', '1|0']:
        return [coverage, random.choice([ref_nt, alt_nt])]
    else:
        return [0, '']



def orient_and_get_count_haploid(haploid_gt,anc_nt,der_nt):
    #print('genotypes:',haploid_genotype,anc_nt,der_nt)
    if haploid_gt==anc_nt:
        if anc_nt+der_nt==('AG' or 'GA' or 'CT' or 'TC'):
            return [1,0,0,0,0,0] # Transition
        else:
            return [0,1,0,0,0,0] # Transversion
    elif haploid_gt==der_nt:
        if anc_nt+der_nt==('AG' or 'GA' or 'CT' or 'TC'):
            return [0,0,1,0,0,0] # Transition
        else:
            return [0,0,0,1,0,0] # Transversion
    elif ((haploid_gt!=anc_nt) and (haploid_gt!=der_nt) and (haploid_gt in nt_set)):
        if anc_nt+der_nt==('AG' or 'GA' or 'CT' or 'TC'):
            return [0,0,0,0,1,0] # Transition
        else:
            return [0,0,0,0,0,1] # Transversion
    else:
        return [0,0,0,0]


def check_if_biallelic(anc_nt,der_nt,ref_nt,alt_nt):
    """
    This function can be used to ensure a site is biallelic with either alt or ref matching the ancestral and/or derived allele
    """
    set1=set([anc_nt,der_nt,ref_nt,alt_nt]).difference('.')
    if set1.issubset(nt_set):
        if len(set1)==2: # this condition means we only consider biallelic sites, where ref or alt allele is same as ancestral or derived
            return 1
        else:
            return 0


def make_out_str(a_list):
    b_str=''
    for x in a_list:
        b_str+=str(x)+'\t'
    return b_str[:-1]

# def check_der(haploid_gt, anc_nt, der_nt):
#     if haploid_gt == anc_nt:
#         pass  # nothing
#     elif haploid_gt == der_nt:
#         print(f"derived allele here", vcf_pos)
#     elif (haploid_gt != anc_nt) and (haploid_gt != der_nt) and (haploid_gt in nt_set):
#         print(f"messed up nt here", vcf_pos)
#     else:
#         pass  # nothing

#############################################################################
#############################################################################
#############################################################################

# Comparing recent individuals' vcfs against the HETPOS file created above

#Throw an error unless you have at least two arguments to command line
if len(sys.argv)<6:
    sys.exit('Input Error: The Anchor continuity test PART2 requires a command line input of format: python script.py the_chr anchor_ind recent_ind path_to_anchor_het path_to_vcf output_directory')

arg_list=sys.argv
the_chr=arg_list[1]
anchorind=arg_list[2]
recentind=arg_list[3]
hetpath=arg_list[4]
vcf_path=arg_list[5]
outdir=arg_list[6]
#print(the_chr)
#print(recentind)
#input()
###########################################################################################
###########################################################################################


NUCL=['A','C','G','T']
nt_set=set(NUCL)


#########################################
#########################################
#########################################

count_dict={the_chr: {}} # Initialize dict

#Create some win_steps to see effect of continuity across and within chromosomes, and also to enable wbj later
win_start=0
win_step=5000000    #(corresponds to about 5 cM)
win_end=win_start+win_step

print(the_chr,win_start,win_end)    #prints to slurm outfile
count_dict[the_chr].update({(win_start,win_end): [0,0,0,0,0,0]}) # In each window total number matching ancestral, derived, and other

# Begin reading HET and VCF file, ensure lines in sync
het_file_path = hetpath

with open(het_file_path, 'r') as het_file:
    with gzip.open(vcf_path, 'rt', encoding='utf-8') as RecentInd:

        # ---- Parse VCF header; find recent individual column ----
        line = RecentInd.readline()
        while line.startswith("##"):
            line = RecentInd.readline()

        # line is now the #CHROM header
        header_fields = line.strip().split('\t')
        sample_names = header_fields[9:]

        if recentind not in sample_names:
            sys.exit(f"ERROR: Recent individual '{recentind}' not found in VCF header.")

        recent_idx = 9 + sample_names.index(recentind)

        # Start reading data lines
        l = RecentInd.readline()
        het_l = het_file.readline()

        while l and het_l:
            vcf_data = l.strip().split()
            het_d = het_l.strip().split()
            vcf_pos = vcf_data[1]
            het_pos = het_d[0]

            # sync positions between VCF and HETPOS
            while vcf_pos != het_pos:
                if int(vcf_pos) < int(het_pos):
                    l = RecentInd.readline()
                else:
                    het_l = het_file.readline()

                if not (l and het_l):
                    break

                vcf_data = l.strip().split()
                vcf_pos = vcf_data[1]
                het_d = het_l.strip().split()
                het_pos = het_d[0]

            if not (l and het_l and het_pos == vcf_pos):
                break

            het_anc = het_d[1]
            het_der = het_d[2]
            vcf_ref = vcf_data[3]
            vcf_alt = vcf_data[4]
            at_pos = int(het_pos)

            while at_pos > win_end:  # Update dict with new window as key with empty value
                win_start += win_step
                win_end += win_step
                print(the_chr, win_start, win_end)
                count_dict[the_chr].update({(win_start, win_end): [0, 0, 0, 0, 0, 0]})

            qual = vcf_data[5]
            if qual != '.':
                flag = vcf_data[6]
                flag_list = ['FAIL', 'FAIL1', 'FAIL2', 'FAIL3']
                if flag not in flag_list:
                    format_str = vcf_data[8]
                    sample_str = vcf_data[recent_idx]

                    if not check_if_missingness(format_str, sample_str):
                        if check_if_biallelic(het_anc, het_der, vcf_ref, vcf_alt):
                            coverage, haploid_genotype = get_haploid_genotype(
                                vcf_ref, vcf_alt, format_str, sample_str
                            )
                            # If we failed to get a usable genotype, skip
                            if coverage > 0 and haploid_genotype in nt_set:
                                match_count = orient_and_get_count_haploid(
                                    haploid_genotype, het_anc, het_der
                                )
                                #check_der(haploid_genotype, het_anc, het_der)
                                count_dict[the_chr][(win_start, win_end)] = [
                                    count_dict[the_chr][(win_start, win_end)][i] + match_count[i]
                                    for i in range(6)
                                ]

            l = RecentInd.readline()
            het_l = het_file.readline()

with open(outdir+'/'+the_chr+'_'+anchorind+'_vs_'+recentind+'.txt','w') as outf:
    for a_tuple in sorted(count_dict[the_chr].keys()):
        out_str=str(a_tuple[0])+','+str(a_tuple[1])+'\t'+make_out_str(count_dict[the_chr][a_tuple])
        outf.write(out_str+'\n')



