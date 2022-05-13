import os,sys
import argparse
import io
import json
import copy
import gzip
import shlex
import shutil
import csv
import vcf
from collections import namedtuple
import subprocess as proc


def complement(a1, a2):
    c = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
    return (c[a1], c[a2])


def genotype(marker_id, a1, a2, ref, alt):
    # sanity check
    if len(ref) == 1 and len(alt) == 1 and (a1 == 'I' or a1 == 'D'):
        sys.exit('marker {} is a SNP but genotype is a INDEL'.format(marker_id))
    if (len(ref) > 1 or len(alt) > 1) and not (a1 == 'I' or a1 == 'D'):
        sys.exit('marker {} is a INDEL but genotype is a SNP'.format(marker_id))

    gt = None
    if len(ref)>len(alt):  # DEL
        if a1=='I' and a2=='I':
            gt = "0/0"
        elif a1=='D' and a2=='D':
            gt = "1/1"
        elif (a1=='I' and a2=='D') or (a1=='D' and a2=='I'):
            gt = "0/1"
        else:
            sys.exit('Marker Error DEL {}: {} {} {} {}'.format(marker_id, a1, a2, ref, alt))
    elif len(ref)<len(alt):  # INS
        if a1=='D' and a2=='D':
            gt = "0/0"
        elif a1=='I' and a2=='I':
            gt = "1/1"
        elif (a1=='I' and a2=='D') or (a1=='D' and a2=='I'):
            gt = "0/1"
        else:
            sys.exit('Marker Error INS {}: {} {} {} {}'.format(marker_id, a1, a2, ref, alt))
    elif len(ref)==1 and len(alt)==1:  # SNP
        if a1==ref and a2==ref:
            gt = "0/0"
        elif a1==alt and a2==alt:
            gt = "1/1"
        elif (a1==ref and a2==alt) or (a1==alt and a2==ref):
            gt = "0/1"
        else: # try complement alleles
            a3, a4 = complement(a1, a2)
            if a3 == ref and a4 == ref:
                gt = "0/0"
            elif a3 == alt and a4 == alt:
                gt = "1/1"
            elif (a3 == ref and a4 == alt) or (a3 == alt and a4 == ref):
                gt = "0/1"
            else:
                sys.exit('Marker Error SNP {}: {} {} {} {}'.format(marker_id, a1, a2, ref, alt))
    else:
        sys.exit('Marker Error REF {}: {} {} {} {}'.format(marker_id, a1, a2, ref, alt))

    if gt is None:
        sys.exit('Marker Error NoMatch {}: {} {} {} {}'.format(marker_id, a1, a2, ref, alt))
    else:
        return gt


def report2vcf(args):
    # read map .gz file
    markers = {}
    with gzip.open(args.gda_marker_path, 'rt') as f:
        for line in f:
            x = line.strip().split()
            markers[x[3]] = x

    # get sample_id from args.genotype_path file
    with open(args.genotype_path) as f:
        for i in range(11):
            line = f.readline()
        for line in f:
            sample_id = line.strip().split()[1]
            break

    # read vcf_b37 header from a template file
    vcf_header = vcf.Reader(filename=args.vcf_header_path)
    vcf_header.samples = [sample_id]

    # write vcf output file: b37
    root_name = '.'.join(os.path.basename(args.genotype_path).split('.')[:-1])
    vcf_out_b37 = os.path.join(args.working_dir, root_name+'.vcf')
    vcf_writer = vcf.Writer(open(vcf_out_b37, 'w'), vcf_header)

    # read genotypes from array report file, and write records to vcf output file
    with open(args.genotype_path) as f:
        for i in range(11):
            line = f.readline()
        for line in f:
            a = line.strip().split()
            marker_id = a[0]
            a1, a2 = a[2], a[3]
            if marker_id not in markers or a1=='-' or a2=='-':
                continue
            m = markers[marker_id]
            chrom, pos, snp_id, ref, alt = m[0], int(m[2]), m[4], m[5], m[6]

            # alt objects: list of objects
            alt_obj = [vcf.model._Substitution(alt)]

            # FORMAT string, and sample format object: named tuple
            samp_fmt_str = 'GT'
            samp_fmt = vcf.model.make_calldata_tuple(samp_fmt_str)

            # dict: {sample_1: 0, sample_2: 1, sample_3: 2, ...}
            sample_indexes = {sample_id: 0}

            # create vcf record object
            record = vcf.model._Record(chrom, pos, snp_id, ref, alt_obj, '', ['PASS'], [], samp_fmt_str, sample_indexes)

            # genotype str: "0/0", "0/1", or "1/1"
            gt = genotype(marker_id, a1, a2, ref, alt)

            # vcf.samples: list of _Call objects
            sampdat = gt.split(':')
            call = vcf.model._Call(record, sample_id, samp_fmt(*sampdat))
            record.samples.append(call)

            # write record to vcf file
            vcf_writer.write_record(record)
    vcf_writer.close()

    # picard lift over Vcf from b37 to hg38
    vcf_out_hg38 = vcf_out_b37.rstrip('.vcf') + '.hg38.vcf'
    reject_file = vcf_out_b37.rstrip('.vcf') + '.reject.vcf'
    cmd_liftover_vcf  = 'java -jar {} LiftoverVcf I={} O={} '.format(args.picard_jar, vcf_out_b37, vcf_out_hg38)
    cmd_liftover_vcf += 'C={} R={} REJECT={}'.format(args.chain_path, args.ref_path, reject_file)
    proc.run(cmd_liftover_vcf, shell=True)

    # annotate vcf with snp_id from SNPweight file
    vcf_out = vcf_out_b37.rstrip('.vcf') + '.hg38.vcf.gz'
    cmd_vcf_annot = "bcftools annotate -c CHROM,FROM,TO,ID -a {0} -Oz -o {1} {2}".format(args.snpwt_bed_path, vcf_out, vcf_out_hg38)
    proc.run(cmd_vcf_annot, shell=True)

    return vcf_out


def report_to_plink(args):
    # convert genotype report to plink format
    base_name = '.'.join(os.path.basename(args.genotype_path).split('.')[:-1])
    plink_base_path = os.path.join(args.working_dir, base_name)

    print('genotype_path: {}'.format(args.genotype_path))
    print('plink_base_path: {}'.format(plink_base_path))

    # create root_path.fam
    cmd_awk1 = '''awk 'NR<12{{next}} !seen[$2]++{{print 0,$2,0,0,0,0}}' {} > {}'''.format(args.genotype_path, plink_base_path+'.fam')
    proc.run(cmd_awk1, shell=True)

    # create root_path.lgen
    cmd_awk2 = '''awk -v OFS='\t' 'NR==FNR{{h[$4]=$5;next}} (FNR>11&&($1 in h)){{$3=($3=="-"?0:$3); $4=($4=="-"?0:$4); \
                  print 0,$2,h[$1],$3,$4}}' {0} {1} > {2}'''.format(args.gda_marker_path, args.genotype_path, plink_base_path+'.lgen')
    proc.run(cmd_awk2, shell=True)

    # create root_path.map
    plink_map = os.path.join(args.working_dir, plink_base_path+'.map')
    if not os.path.islink(plink_map):
        os.symlink(args.gda_map_path, plink_map)

    # plink recode
    plink_out_path = plink_base_path + '_out'
    cmd_plink1 = "plink --lfile {0} --recode --out {1}".format(plink_base_path, plink_out_path)
    proc.run(cmd_plink1, shell=True)

    # exclude markers which has string length greater than 39. Eigensoft/convertf does not work with long marker id
    plink_filt_path = plink_base_path + '_filt'
    cmd_plink2 = "plink --file {0} --exclude {1} --chr 1-22 --geno 0.0 --recode --out {2}".format(plink_out_path, args.gda_long_snpid, plink_filt_path)
    proc.run(cmd_plink2, shell=True)
    # plink --bfile GDA_b37 --exclude snp.exclude --chr 1-22 --geno 0.1 --recode --out GDA_b37

    return plink_filt_path


def plink_to_eigenstrat(plink_base_path):
    plink2eigenstrat_par = plink_base_path+'.plink2eigenstart.par'
    with open(plink2eigenstrat_par, 'w') as f:
        f.write('genotypename:\t{}.ped\n'.format(plink_base_path))
        f.write('snpname:\t{}.map\n'.format(plink_base_path))
        f.write('indivname:\t{}.ped\n'.format(plink_base_path))
        f.write('outputformat:\tEIGENSTRAT\n')
        f.write('genotypeoutname:\t{}.geno\n'.format(plink_base_path))
        f.write('snpoutname:\t{}.snp\n'.format(plink_base_path))
        f.write('indivoutname:\t{}.ind\n'.format(plink_base_path))
        f.write('familynames:\tNO\n')
    cmd_eigenstrat = "convertf -p {}".format(plink2eigenstrat_par)
    proc.run(cmd_eigenstrat, shell=True)

    geno_base_path = plink_base_path
    return geno_base_path


def call_variant(bam_path, ref_file, annot_file, working_dir):
    bam_name = os.path.basename(bam_path)
    base_name = '.'.join(bam_name.split('.')[:-1])
    vcf_name = '.'.join([base_name, 'vcf'])
    vcf_path = os.path.join(working_dir, vcf_name)

    cmd_vcf  = "samtools mpileup -q 30 -Q 20 -v -f {0} -l {1} {2} | ".format(ref_file, annot_file, bam_path)
    cmd_vcf += "bcftools call -c -Ov | bcftools filter -e 'ALT=\".\"' | "
    cmd_vcf += "bcftools annotate -c CHROM,FROM,TO,ID -a {0} > {1}".format(annot_file, vcf_path)

    # run cmd_vcf
    if not os.path.isfile(vcf_path):
        proc.run(cmd_vcf, shell=True)

    if os.path.isfile(vcf_path):
        return vcf_path
    else:
        sys.exit('Function call_variant failed')


def subset_variant(vcf_path, snpwt_ref, annot_file, working_dir):
    # return subset of vcf_path: vcf2_path
    sample_key = os.path.basename(vcf_path).split('.')[0]
    vcf2_path = os.path.join(working_dir, '.'.join([sample_key, snpwt_ref, 'vcf']))
    cmd_vcf2  = "bcftools view {0} -v snps -f PASS -T {1} -Ov | ".format(vcf_path, annot_file)
    cmd_vcf2 += "bcftools annotate -c CHROM,FROM,TO,ID -a {0} > {1}".format(annot_file, vcf2_path)

    # run cmd_vcf2
    if not os.path.isfile(vcf2_path):
        proc.run(cmd_vcf2, shell=True)

    if os.path.isfile(vcf2_path):
        return vcf2_path
    else:
        sys.exit('Function subset_variant failed')


def vcf_to_eigenstrat(vcf2eigenstart_py, vcf_path):
    geno_base_path = '.'.join(vcf_path.split('.')[:-1])
    cmd_geno  = "python2 {0} -v {1} -o {2}".format(vcf2eigenstart_py, vcf_path, geno_base_path)

    # run cmd_geno
    if not os.path.isfile(geno_base_path+'.geno'):
        proc.run(cmd_geno, shell=True)

    if os.path.isfile(geno_base_path+'.geno'):
        return geno_base_path
    else:
        sys.exit('Function vcf_to_eigenstrat failed')


def infer_ancestry(ancestry_py, geno_base_path, snpwt_path):
    par_file = geno_base_path+'.par'
    with open(par_file,'w') as f:
        f.write('geno:\t{}.geno\n'.format(geno_base_path))
        f.write('snp:\t{}.snp\n'.format(geno_base_path))
        f.write('ind:\t{}.ind\n'.format(geno_base_path))
        f.write('snpwt:\t{}\n'.format(snpwt_path))
        f.write('predpcoutput:\t{}.predpc\n'.format(geno_base_path))

    cmd_snpwt = "python2 {0} --par {1}".format(ancestry_py, par_file)

    # run cmd_snpwt
    if not os.path.isfile(geno_base_path+'.predpc'):
        proc.run(cmd_snpwt, shell=True)

    if os.path.isfile(geno_base_path+'.predpc'):
        return geno_base_path+'.predpc'
    else:
        sys.exit('Function infer_ancestry failed')


def main(args):
    # input: Illumina GDA array report file
    if args.genotype_path is not None:
        plink_base_path = report_to_plink(args)
        # convert plink to eigenstrat format
        geno_base_path = plink_to_eigenstrat(plink_base_path)
        # infer ancestry
        predpc_path = infer_ancestry(args.ancestry_exe, geno_base_path, args.snpwt_path)
        print(predpc_path)

    # input: bam_path, and it needs to call variants
    if args.bam_path is not None:
        bam_name = os.path.basename(bam_path)
        base_name = '.'.join(bam_name.split('.')[:-1])
        vcf_name = '.'.join([base_name, 'vcf'])
        vcf_path = os.path.join(working_dir, vcf_name)
        if not os.path.isfile(vcf_path):  # call variants once
            vcf_path = call_variant(args.bam_path, args.ref_path, args.snpwt_bed_path, args.working_dir)
        if args.snpwt_path.endswith('.AS'):
            # subset for snpwt.AS: ASIAN reference panel, order: SAS EAS AFR EUR
            snpwt_ref = 'AS'
            vcf2_path = subset_variant(vcf_path, snpwt_ref, args.snpwt_bed_path, args.working_dir)
            # convert vcf to eigenstrat
            geno_base_path = vcf_to_eigenstrat(args.vcf2eigenstrat_exe, vcf2_path)
            # infer ancestry
            predpc_path = infer_ancestry(args.ancestry_exe, geno_base_path, args.snpwt_path)
        else:
            # subset for snpwt.NA: Native American reference panel, order: AFR EUR EAS NAT
            snpwt_ref = 'NA'
            vcf2_path = subset_variant(vcf_path, snpwt_ref, args.snpwt_bed_path, args.working_dir)
            # convert vcf to eigenstrat
            geno_base_path = vcf_to_eigenstrat(args.vcf2eigenstrat_exe, vcf2_path)
            # infer ancestry
            predpc_path = infer_ancestry(args.ancestry_exe, geno_base_path, args.snpwt_path)
        print(predpc_path)

    # input: vcf_path from phoenix pipeline, and it needs to extract subset/vcf: snpwt_bed
    if args.vcf_path is not None:
        if args.snpwt_path.endswith('.AS'):
            # subset for snpwt.AS: ASIAN reference panel, order: SAS EAS AFR EUR
            snpwt_ref = 'AS'
            vcf2_path = subset_variant(args.vcf_path, snpwt_ref, args.snpwt_bed_path, args.working_dir)
            # convert vcf to eigenstrat
            geno_base_path = vcf_to_eigenstrat(args.vcf2eigenstrat_exe, vcf2_path)
            # infer ancestry
            predpc_path = infer_ancestry(args.ancestry_exe, geno_base_path, args.snpwt_path)
        else:
            # subset for snpwt.NA: Native American reference panel, order: AFR EUR EAS NAT
            snpwt_ref = 'NA'
            vcf2_path = subset_variant(args.vcf_path, snpwt_ref, args.snpwt_bed_path, args.working_dir)
            # convert vcf to eigenstrat
            geno_base_path = vcf_to_eigenstrat(args.vcf2eigenstrat_exe, vcf2_path)
            # infer ancestry
            predpc_path = infer_ancestry(args.ancestry_exe, geno_base_path, args.snpwt_path)
        print(predpc_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="snpwt_ancestry", description="infer ancestry with germline bam file for WES/WGS.")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--genotype_path', help="An array genotype report file from Illumina Genome Studio")
    group.add_argument('--bam_path', help="A alignment file in bam or cram format")
    group.add_argument('--vcf_path', help="A variant call file in vcf format")
    parser.add_argument('--ref_path', default="/home/tgenref/homo_sapiens/grch38_hg38/hg38tgen/genome_reference/GRCh38tgen_decoy_alts_hla.fa")
    parser.add_argument('--out_dir', default="/home/gzhang/scratch/rodriguez/out", help="Output directory")
    parser.add_argument('--working_dir', default="/home/gzhang/scratch/rodriguez/scratch", help="Scratch directory")
    parser.add_argument('--cohort', default="COHORIEN", help="Cohort name")
    parser.add_argument('--ancestry_exe', default="/home/gzhang/compute/github/Ancestry-SNPweights/SNPweights2.1/inferancestry.py", help="infer ancestry executable")
    parser.add_argument('--convertf_exe', default="/home/gzhang/compute/github/Ancestry-SNPweights/eigensoft/convertf", help="eigensoft convertf executable")
    parser.add_argument('--vcf2eigenstrat_exe', default="/home/gzhang/compute/github/Ancestry-SNPweights/gdc/vcf2eigenstrat.py", help="convert vcf to eigenstrat format")
    parser.add_argument('--snpwt_path', default="/home/gzhang/compute/github/Ancestry-SNPweights/data/snpwt.NA", help="snp weights file for Native Aemrican and more")
    parser.add_argument('--snpwt_bed_path', default="/home/gzhang/compute/github/Ancestry-SNPweights/data/snpwt_NA.hg38.bed.gz", help="snp weights bed file for NA in hg38")
    parser.add_argument('--gda_marker_path', default="/home/gzhang/compute/github/Ancestry-SNPweights/data/GDA_GRCh37.dbSNP151.genotype.bed", help="marker file in GRCh37 for GDA array")
    parser.add_argument('--gda_map_path', default="/home/gzhang/compute/github/Ancestry-SNPweights/data/GDA_GRCh37.map", help="plink map file in GRCh37 for GDA array")
    parser.add_argument('--gda_long_snpid', default="/home/gzhang/compute/github/Ancestry-SNPweights/data/GDA_GRCh37.long39.snpid", help="list of SNP IDs with string length > 39")
    parser.add_argument('--vcf_header_path', default="/home/gzhang/compute/github/Ancestry-SNPweights/data/header.b37.vcf", help="vcf header file in GRCh37 for GDA array")
    parser.add_argument('--samtools_exe', default="/home/gzhang/bin/samtools", help="samtools executable")
    parser.add_argument('--bcftools_exe', default="/home/gzhang/bin/bcftools", help="bcftools executable")
    parser.add_argument('--plink_exe', default="/home/gzhang/bin/plink", help="plink executable")
    parser.add_argument('--picard_jar', default="/home/gzhang/bin/picard.jar", help="picard jar binary")
    parser.add_argument('--chain_path', default="/home/gzhang/compute/github/Ancestry-SNPweights/data/b37ToHg38.over.chain", help="liftover chain file")

    args, remaining_argv = parser.parse_known_args()

    # add executables to env PATH
    if os.path.dirname(args.samtools_exe) not in os.environ['PATH']:
        os.environ['PATH'] += ":"+os.path.dirname(args.samtools_exe)
    if os.path.dirname(args.bcftools_exe) not in os.environ['PATH']:
        os.environ['PATH'] += ":"+os.path.dirname(args.bcftools_exe)

    main(args)
