import os,sys
import argparse
import shutil
import csv
import subprocess as proc


def parse_output(file):
    with open(file) as f:
        for line in f:
            a = line.strip().split()
            #return [a[6],a[7],a[8],a[9]]
            # [:-num_pop]
            return a[-4:]

def main(args):
    bam_path = args.bam_path
    if not os.path.isfile(bam_path):
        sys.exit('bam/cram file not exist: {}'.format(bam_path))
    bam_name = os.path.basename(bam_path)

    # infer ancestry using 4 Population module 'snpwt.NA' including Native American: 'NA'
    predpc_path = os.path.join(args.working_dir, '.'.join([bam_name.split('.')[0], 'NA', 'predpc']))
    cmd_infer_ancestry = 'python {} --bam_path {} --working_dir {}'.format(args.ancestry_script, bam_path, args.working_dir)
    # run cmd_infer_ancestry
    #if not os.path.isfile(predpc_path) or os.path.getsize(predpc_path) == 0:
    proc.run(cmd_infer_ancestry, shell=True)
    if os.path.isfile(predpc_path) and os.path.getsize(predpc_path) > 0:
        # predpc_na = ['AFR', 'EUR', 'EAS', 'NAT']
        predpc_na = parse_output(predpc_path)
    else:
        sys.exit('infer_ancestry {} failed'.format('NA module'))

    # infer ancestry using 4 population module 'snpwt.AS' including ASIAN: 'AS'
    predpc_path = os.path.join(args.working_dir, '.'.join([bam_name.split('.')[0], 'AS', 'predpc']))
    cmd_infer_ancestry = 'python {} --bam_path {} --working_dir {} --snpwt_path {} --snpwt_bed_path {}'.format(args.ancestry_script,
                                      bam_path, args.working_dir, args.snpwt_as_path, args.snpwt_as_bed_path)
    # run cmd_infer_ancestry
    #if not os.path.isfile(predpc_path) or os.path.getsize(predpc_path) == 0:
    proc.run(cmd_infer_ancestry, shell=True)
    if os.path.isfile(predpc_path) and os.path.getsize(predpc_path) > 0:
        # predpc_as = ['SAS', 'EAS', 'AFR', 'EUR']
        predpc_as = parse_output(predpc_path)
    else:
        sys.exit('infer_ancestry {} failed'.format('AS module'))

    # write out csv file
    out_file = '.'.join(bam_name.split('.')[:1]+['predpc.csv'])
    f_out = open(out_file, 'w')
    # header line
    header = ['AFR', 'EUR', 'SAS', 'EAS', 'NAT', 'Ancestry']
    print('\t'.join(header))
    f_out.write(','.join(header) + '\n')

    anc = {}
    # 4-POP 'AS' parsed order: SAS EAS AFR EUR
    anc['SAS'], anc['EAS'], anc['AFR'], anc['EUR'] = [round(float(x),3) for x in predpc_as]
    anc['NAT'] = '.'
    if anc['AFR'] >= 0.8:
        ancestry = 'AFR'
    elif anc['EUR'] >= 0.8:
        ancestry = 'EUR'
    elif anc['SAS'] >= 0.8:
        ancestry = 'SAS'
    elif anc['EAS'] >= 0.8:
        ancestry = 'EAS'
    else:
        # 4-POP 'NA' parsed order: AFR EUR EAS NAT
        anc = {}
        anc['AFR'], anc['EUR'], anc['EAS'], anc['NAT'] = [round(float(x), 3) for x in predpc_na]
        if anc['AFR'] >= 0.8:
            ancestry = 'AFR'
        elif anc['EUR'] >= 0.8:
            ancestry = 'EUR'
        elif anc['EAS'] >= 0.8:
            ancestry = 'EAS'
        elif anc['NAT'] >= 0.8:
            ancestry = 'NAT'
        else:
            anc_sorted = sorted(anc, key=anc.get, reverse=True)
            ancestry = 'Mix({}+{})'.format(anc_sorted[0], anc_sorted[1])
        anc['SAS'] = '.'

    # 5-POP output order: 'AFR', 'EUR', 'SAS', 'EAS', 'NAT'
    pop5 = [anc['AFR'], anc['EUR'], anc['SAS'], anc['EAS'], anc['NAT'], ancestry]
    line = [str(x) for x in pop5]
    print('\t'.join(line))
    f_out.write(','.join(line) + '\n')
    f_out.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="infer_ancestry_bam.py", description="infer ancestry using a bam/cram file.")
    parser.add_argument('-b', '--bam_path', required=True, help="A bam file path")
    parser.add_argument('-o', '--out_dir', default="/home/gzhang/scratch/rodriguez/out", help="Output directory")
    parser.add_argument('-a', '--ancestry_script', default="/home/gzhang/compute/github/Ancestry-SNPweights/snpwt_ancestry.py",
                        help="ancestry estimate script")
    parser.add_argument('-w', '--working_dir',  default="/home/gzhang/scratch/ancestry", help="Scratch directory")
    parser.add_argument('--snpwt_as_path', default="/home/gzhang/compute/github/Ancestry-SNPweights/data/snpwt.AS",
                        help="snp weights file for AS module")
    parser.add_argument('--snpwt_as_bed_path', default="/home/gzhang/compute/github/Ancestry-SNPweights/data/snpwt_AS.hg38.bed.gz",
                        help="snp weights bed file for AS in hg38")

    args, remaining_argv = parser.parse_known_args()

    main(args)
