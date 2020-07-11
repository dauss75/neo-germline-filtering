#!/isilon/prod2/bcbio/anaconda/bin/python
#/usr/bin/python
#/isilon/prod2/bcbio/anaconda/bin/python --> version 3

'''
nohup python /isilon/RnD/tools/custom_script/purecn/bin/run_pureCN_solid_commands_v0.2.py
-t ../../../data/bam/tmb_msi_T_bam
-n ../../../data/bam/normal
-v ../../../data/vcf/T_new_mutect2
-p /isilon/RnD/tools/custom_script/purecn/pon/qiaseq323/pon_qiaseq323_14normals_min1.vcf.gz
-o 0630_mutect_new_T &
'''

import os
import os.path
import subprocess
import argparse
from os import listdir
from os.path import isfile, join
from joblib import Parallel, delayed
import pyranges as pr   # header needs to match as Chromosome Start End
import pandas as pd
import glob
from datetime import date

NCPU = 20
CHUNK_SIZE = 2**20


def parallel_run(output_dir, sn, rscript, purecn_exe, bam_files_dir, sn_bam_dict, intervals_output_file):
    """
    creates files:
    _coverage_loess.png
    _coverage_loess.txt
    _coverage_loess_qc.txt
    _coverage.txt

    """
    coverage_command = [rscript, purecn_exe + '/Coverage.R',
                       '--force',
                       '--keepduplicates',
                       '--outdir', output_dir + '/' + sn,
                       '--bam', '/'.join([bam_files_dir, sn_bam_dict[sn]]),
                       '--intervals', intervals_output_file]
    print(' '.join(coverage_command))

    p = subprocess.Popen(" ".join(coverage_command), stdout=subprocess.PIPE, shell=True)
    out = p.communicate()
    return_code = p.wait()

# def parallel_run_pureCN(output_dir, sn, rscript, purecn_exe, tumorcov, normaldb_file, mapping_bias, genome, tumorvcf, intervals_output_file):
#
#     pureCN_command = [rscript, purecn_exe + '/PureCN.R',
#                        '--force',
#                        '--out', output_dir + '/' + sn + '/',
#                        '--sampleid', sn,
#                        '--tumor', tumorcov,
#                        '--vcf', tumorvcf,
#                        '--mappingbiasfile', mapping_bias,
#                        '--normaldb', normaldb_file,
#                        '--intervals', intervals_output_file,
#                        '--force --outvcf --postoptimize --seed 123',
#                        '--funsegmentation PSCBS',
#                        '--genome', genome]
#
#     print(' '.join(pureCN_command))
#
#     p = subprocess.Popen(" ".join(pureCN_command), stdout=subprocess.PIPE, shell=True)
#     out = p.communicate()
#     return_code = p.wait()

def parallel_run_pureCN_SGZ(output_dir, sn, rscript, purecn_exe, tumorcov, normaldb_file, mapping_bias, genome, tumorvcf, intervals_output_file):

    pureCN_command = [rscript, purecn_exe + '/PureCN.R',
                       '--force',
                       '--out', output_dir + '/' + sn + '/',
                       '--sampleid', sn,
                       '--tumor', tumorcov,
                       '--vcf', tumorvcf,
                       '--mappingbiasfile', mapping_bias,
                       '--normaldb', normaldb_file,
                       '--intervals', intervals_output_file,
                       '--force --outvcf --postoptimize --seed 123',
                       '--funsegmentation PSCBS',
                       '--genome', genome]

    print(' '.join(pureCN_command))

    p = subprocess.Popen(" ".join(pureCN_command), stdout=subprocess.PIPE, shell=True)
    out = p.communicate()
    return_code = p.wait()

    sgz_command = [rscript, purecn_exe + '/TumorOnlySGZ_v01.R',
                       '--rds', output_dir + '/' + sn + '/' + sn +'.rds',
                       '--out', output_dir + '/' + sn + '/' + sn,
                       '--vcf', tumorvcf]
    p = subprocess.Popen(" ".join(sgz_command), stdout=subprocess.PIPE, shell=True)
    out = p.communicate()
    return_code = p.wait()

def __main__():

    parser = argparse.ArgumentParser(description='Parse Carlsbad Archer VCs and output into single file')
    parser.add_argument('-t', '--tumor_bam_dir',  help="Directory containing tumor bam files", required=True)
    parser.add_argument('-n', '--normal_bam_dir',  help="Directory containing normal bam files", required=False)
    parser.add_argument('-v', '--tumor_vcf_dir',  help="Directory containing tumor vcf files", required=True)
    parser.add_argument('-p', '--normal_panel',  help="normal vcf file", required=False)
    parser.add_argument('-o', '--output_dir',  help="Directory pipe outputs", required=True)

    args = parser.parse_args()

    rscript = '/isilon/RnD/tools/bcbio_dev_081519/bin/Rscript'

    input_bed_file = '/isilon/RnD/tools/custom_script/cnvkit/bed/QIAseq323.CDHS-24104Z-14204.refseq-anno.roi.exons.nolowcovrois_genes_only.bed' #Note, had to remove the header line to make work
    ref_fa = '/isilon/RnD/bcbio-tool-data-dev/genomes/Hsapiens/GRCh37/seq/GRCh37.fa'
    genome = 'hg19'
    mappability_file = '/isilon/R_and_D/user_folders/sjung/project/msi/pureCN/reference_files/GRCh37_mappability/GRCh37_100.bw'
    intervals_output_file =args.output_dir +'/baits_' + genome + '_intervals.txt'
    purecn_exe = '/isilon/RnD/bcbio-tool-data-dev/anaconda/lib/R/library/PureCN/extdata'

    if not os.path.exists(args.output_dir):
        os.mkdir(args.output_dir)

    print('--------------------------------\n| generate interval |\n--------------------------------\n')
    prepare_command = [rscript, purecn_exe + '/IntervalFile.R',
                        '--force',
                        '--infile', input_bed_file,
                        '--fasta', ref_fa,
                        '--outfile', intervals_output_file,
                        '--genome', genome,
                        '--export', args.output_dir + '/baits_optimized' + genome + '.bed',
                        '--mappability', mappability_file]

    cmd1 = " ".join(prepare_command)
    print(cmd1)
    p = subprocess.Popen(cmd1, stdout=subprocess.PIPE, shell=True)
    out = p.communicate()
    return_code = p.wait()

    # 3. PureCN Create VCFs
    # VCFs were generated by ensemble
    EXT = (".vcf.gz",".vcf")
    if args.tumor_vcf_dir:
        tumorvcffiles = [f for f in sorted(listdir(args.tumor_vcf_dir)) if isfile(join(args.tumor_vcf_dir, f)) and f.endswith(EXT)]

    print('--------------------------------\n| generate tumor coverage |\n--------------------------------\n')
    # 4. PureCN Coverage
    if args.tumor_bam_dir:
        tumorbamfiles = [f for f in sorted(listdir(args.tumor_bam_dir)) if isfile(join(args.tumor_bam_dir, f)) and f.endswith('.bam')]
        tumor_sample_list = []
        tumor_sn_bam_dict = {}
        tumor_coverage_list = []

        for bam in tumorbamfiles:
            sn = os.path.basename(bam).split('.')[0]
            tumor_sample_list.append(sn)
            tumor_sn_bam_dict[sn] = bam
            if not os.path.exists(args.output_dir + '/' + sn):
                os.makedirs(args.output_dir + '/' + sn)

            tumor_coverage_list.append(args.output_dir + '/' + sn + '/' + sn + '_coverage_loess.txt.gz')

        Parallel(n_jobs=NCPU)(delayed(parallel_run)(output_dir=args.output_dir, \
                              sn=tumor_sample_list[i], rscript=rscript, purecn_exe=purecn_exe, bam_files_dir=args.tumor_bam_dir, \
                              sn_bam_dict=tumor_sn_bam_dict, intervals_output_file=intervals_output_file) for i in range(0,len(tumor_sample_list)))

    print('--------------------------------\n| generate normal coverage |\n--------------------------------\n')
    if args.normal_bam_dir:
        normalbamfiles = [f for f in sorted(listdir(args.normal_bam_dir)) if isfile(join(args.normal_bam_dir, f)) and f.endswith('.bam')]
        normal_sample_list = []
        normal_sn_bam_dict = {}
        normal_coverage_list = []
        for bam in normalbamfiles:
            sn = os.path.basename(bam).split('.')[0]
            normal_sample_list.append(sn)
            normal_sn_bam_dict[sn] = bam
            if not os.path.exists(args.output_dir + '/' + sn):
                os.makedirs(args.output_dir + '/' + sn)

            normal_coverage_list.append(args.output_dir + '/' + sn + '/' + sn + '_coverage_loess.txt.gz')

        if normal_coverage_list:
            normal_coverage_list_file_name = args.output_dir + '/normal_coverage_list_file.txt'
            normal_coverage_list_file = open(normal_coverage_list_file_name, 'w')
            for i in normal_coverage_list:
                normal_coverage_list_file.write(i + '\n')
            normal_coverage_list_file.flush()
            normal_coverage_list_file.close()
        Parallel(n_jobs=NCPU)(delayed(parallel_run)(output_dir=args.output_dir, \
                              sn=normal_sample_list[i], rscript=rscript, purecn_exe=purecn_exe, bam_files_dir=args.normal_bam_dir, \
                              sn_bam_dict=normal_sn_bam_dict, intervals_output_file=intervals_output_file) for i in range(0,len(normal_sample_list)))
        print('--------------------------------\n| generate normalDB |\n--------------------------------\n')
        normal_panel = args.normal_panel    # normal vcf

        cmd2 = '{} {}/NormalDB.R --outdir {} --coveragefiles {}/normal_coverage_list_file.txt --force --genome hg19 --normal_panel {}'.format(rscript, purecn_exe, args.output_dir, args.output_dir, normal_panel)
        print('cmd2: ' + cmd2)
        p = subprocess.Popen(cmd2, stdout=subprocess.PIPE, shell=True)
        out = p.communicate()
        return_code = p.wait()

        normaldb_file = args.output_dir + '/' + 'normalDB_hg19.rds'
        mapping_bias =  args.output_dir + '/' + 'mapping_bias_hg19.rds'
    else:
        normaldb_file = '/isilon/RnD/tools/custom_script/purecn/mapping_bias/qiaseq323/pon_min1_qiaseq323_14normals_normalDB_hg19.rds'
        mapping_bias =  '/isilon/RnD/tools/custom_script/purecn/mapping_bias/qiaseq323/pon_min1_qiaseq323_14normals_mapping_bias_hg19.rds'

    # # 5. Run PureCN for Tumor Files
    tumorvcfs = []
    for v in tumorvcffiles:
        tumorvcfs.append(args.tumor_vcf_dir+'/'+ v)
    print('--------------------------------\n| run PureCN.R |\n--------------------------------\n')
    if len(tumorbamfiles) == len(tumorvcffiles):
        Parallel(n_jobs=NCPU)(delayed(parallel_run_pureCN_SGZ)(output_dir=args.output_dir, \
                          sn=os.path.basename(tumor_coverage_list[i]).split('_coverage')[0], rscript=rscript, purecn_exe=purecn_exe, \
                          tumorcov=tumor_coverage_list[i], normaldb_file=normaldb_file, mapping_bias=mapping_bias, genome=genome, tumorvcf=tumorvcfs[j], \
                          intervals_output_file=intervals_output_file) for i,j in zip(range(0,len(tumor_coverage_list)),range(0,len(tumorvcfs))))


        '''
        create a final report
        '''
        fields1 = ['Sampleid','chr', 'start','end','arm','C','M','type','seg.mean','M.flagged']
        fields2 = ['Sampleid','chr', 'start','end','ID','REF','ALT','ML.SOMATIC','POSTERIOR.SOMATIC','ML.LOH','log.ratio','depth','gene.symbol']
        outputF = '{}/{}_purecn_final_report.xlsx'.format(args.output_dir,args.output_dir)
        x = []

        flag=0
        print('--------------------------------\n| generate a final report |\n--------------------------------\n')
        with pd.ExcelWriter(outputF, mode='w') as writer:
            # lohFiles=sorted(glob.glob('{}/[NTP,MOL]*/*_loh.csv'.format(args.output_dir)))
            # genesFiles=sorted(glob.glob('{}/[NTP,MOL]*/*_variants.csv'.format(args.output_dir)))
            '''
            R&D only
            '''
            lohFiles=sorted(glob.glob('{}/[Acc,NTP,MOL]*/*_loh.csv'.format(args.output_dir)))
            genesFiles=sorted(glob.glob('{}/[Acc,NTP,MOL]*/*_variants.csv'.format(args.output_dir)))

            for lohF, geneF in zip(lohFiles,genesFiles):
                sn1 = '-'.join(os.path.basename(geneF).split('-')[0:3])
                try:
                    df1 = pd.read_csv(lohF, delimiter=',', comment='#', usecols=fields1, skip_blank_lines=True,encoding='utf-8')
                    df2 = pd.read_csv(geneF, delimiter=',', comment='#', usecols=fields2, skip_blank_lines=True,encoding='utf-8')

                    if df1.shape[0]!=0:
                        df1 = df1[(df1['C']!=2) | (df1['M']==0)]  # exclude CN=2 and yet keep whoe neutral arm
                        df1.sort_values(by=['chr'])

                        # print(df1)
                        x.append(df1)
                        if flag==0:
                            df1.to_excel(writer,sheet_name='master',index=False)
                            flag=flag+1
                        #df1.to_excel(writer, sheet_name=sn1,index=False)

                        # create PyRanges-objects from the dfs
                        # pr requires a certain format with column names Chromosome Start End
                        df1.rename(columns={'chr':'Chromosome','start':'Start','end':'End'},inplace=True)
                        df2.rename(columns={'chr':'Chromosome','start':'Start','end':'End'},inplace=True)

                        gr1, gr2 = pr.PyRanges(df1),pr.PyRanges(df2)

                        # intersect the two
                        gr = gr2.intersect(gr1)

                        gr.df.to_excel(writer, sheet_name=sn1, index=False)
                except:
                    print(sn1)
            # print(x)
            x1 = pd.concat(x)
            x1.to_excel(writer, sheet_name='master',index=False)

if __name__=="__main__": __main__()
