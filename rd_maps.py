#!/usr/bin/env python
import sys
import os
import glob
import subprocess as sub
import time
import csv
import argparse

"""
Re-implementation of:
https://github.com/rob123king/EMS_test_scripts/blob/master/Dragen_VCF_filtering2.pl

Helpful reading:
https://www.biostars.org/p/84951/
http://bedtools.readthedocs.io/en/latest/content/tools/coverage.html
http://bedtools.readthedocs.io/en/latest/content/overview.html
https://stackoverflow.com/questions/5466451/how-can-i-print-literal-curly-brace-characters-in-python-string-and-also-use-fo
https://stackoverflow.com/questions/6997430/looping-over-input-fields-as-array
https://stackoverflow.com/questions/32481877/what-is-nr-fnr-in-awk
"""

SBATCH = 'sbatch -p {} --mem {} -c {} --wrap="{}"'
SBATCH_ARRAY = 'sbatch -p {} --mem {} -c {} --wrap="{}" --array=1-{}'

VCFMERGE = 'source vcftools-0.1.13; vcf-merge {} | grep -v \'#\' | cut -f 1,2 | awk -v OFS=\'\\t\' \'{{ print \$1,\$2-1,\$2; }}\' > {}; touch {};'
MULTICOV_ARRAY = 'source bedtools-2.26.0; bed=\$(head -n \$SLURM_ARRAY_TASK_ID {} | tail -n 1); bedtools multicov -q 1 -p -bams {} -bed \$bed | awk -v OFS=\'\\t\' -v mins={} \'{{ c=0; for(i = 4; i <= NF; i++) {{ if (\$i > 0) c++; }}; if (c >= mins) print \$0; }}\' > \$bed.mincov; touch \$bed.done;'
#DIFILTER = 'source vcftools-0.1.13; gzip -dc {0} | awk -F\'\\t\' \'/^#/{{ print \$0; next; }}; NR==FNR{{c[\$1\$3]++; next; }}; c[\$1\$2] > 0\' {1} - | awk -v OFS=\'\\t\' \'/^#/{{print \$0; next;}} {{ split(\$10, data, \\":\\"); split(data[2], ad, \\",\\"); if (data[1] == \\"0/1\\" && (ad[1] > 0 || ad[2] > 0)) {{ if ( ad[2] / (ad[1] + ad[2]) > 0.1999 && ad[2] >= 5) print \$0; }} else if (data[1] == \\"1/1\\" && ad[2] >= 3) print \$0; }}\' | bgzip -c > {2} && tabix -p vcf {2}; touch {3};'
DIFILTER = 'source vcftools-0.1.13; gzip -dc {0} | awk -F\'\\t\' \'/^#/{{ print \$0; next; }}; NR==FNR{{c[\$1\$3]++; next; }}; c[\$1\$2] > 0\' {1} - | awk -v OFS=\'\\t\' \'/^#/{{print \$0; next;}} {{ split(\$10, data, \\":\\"); split(data[2], ad, \\",\\"); if (data[1] == \\"0/1\\" && (ad[1] > 0 || ad[2] > 0)) {{ if ( ad[2] / (ad[1] + ad[2]) > 0.1999 && ad[2] >= {4}) print \$0; }} else if (data[1] == \\"1/1\\" && ad[2] >= {5}) print \$0; }}\' | bgzip -c > {2} && tabix -p vcf {2}; touch {3};'
INTERSECT = 'source vcftools-0.1.13; vcf-isec -f -c {} {} > {}; touch {};'

"""
obsolete, potentially defunct:
MULTICOV = 'source bedtools-2.26.0; bedtools multicov -q 1 -p -bams {} -bed {} | awk -v OFS=\'\\t\' -v mins={} \'{{ c=0; for(i = 4; i <= NF; i++) {{ if (\$i > 0) c++; }}; if (c >= mins) print \$0; }}\' > {}; touch {};' 

VCFFILTER = 'gzip -dc {0} | grep \'^#\' | gzip > {1}; gzip -dc {0} | grep -v \'#\' | awk -F\'\\t\' \'NR==FNR{{c[\$1\$3]++;next}};c[\$1\$2] > 0\' {3} - | gzip >> {1}; touch {2};'

HHFILTER = 'gzip -dc {} | awk -v OFS=\'\\t\' \'{{ split(\$10, data, \\":\\"); split(data[2], ad, \\",\\"); if (data[1] == \\"0/1\\" && (ad[1] > 0 || ad[2] > 0)) {{ if ( ad[2] / (ad[1] + ad[2]) > 0.1999 && ad[2] >= 5) print \$0; }} else if (data[1] == \\"1/1\\" && ad[2] >= 3) print \$0; }}\' | gzip > {}; touch {};'
"""

def splitInputFile(inputfile):
    nlines = sum(1 for line in open(inputfile))        
    with open(inputfile) as _in:
        part = -1
        outf = None
        for i, line in enumerate(_in):
            current = i // 20000
            if current != part:
                part = current
                if outf is not None:
                    outf.close()
                f = inputfile + '.part{}'.format(current)
                yield f
                outf = open(f, 'w')
            print(line.strip(), file=outf)
        if outf is not None:
            outf.close()
    pass


def slurmJob(cmd, sentinels, p='ei-medium', mem='8GB', c=1, arraylen=0, timelimit_s=24*3600, wait_s=180):
    if arraylen:
        sbatch = SBATCH_ARRAY.format(p, mem, c, cmd, arraylen)
    else:
        sbatch = SBATCH.format(p, mem, c, cmd)
    pr = sub.Popen(sbatch, stdout=sub.PIPE, stderr=sub.PIPE, stdin=sub.PIPE, shell=True)        
    err, out = pr.communicate()
    tstart = time.time()
    if arraylen:
        missing = list(sentinels)
        while True:
            missing = list(f for f in missing if not os.path.exists(f))
            if not missing or time.time() - tstart > timelimit_s:
                break
            time.sleep(wait_s)
    else:
        while not os.path.exists(sentinels[0]) and time.time() - tstart < timelimit_s:
            time.sleep(wait_s)

def slurmMultiJob(cmds, sentinels, p='ei-medium', mem='8GB', c=1, timelimit_s=24*3600, wait_s=180):
    for cmd in cmds:
        sbatch = SBATCH.format(p, mem, c, cmd)
        pr = sub.Popen(sbatch, stdout=sub.PIPE, stderr=sub.PIPE, stdin=sub.PIPE, shell=True)
        err, out = pr.communicate()
    tstart = time.time()
    missing = list(sentinels)
    while True:
        missing = list(f for f in missing if not os.path.exists(f))
        if not missing or time.time() - tstart > timelimit_s:
            break
        time.sleep(wait_s)


if __name__ == '__main__':

    ap = argparse.ArgumentParser()
    ap.add_argument('--outdir', '-o', default='rd_maps_out')
    ap.add_argument('--force-overwrite', '-f', action='store_true')
    ap.add_argument('--continue-previous', '-p', action='store_true')
    ap.add_argument('--verbose', '-v', action='store_true')
    ap.add_argument('--min-libs', type=int, default=24)
    ap.add_argument('--min-het', type=int, default=5)
    ap.add_argument('--min-hom', type=int, default=3)
    ap.add_argument('sampledirs_file')
    args = ap.parse_args()
   
    print(args)
    assert 0 < args.min_libs
    assert 0 < args.min_het
    assert 0 < args.min_hom
    assert os.path.exists(args.sampledirs_file)
    assert not os.path.exists(args.outdir) or args.force_overwrite or args.continue_previous
    assert (args.force_overwrite and not args.continue_previous) or (not args.force_overwrite and args.continue_previous)
    if os.path.exists(args.outdir) and args.force_overwrite: #  not args.continue_previous:
        import shutil
        shutil.rmtree(args.outdir)
    if not os.path.exists(args.outdir):
        os.mkdir(args.outdir) 

    workdir = args.outdir
    VCFMERGE_DONE = os.path.join(workdir, 'vcfmerge.done')
    VCFMERGE_OUTPUT = os.path.join(workdir, 'merged.bed')
    MULTICOV_DONE = os.path.join(workdir, 'multicov.done')
    MULTICOV_OUTPUT = os.path.join(workdir, 'mincov.bed')
    MULTICOV_INPUT = os.path.join(workdir, 'split_mbedfiles_in')
    VCFFILTER_DONE = os.path.join(workdir, 'vcffilter.done')
    HHFILTER_DONE = os.path.join(workdir, 'hhfilter.done')
    DIFILTER_DONE = os.path.join(workdir, 'difilter.done')
    INTERSECT_DONE = os.path.join(workdir, 'intersect.done')
    ALL_DONE = os.path.join(workdir, 'all.done')
    FINAL_OUTPUT = os.path.join(workdir, 'unique_snps.vcf')

    with open(args.sampledirs_file) as _in:
        directories = list(line.strip() for line in _in)

    vcffiles = list(os.path.join('/'.join(path.strip().split('/')), path.strip().split('/')[-1] + '.iwgsc10.vcf.nod.gz') for path in directories)
    assert all(map(os.path.exists, vcffiles))
    bamfiles = list(os.path.join('/'.join(path.strip().split('/')), path.strip().split('/')[-1] + '.iwgsc10.bam') for path in directories)
    assert all(map(os.path.exists, bamfiles))
    if args.verbose:
        print(bamfiles, vcffiles) 

    if not os.path.exists(VCFMERGE_DONE): 
        cmd = VCFMERGE.format(' '.join(vcffiles), VCFMERGE_OUTPUT, VCFMERGE_DONE)
        slurmJob(cmd, [VCFMERGE_DONE], mem='64GB')

    if not os.path.exists(MULTICOV_DONE):
        mbedfiles = list(splitInputFile(VCFMERGE_OUTPUT))
        with open(MULTICOV_INPUT, 'w') as mbedfiles_f:
            print(*mbedfiles, sep='\n', file=mbedfiles_f)
        donefiles = list(f + '.done' for f in mbedfiles)
        mincovfiles = list(f + '.mincov' for f in mbedfiles)
        cmd = MULTICOV_ARRAY.format(MULTICOV_INPUT, ' '.join(bamfiles), args.min_libs)
        slurmJob(cmd, donefiles, arraylen=len(mbedfiles))
        
        with open(MULTICOV_DONE, 'w'):
            pass
    else:
        mincovfiles = glob.glob(os.path.join(workdir, '*.mincov'))
    
    cmd = 'cat {} > {};'.format(' '.join(mincovfiles), MULTICOV_OUTPUT)
    pr = sub.Popen(cmd, shell=True, stdin=sub.PIPE, stdout=sub.PIPE, stderr=sub.PIPE)
    out, err = pr.communicate()

    if not os.path.exists(DIFILTER_DONE):
        donefiles, fvcffiles, jobs = list(), list(), list()
        for vf in vcffiles:
            prefix = os.path.join(workdir, os.path.basename(vf))
            donefiles.append(prefix + '.done')
            fvcffiles.append(prefix.replace('.vcf.nod.gz', '.filtered.vcf.gz'))
            cmd = DIFILTER.format(vf, MULTICOV_OUTPUT, fvcffiles[-1], donefiles[-1], args.min_het, args.min_hom)
            jobs.append(cmd)
        slurmMultiJob(jobs, donefiles)

        with open(DIFILTER_DONE, 'w'):
            pass 	

    else:
        fvcffiles = list(os.path.join(workdir, os.path.basename(vf)).replace('.vcf.nod.gz', '.filtered.vcf.gz') for vf in vcffiles) 

    if not os.path.exists(INTERSECT_DONE):
        donefiles, uniqfiles, jobs = list(), list(), list()
        isect_sets = set(fvcffiles)
        for vf in fvcffiles:
            donefiles.append(vf + '.done')
            uniqfiles.append(vf.replace('.filtered.vcf.gz', '.uniq.vcf'))
            cmd = INTERSECT.format(vf, ' '.join(isect_sets.difference([vf])), uniqfiles[-1], donefiles[-1])
            jobs.append(cmd)
        slurmMultiJob(jobs, donefiles)

        with open(INTERSECT_DONE, 'w'):
            pass

    else:
        uniqfiles = list(vf.replace('.filtered.vcf.gz', '.uniq.vcf') for vf in fvcffiles)

    snps = list()
    for uf in uniqfiles: 
        sampleID = os.path.basename(uf).split('.')[0]
        with open(uf) as _in:
            for row in csv.reader(_in, delimiter='\t'):
                if not row[0].startswith('#'):
                    row[1], row[2] = int(row[1]), sampleID
                    snps.append(row)
    with open(FINAL_OUTPUT, 'w') as final_out:
        for line in open(uniqfiles[0]):
            if line.startswith('#'):
                if not line.startswith('##reference') and not line.startswith('##source'):
                    print(line.strip(), file=final_out)
            else:
                break
        for snp in sorted(snps, key=lambda x:(x[0], x[1])):
            print(*snp, sep='\t', file=final_out)

    with open(ALL_DONE, 'w'):
        pass

