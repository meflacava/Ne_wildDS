#!/bin/bash
#SBATCH -e call_genotypes.j%j.err
#SBATCH -o call_genotypes.j%j.out
#SBATCH -J angsd
#SBATCH --time=2-00:00:00
#SBATCH --ntasks=24
#SBATCH --mem=128G
#SBATCH -p bigmemh
#SBATCH --mail-type=END,FAIL

# call script with:
# sbatch call_genotypes_withDepthfilter.sh <bamlist> <output name>


module load angsd

nInd=$(wc -l $1 | awk '{print $1}')
mInd=$(echo "$nInd * 0.75" | bc | awk '{print int($1+0.5)}')

angsd \
  -bam "$1" \
  -out "$2" \
  -rf Hyp_tra_F_20210429.loci \
  -nThreads 24 \
  -nQueueSize 2 \
  -minMapQ 20 \
  -minQ 20 \
  -minInd ${mInd} \
  -GL 1 \
  -doGLF 2 \
  -doMaf 2 \
  -doPost 1 \
  -postCutoff 0.85 \
  -minMaf 0.02 \
  -SNP_pval 1e-6 \
  -doIBS 1 \
  -doCounts 1 \
  -dumpCounts 2 \
  -doMajorMinor 1 \
  -makeMatrix 1 \
  -doCov 1 \
  -doHWE 1 \
  -minHWEpval 0.01 \
  -doGeno 4 \
  -geno_minDepth 5 \
  -geno_maxDepth 50
