#!/bin/bash
#SBATCH -J setup_angsd
#SBATCH -o setup_angsd.j%j.out
#SBATCH --time=00:05:00
#SBATCH -p low

# Usage:
# sbatch multiANGSD_minMAF0.05_Dp5-100.sh list_of_bamlists

list=$1
wc=$(wc -l < "${list}")
x=1

while [[ $x -le $wc ]]; do
    bamlist=$(sed -n "${x}p" "${list}")
    yr="${bamlist%%.*}"
    outname="${yr}_minMAF0.02_Dp5-50"
    scriptname="${yr}_minMAF0.02_Dp5-50_call_genotypes.sh"

    echo "#!/bin/bash" > "$scriptname"
    echo "#SBATCH -e ${outname}.j%j.err" >> "$scriptname"
    echo "#SBATCH -o ${outname}.j%j.out" >> "$scriptname"
    echo "#SBATCH -J ${outname}" >> "$scriptname"
    echo "#SBATCH --time=01:00:00" >> "$scriptname"
    echo "#SBATCH --ntasks=12" >> "$scriptname"
    echo "#SBATCH -p med" >> "$scriptname"

    echo "module load angsd" >> "$scriptname"

    nInd=$(wc -l < "$bamlist")
    mInd=$((nInd / 2))

    echo "nInd=${nInd}" >> "$scriptname"
    echo "mInd=${mInd}" >> "$scriptname"

    echo "angsd \\" >> "$scriptname"
    echo "  -bam ${bamlist} \\" >> "$scriptname"
    echo "  -out ${outname} \\" >> "$scriptname"
    echo "  -rf Hyp_tra_F_20210429.loci \\" >> "$scriptname"
    echo "  -nThreads 12 \\" >> "$scriptname"
    echo "  -nQueueSize 2 \\" >> "$scriptname"
    echo "  -minMapQ 20 \\" >> "$scriptname"
    echo "  -minQ 20 \\" >> "$scriptname"
    echo "  -minInd ${mInd} \\" >> "$scriptname"
    echo "  -GL 1 \\" >> "$scriptname"
    echo "  -doGLF 2 \\" >> "$scriptname"
    echo "  -doMaf 2 \\" >> "$scriptname"
    echo "  -doPost 1 \\" >> "$scriptname"
    echo "  -postCutoff 0.85 \\" >> "$scriptname"
    echo "  -minMaf 0.02 \\" >> "$scriptname"
    echo "  -SNP_pval 1e-6 \\" >> "$scriptname"
    echo "  -doIBS 1 \\" >> "$scriptname"
    echo "  -doCounts 1 \\" >> "$scriptname"
    echo "  -dumpCounts 2 \\" >> "$scriptname"
    echo "  -doMajorMinor 1 \\" >> "$scriptname"
    echo "  -makeMatrix 1 \\" >> "$scriptname"
    echo "  -doCov 1 \\" >> "$scriptname"
    echo "  -doHWE 1 \\" >> "$scriptname"
    echo "  -minHWEpval 0.01 \\" >> "$scriptname"
    echo "  -doGeno 4 \\" >> "$scriptname"
    echo "  -geno_minDepth 5 \\" >> "$scriptname"
    echo "  -geno_maxDepth 50" >> "$scriptname"

    sbatch "$scriptname"
    x=$((x + 1))
done
