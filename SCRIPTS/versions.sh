

python --version 2>&1 | sed -e "s/^Python/Python:/"
bwa 2>&1 | grep "^Version" | sed -e "s/Version/bwa/"
samtools --version 2>&1 | head -1 | sed -e "s/samtools/samtools:/"
fastqc --version 2>&1 | sed -e "s/FastQC/FastQC:/"
multiqc --version 2>&1 | sed -e "s/multiqc, version/multiQC:/"

source activate dgeDE
Rscript --version 2>&1 | sed -e "s/.*version/R:/"
Rscript -e "packageVersion('DESeq2')" | sed -e "s/\[1\]/DESeq2:/"
Rscript -e "packageVersion('limma')" | sed -e "s/\[1\]/limma:/"
