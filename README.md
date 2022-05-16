# Scripts I have written or optimized for bulk and single-cell sequencing analysis
## 1. 220408_STARsoloOutputSPLiTseq_MergePrimers.py
Script to merge poly-T and random hexamer barcodes that were pooled together in one well during the first barcoding step (RT reaction) of the SPLiTseq protocol as they probably originate from the same cell. The script also assigns each read to its matching sample sam file.\
Test SAM file was created by running STARsolo with the following command: \
STAR --runMode alignReads --runThreadN 4 --genomeDir /home/jenny/Documents/Annotation/m39_STAR --readFilesIn S91_SL1_R1.fastq.gz S91_SL1_R2.fastq.gz --readFilesCommand zcat --outFileNamePrefix ./Bl6prion_S91_SL1. --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM --outSAMtype BAM SortedByCoordinate --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 0 --soloType CB_UMI_Complex --soloCBwhitelist /home/jenny/Documents/barcodeList.txt /home/jenny/Documents/barcodeList.txt /home/jenny/Documents/barcodeList.txt --soloBarcodeReadLength 0 --soloCBposition 0_10_0_17 0_48_0_55 0_86_0_93 --soloUMIposition 0_0_0_9 --soloCBmatchWLtype EditDist_2 --soloFeatures Gene GeneFull Velocyto --soloMultiMappers EM --soloCellReadStats Standard\
I also included an example 1st barcode-sample association file for review.\
