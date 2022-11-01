samtools view QC.sor.bam | grep "AS:" | grep –v "XS:" > unique_alignments.sam
# http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml
# AS:i:<N>  Alignment score. Can be negative. Can be greater than 0 in --local mode (but not in --end-to-end mode). Only present if SAM record is for an aligned read.
# AS:i:<N> 只有比对上的read才有该标签

# XS:i:<N> Alignment score for the best-scoring alignment found other than the alignment reported. Can be negative. Can be greater than 0 in --local mode (but not in --end-to-end mode). Only present if the SAM record is for an aligned read and more than one alignment was found for the read. Note that, when the read is part of a concordantly-aligned pair, this score could be greater than AS:i.
# XS:i:<N> 当read有多个比对位置时候，才出现该标签
