(samtools view -H QC.sort.bam; samtools view QC.sort.bam | grep -P "\tNH:i:1\t|\tNH:i:1$" | grep -v "ZS:i" ) | samtools view -bS - > unique.bam
# samtools view -H QC.sort.bam 
# 保留bam的头文件

# NH:i:<N> : The number of mapped locations for the read or the pair.
# http://daehwankimlab.github.io/hisat2/manual/
# NH:i:1 表示 read比对上的位置的数目只有一个,即唯一比对

# ZS:i:<N> : Alignment score for the best-scoring alignment found other than the alignment reported.
# ZS:i : 当有多个比对位置时候才出现

