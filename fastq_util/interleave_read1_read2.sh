# 将read1和read2交互排列成一个文件
paste <(paste - - - - < reads1.fastq) \
      <(paste - - - - < reads2.fastq) \
    | tr '\t' '\n' \
    > reads_interleave.fastq


# 分离成两个文件
paste - - - - - - - - < reads-interleave.fastq \
    | tee >(cut -f 1-4 | tr '\t' '\n' > reads1.fastq) \
    | cut -f 5-8 | tr '\t' '\n' > reads2.fastq
    
