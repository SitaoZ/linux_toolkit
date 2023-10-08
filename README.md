## Linux_toolkit
Linux commands and tricks in bioinformatics


### One-liner
```bash
$ history | awk '{a[$2]++}END{for(i in a){print a[i]" "i}}' | sort -rn | head # 列出常用的命令
$ wtach vmstat -sSM # 实时监控
$ vmstat -sSM # 监控一次
$ du -h -d 1 | sort -rh # 找出最大文件夹
```
### Conda tips
1. 清理安装包的缓存
```bash
$ conda clean -h
$ conda clean -p #
$ conda clean -t #
```

### BAM
```bash
$ samtools view QC.sort.bam | grep "XM:i:0" > noMismatch.sam # 找出没有mismatch的比对
$ samtools view QC.sor.bam | grep "AS:" | grep –v "XS:" > unique_alignments.sam # 从bowtie2筛选唯一比对
$ samtools view reads.bam | grep 'XT:A:U' | samtools view -bS -T referenceSequence.fa - > reads.uniqueMap.bam # 从bwa比对文件中筛选唯一比对
$ (samtools view -H QC.sort.bam; samtools view QC.sort.bam | grep -P "\tNH:i:1\t|\tNH:i:1$" | grep -v "ZS:i" ) | samtools view -bS - > unique.bam # 从hisat2中筛选唯一比对
```
