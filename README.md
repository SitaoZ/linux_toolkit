## Linux_toolkit
Linux commands and tricks in bioinformatics


### One-liner
```bash
$ history | awk '{a[$2]++}END{for(i in a){print a[i]" "i}}' | sort -rn | head # 列出常用的命令
$ wtach vmstat -sSM     # 实时监控
$ vmstat -sSM           # 监控一次
$ du -h -d 1 | sort -rh # 找出最大文件夹
```
### Conda tips
1.创建环境
```bash
$ conda --version          # 显示conda版本
$ conda create -n env_name python=3.7.2 # 创建环境
$ conda activate env_name  # 激活环境
$ conda deactivate         # 退出环境
$ conda env list           # 显示当前所有环境
$ conda info --env         # 显示当前所有环境
```

2. 清理安装包的缓存
```bash
$ conda clean -h
$ conda clean -a # 快速删除
$ conda clean -p # 从可写包缓存中删除没有使用的包，但是不会检查已经安装的包是否有软连接到其中的缓存包
$ conda clean -t # 一键删除anaconda pkgs下面的压缩包
```

### BAM
```bash
$ samtools view QC.sort.bam | grep "XM:i:0" > noMismatch.sam # 找出没有mismatch的比对
$ samtools view QC.sor.bam | grep "AS:" | grep –v "XS:" > unique_alignments.sam # 从bowtie2筛选唯一比对
$ samtools view reads.bam | grep 'XT:A:U' | samtools view -bS -T referenceSequence.fa - > reads.uniqueMap.bam # 从bwa比对文件中筛选唯一比对
$ (samtools view -H QC.sort.bam; samtools view QC.sort.bam | grep -P "\tNH:i:1\t|\tNH:i:1$" | grep -v "ZS:i" ) | samtools view -bS - > unique.bam # 从hisat2中筛选唯一比对
```
