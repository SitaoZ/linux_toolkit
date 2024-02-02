## Table of content
* [One-liner](#One-liner)
*  
## One-liner
### 基本文件处理
```bash
$ history | awk '{a[$2]++} END{for(i in a){print a[i]" "i}}' | sort -rn | head # 列出常用的命令

$ watch vmstat -sSM     # 实时监控
$ vmstat -sSM           # 监控一次

$ du -h -d 1 | sort -rh # 找出最大文件夹

$ cat file.txt | sort | uniq -c | sort -k1nr | head # 排序找出现最多的

$ echo $PATH | tr ":" "\n" | nl # 打印全部路径按行排列

$ sed '/^$/d' file.txt  # 去除空白行
$ grep -v '^$' file.txt # 去除空白行
$ awk '/./' file.txt    # 去除空白行
$ cat file.txt | tr -s "\n" # 去除空白行

$ echo 'ATTGCTATGCTNNNT' | rev | tr 'ACTG' 'TGAC' # 反向互补序列

$ cat a b | sort | uniq -d | wc # 求两个文件的交集
$ cat a b | sort | uniq | wc    # 两个文件的并集
$ cat a b b | sort | uniq -u    # a - b 差集
$ cat a b a | sort | uniq -u    # b - a 差集

$ 打印文件第一行
$ cat file.txt | sed -n '1p' # -n 表示silence模式，只有命令中指定的行才会被打印,1..n都可以
$ cat file.txt | awk -F',' 'NR==1{print;next} {print $1}'
```
### fastq处理
```
$ # fastq长度分布
$ zcat file.fastq.gz | awk 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}'  
$ # 过滤掉小片段fastq
$ awk 'BEGIN {OFS = "\n"} {header = $0 ; getline seq ; getline qheader ; getline qseq ; \
  if (length(seq) >= 10000) {print ">"header, seq}}' < input_reads.fastq > filtered_gt10kb.fasta 
$ # 保留区间的长度
$ awk 'BEGIN {OFS = "\n"} {header = $0 ; getline seq ; getline qheader ; getline qseq ; \
  if (length(seq) >= 10000 && length(seq) <= 20000) {print header, seq, qheader, qseq}}' < input.fastq > filtered_10kb-20kb.fastq
$ # 计算fastq的碱基数
$ awk 'BEGIN{sum=0;}{if(NR%4==2){sum+=length($0);}}END{print sum;}' sequences.fastq 
$ 
$ 交错排布read1和2
$ paste <(paste - - - - < reads-1.fastq) \
      <(paste - - - - < reads-2.fastq) \
    | tr '\t' '\n' \
    > reads-int.fastq
$ 分开read1和2
$ paste - - - - - - - - < reads-int.fastq \
    | tee >(cut -f 1-4 | tr '\t' '\n' > reads-1.fastq) \
    | cut -f 5-8 | tr '\t' '\n' > reads-2.fastq

$ 合并单细胞数据，合并cell barcodes 和 UMI, 名称 和 + 只用R2的，序列和质量值合并
$ paste <(zcat Sample01_S1_R2_001.fastq.gz) \
      <(zcat Sample01_S1_R3_001.fastq.gz) | \
      awk -F '\t' '{ if(NR%4==1||NR%4==3) {print $1} else {print $1 $2} }' | \
      gzip > Sample01_S1_CB_UMI.fastq.gz
```

### fasta 处理
```bash
$ # 将染色体分开
$ cat Homo_sapiens.GRCh38.chr.dna.toplevel.fa | awk '{
        if (substr($0, 1, 1)==">") {filename=(substr($1,2) ".fasta")}
        print $0 >> filename
        close(filename)
}'

```
- seqkit   
seqkit工具可以快速处理fasta文件 [seqkit](https://bioinf.shenwei.me/seqkit/usage/)。
```bash
$ # 按照长度过滤,选取长度小于300bp的fasta子集
$ seqkit seq -M 300 refer.fasta -o lt300.fa
$ # 使用seqkit 选取特定的子集，使用grep子命令 
$ seqkit grep -n -f wanted_gene.csv refer.fasta -o wanted.fa

$ # 全称匹配才行
$ seqkit grep -n -f id.txt swissprot/swissprot.fa -o result_swiss.fa

$ # 正则匹配
$ seqkit grep -r -f id.txt TrEMBL/uniprot_trembl.fasta -o result.fa
```

```bash
$ # faFilter obtained from http://hgdownload.soe.ucsc.edu/admin/exe/
$ faFilter -minSize=N -maxSize=N in.fa out.fa
```




### samtools处理
SAMtools是li heng开发的用于比对文件处理的利器[samtools](http://www.htslib.org/)。
```bash
$ # 不同比对软件的tag有细微差异，注意区分
$ samtools view QC.sort.bam | grep "XM:i:0" > noMismatch.sam # 找出没有mismatch的比对
$ samtools view QC.sor.bam | grep "AS:" | grep –v "XS:" > unique_alignments.sam # 从bowtie2筛选唯一比对
$ samtools view reads.bam | grep 'XT:A:U' | samtools view -bS -T referenceSequence.fa - > reads.uniqueMap.bam # 从bwa比对文件中筛选唯一比对
$ (samtools view -H QC.sort.bam; samtools view QC.sort.bam | grep -P "\tNH:i:1\t|\tNH:i:1$" | grep -v "ZS:i" ) | samtools view -bS - > unique.bam # 从hisat2中筛选唯一比对
$
$ # 计算基因组覆盖度
$ samtools depth my.bam | awk '{sum+=$3} END { print "Average = ",sum/NR}'
$ # 计算基因组大小
$ samtools view -H my.bam | grep -P '^@SQ' | cut -f 3 -d ':' | awk '{sum+=$1} END {print sum}'
$ # 添加MD标签, MD标签用于SNP calling
$ samtools calmd -b <my.bam> <ref_genome.fasta> > my_md.bam
```
[MD tag](https://vincebuffalo.com/notes/2014/01/17/md-tags-in-bam-files.html)

### 链接处理

```bash
$ readlink -f/--canonicalize # 检查链接地址，且地址不能失效，失效报错
$ readlink -m/--canonicalize-missing # 返回地址，且地址可以失效
$ # 循环更改失效路径
$ for i in input ip;
  do
   for j in 1 2 3;
    do
		   ln -s `readlink -m ${i}${j}.bam | sed 's/02.hisat2/03.hisat2/';rm ${i}${j}.bam` ${i}${j}.bam 
	   done
  done
```

### gtf & gff
```bash
$ # gff -> gtf
$ gffread in.gff -T -o out.gtf

$ # gtf -> gff
$ gffread in.gtf -o out.gff

$ # 提取转录本序列
$ gffread -w transcripts.fa -g mm10.fa mm10.gtf

$ # 提取CDS序列
$ gffread -x cds.fa -g mm10.fa mm10.gtf
```
