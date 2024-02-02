## One-liner
* [基本文件处理](#基本文件处理)
* [fastq处理](#fastq处理)
* [fasta处理](#fasta处理)
* [bam处理](#bam处理)
* [链接处理](#链接处理)
* [gtf_gff](#gtf_gff)
* [vcf](#vcf)

### 基本文件处理
- 找出常用命令
```bash
$ history | awk '{a[$2]++} END{for(i in a){print a[i]" "i}}' | sort -rn | head # 列出常用的命令
```

- 内存监控
```
$ watch vmstat -sSM     # 实时监控
$ vmstat -sSM           # 监控一次
```
- 找出最大文件夹
```bash 
$ du -h -d 1 | sort -rh # 找出最大文件夹
```
- 找出文件中出现次数最多的
```bash
$ cat file.txt | sort | uniq -c | sort -k1nr | head # 排序找出现最多的
```
- 打印全部路径
```bash
$ echo $PATH | tr ":" "\n" | nl # 打印全部路径按行排列
```


- 打印文件第一行
```bash
$ cat file.txt | sed -n '1p' # -n 表示silence模式，只有命令中指定的行才会被打印,1..n都可以
$ cat file.txt | awk -F',' 'NR==1{print;next} {print $1}'
```

- 去除空白行
```bash
$ sed '/^$/d' file.txt  # 去除空白行
$ grep -v '^$' file.txt # 去除空白行
$ awk '/./' file.txt    # 去除空白行
$ cat file.txt | tr -s "\n" # 去除空白行
```

- 打印所有3mer的DNA组合
```bash
echo {A,C,T,G}{A,C,T,G}{A,C,T,G}
```
- 序列反向互补
```bash
echo 'AGTCATGCAGTGCNNNNT' | rev | tr 'ACTG' 'TGAC'

echo 'AGTCATGCAGTGCNNNNT' | python -c "import sys;from Bio.Seq import Seq;a = [print(Seq(i.strip()).reverse_complement()) for i in sys.stdin];"
```

- 文件的交集、并集、差集
```bash
#  使用sort和uniq快速求文件的交集、并集、差集
# 主要使用的uniq 
# uniq -d (only print duplicate lines, one for each group)
# uniq -u (only print unique lines)
# 并集
sort a.txt b.txt | uniq | wc

# 交集
sort a.txt b.txt | uniq -d | wc

# 差集 1
# a.txt - b.txt 
sort a.txt b.txt b.txt | uniq -u | wc

# b.txt - a.txt 
sort b.txt a.txt a.txt | uniq -u | wc
```
- 删除文件
```bash
# wildcard
rm *bam

# find and xargs
find . -name "*.bam" | xargs rm
```
- 文件检查
```bash
# 生成一个文件的MD5值
md5sum xxx.csv > md5.txt

# 检查文件是否完整
md5sum -c md5.txt 

# or 
md5sum -c --status md5.txt
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

### fasta处理
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

- split_multi_fasta.sh
```bash
$ awk '/^>/{s=++d".fa"}{print > s}' xxx.fa
# /^>/ 匹配fasta 文件名
# s=++d".fa" 生成文件名，一每一个fasta格式自追加编号
# {print > s} 输出到该文件名

$ csplit -z -q -n 4 -f sequence_ sequences.fasta /\>/ {*}  
# /[正则表达式]/   #匹配文本样式，比如/SERVER/，从第一行到包含SERVER的匹配行。
# {*}     #表示根据匹配重复执行分割，直到文件尾停止，使用{整数}的形式指定分割执行的次数。
# -z或--elide-empty-files 删除长度为0 Byte文件。
# -q或-s或--quiet或--silent 不显示指令执行过程。
# -n      #指定分割后的文件名后缀的数字个数。比如01、02、03等。
# -f      #指定分割后的文件名前缀。
```
- 去除空白
```bash
$ sed /^$/d
```


### bam处理
SAMtools是li heng开发的用于比对文件处理的利器[samtools](http://www.htslib.org/)。

- bowtie2提取唯一比对文件
```bash
$ samtools view QC.sort.bam | grep "AS:" | grep –v "XS:" > unique_alignments.sam
# http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml
# AS:i:<N>  Alignment score. Can be negative. Can be greater than 0 in --local mode (but not in --end-to-end mode). Only present if SAM record is for an aligned read.
# AS:i:<N> 只有比对上的read才有该标签

# XS:i:<N> Alignment score for the best-scoring alignment found other than the alignment reported. Can be negative. Can be greater than 0 in --local mode (but not in --end-to-end mode). Only present if the SAM record is for an aligned read and more than one alignment was found for the read. Note that, when the read is part of a concordantly-aligned pair, this score could be greater than AS:i.
# XS:i:<N> 当read有多个比对位置时候，才出现该标签
```
- hisat2提取唯一比对文件
```bash
$ (samtools view -H QC.sort.bam; samtools view QC.sort.bam | grep -P "\tNH:i:1\t|\tNH:i:1$" | grep -v "ZS:i" ) | samtools view -bS - > unique.bam
# samtools view -H QC.sort.bam 
# 保留bam的头文件

# NH:i:<N> : The number of mapped locations for the read or the pair.
# http://daehwankimlab.github.io/hisat2/manual/
# NH:i:1 表示 read比对上的位置的数目只有一个,即唯一比对

# ZS:i:<N> : Alignment score for the best-scoring alignment found other than the alignment reported.
# ZS:i : 当有多个比对位置时候才出现
```
- bwa提取唯一比对文件
```bash
$ samtools view reads.bam | grep 'XT:A:U' | samtools view -bS -T referenceSequence.fa - > reads.uniqueMap.bam
```

- 找出没有mismatch的比对
```bash
$ samtools view QC.sort.bam | grep "XM:i:0" > noMismatch.sam
```
- 计算基因组覆盖度
```bash 
$ samtools depth my.bam | awk '{sum+=$3} END { print "Average = ",sum/NR}'
```
- 计算基因组大小
```bash 
$ samtools view -H my.bam | grep -P '^@SQ' | cut -f 3 -d ':' | awk '{sum+=$1} END {print sum}'
```
-  添加MD标签, MD标签用于SNP calling
```bash
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

### gtf_gff
- gffread处理
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
- 提取TSS信息
```bash
# tss
cat xxx.gff3 | awk '$3=="gene" {print $0}' | grep protein_coding | awk -v OFS="\t" '{if ($7=="+") {print $1, $4, $4+1}} else {print $1, $5-1, $5}' > tss.bed

# promoter 5k upstream from tss
cat xxx.gff3 | awk '$3=="gene" {print $0}' | grep protein_coding | awk -v OFS="\t" '{if ($7=="+") {print $1, $4, $4+5000} else {print $1, $5-5000, $5}}' > promoters.bed
```

### vcf 
- 按染色体提取vcf
```bash
for chr in {1..22}; do awk -v chr=$chr -v threshold=10 '$1 == chr && $6 > threshold' xxx.vcf > xxx_chr$chr.vcf; done
# https://samtools.github.io/hts-specs/VCFv4.2.pdf
# shell中的变量在awk中无法识别，可以使用-v指定变量传入awk
# -v chr=$chr, 将shell中的变量$chr传入awk


awk '$1 == 5 && $7 == "PASS"' xxx.vcf > xxx_filter.vcf
# 挑选五号染色体上的SNP
# 过滤值为PASS
```
