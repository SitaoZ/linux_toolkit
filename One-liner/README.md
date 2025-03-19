## One-liner
* [基本文件处理](#基本文件处理)
* [fastq处理](#fastq处理)
* [fasta处理](#fasta处理)
* [bam处理](#bam处理)
* [链接处理](#链接处理)
* [gtf_gff](#gtf_gff)
* [vcf](#vcf)
* [table](#table)
* [http](#http)
* [命令执行](#命令执行)
* [服务器配置信息](#服务器配置信息)
### 基本文件处理
- 找出常用命令
```bash
$ history | awk '{a[$2]++} END{for(i in a){print a[i]" "i}}' | sort -rn | head # 列出常用的命令
```

- 批量给文件末尾添加一个换行符
```
$ for i in `ls *txt`; do echo >> $i  ;done

$ sed -i '$a\' filename.sh

$ printf "\n" >> filename.sh

```

- 循环脚本
```bash
$ for f in $(ls *.sh); do wc ${f}; done
```
- 一次性杀掉该用户的全部任务
```bash
$ ps -ef | grep zhusitao | awk '{ print $2 }'| xargs kill -9
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

- tab分隔文件排序
```bash
$ sort -t $'\t' -rnk2 xxx.txt
```

- 前一个命令的输出作为后一个命令的参数
```bash
$ # 1. `` 反斜杠
$ echo `ls *bam`

$ # 2. $() 效果等同于反斜杠
$ echo $(ls *bam)

$ # 3. xargs ,缺点在于参数过多(成千上万),导致后面的命令执行失败
$ find . -name "*gz" | xargs rm
$ # 批量查看文件是否存在
$ cat KO_WT_SampleList.txt | cut -d " "  -f 2 | xargs ls -l

$ # 4. find -exec执行，不受参数数量的影响，克服xargs的缺点
$ find . -name "*gz" -exec ls {} \;

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
$ echo {A,C,T,G}{A,C,T,G}{A,C,T,G}
```
- 序列反向互补
```bash
$ echo 'AGTCATGCAGTGCNNNNT' | rev | tr 'ACTG' 'TGAC'

$ echo 'AGTCATGCAGTGCNNNNT' | python -c "import sys;from Bio.Seq import Seq;a = [print(Seq(i.strip()).reverse_complement()) for i in sys.stdin];"
```

- 文件的交集、并集、差集
```bash
#  使用sort和uniq快速求文件的交集、并集、差集
# 主要使用的uniq 
# uniq -d (only print duplicate lines, one for each group)
# uniq -u (only print unique lines)
# 并集
$ sort a.txt b.txt | uniq | wc

# 交集
$ sort a.txt b.txt | uniq -d | wc

# 差集 1
# a.txt - b.txt 
$ sort a.txt b.txt b.txt | uniq -u | wc

# b.txt - a.txt 
$ sort b.txt a.txt a.txt | uniq -u | wc
```
- 删除文件
```bash
# wildcard
$ rm *bam

# find and xargs
$ find . -name "*.bam" | xargs rm
$ find . -name "*.sam" | xargs rm -rf 
```
- 文件检查
```bash
# 生成一个文件的MD5值
$ md5sum xxx.csv > md5.txt

# 检查文件是否完整
$ md5sum -c md5.txt 

# or 
$ md5sum -c --status md5.txt
```
- 子表格生成
```bash
$ cat <(less total.csv | head -n 1) <(less total.csv| grep AT3G49430) > result.csv
# < 表示重定向输入符号,usage: command < file 
# () 表示整体执行括弧中的命令
# head -n 1 表示表头
# <() 接 grep 表示另外一个筛选的内容
# cat 合并; > 重定向文件
```
- 找到含有特定字符的文件名
```bash
$ grep -lir "All_isoforms_peak_process" * # 在当前目录下找出含有All_isoforms_peak_process字符的文件名
$ # -l --files-with-matches 打印匹配的文件名
$ # -i --ignore-case 忽略大小写
$ # -r --recursive 循环当前目录所有文件
```

- 输出上一次执行的命令，并执行
```bash
$ !!
```

- 输出上一次命令执行的参数，并执行
```bash
$ !*
```

- 删除或者列出正则匹配的文件
```bash
$ ls !(*.gz)        # 列出没有匹配的gz的文件
$ rm !(*.foo|*.bar) # 删除没有匹配的*foo或者*bar
```

- fc 固定命令 fix command，和history类似，多了编辑功能
```bash
$ # -e ename 指定的编辑器
$ # -l list 列出的命令行数
$ # -n 去除行号
$ # -r 倒序列出命令，需要配合-l使用

$ fc -l # 列出最后执行的16个命令
$ fc -l -3 # 列出最后执行的3个命令
$ fc -ln # 列出命令，但不现实行号
$ 
```
- 原样结构目录备份脚本和程序

将下面的脚本写入 backup.sh; 执行 `bash backup.sh dirname`; 执行当前目录备份

```bash
for f in `find ../$1/ -name "*.py" | sed 's#\.\.\/#\.\/#g'`
do 
    DIRNAME=`dirname $f`
    BASENAME=`basename $f`
    mkdir -p $DIRNAME
    cp ../$DIRNAME/$BASENAME $DIRNAME/$BASENAME
done

# shell scripts 

for f in `find ../$1/ -name "*.sh" | sed 's#\.\.\/#\.\/#g'`
do
    DIRNAME=`dirname $f`
    BASENAME=`basename $f`
    mkdir -p $DIRNAME
    cp ../$DIRNAME/$BASENAME $DIRNAME/$BASENAME
done
```

### fastq处理
- read长度分布
```
$ # fastq长度分布
$ zcat xx.fq.gz | awk 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}' > reads.dist

# read species
zcat xx.fq.gz | awk 'NR%4 == 2 {species[$0]++} END {for (s in species) {print s, species[s]}}' > reads.count
```
- 过滤掉小片段fastq
```
$ # 过滤掉小片段fastq
$ awk 'BEGIN {OFS = "\n"} {header = $0 ; getline seq ; getline qheader ; getline qseq ; \
  if (length(seq) >= 10000) {print ">"header, seq}}' < input_reads.fastq > filtered_gt10kb.fasta
```
- 保留区间的长度
```
$ awk 'BEGIN {OFS = "\n"} {header = $0 ; getline seq ; getline qheader ; getline qseq ; \
  if (length(seq) >= 10000 && length(seq) <= 20000) {print header, seq, qheader, qseq}}' < input.fastq > filtered_10kb-20kb.fastq
```

- 计算fastq的碱基数
```bash
$ awk 'BEGIN{sum=0;}{if(NR%4==2){sum+=length($0);}}END{print sum;}' sequences.fastq 
```
- fastq 变成一行
```bash
$ zcat sample.fastq.gz | paste - - - - | gzip > sample.one_line.fastq.gz
```
- fastq to fasta
```bash
$ sed -n '1~4s/^@/>/p;2~4p' file.fq > file.fa
```

- 交错排布read1和2
```bash
$ paste <(paste - - - - < reads-1.fastq) \
      <(paste - - - - < reads-2.fastq) \
    | tr '\t' '\n' \
    > reads-int.fastq
```
- 分开read1和2
```bash
$ paste - - - - - - - - < reads-int.fastq \
    | tee >(cut -f 1-4 | tr '\t' '\n' > reads-1.fastq) \
    | cut -f 5-8 | tr '\t' '\n' > reads-2.fastq
```

- seqkit 按照`read id`去除fastq中的read
 ```bash
$ seqkit grep -v -f exclude_list.txt -n -i ../01.fastp/L1_A_R1.clean.fq.gz -o L1_A_R1.clean.remove.fq.gz
# -v 取反
# -f 输入匹配的字符文本
# -n --by-name
# -i --ignore-case 忽略大小写
# -o --out-file 
```

- 合并单细胞数据，合并cell barcodes 和 UMI, 名称 和 + 只用R2的，序列和质量值合并
```
$ paste <(zcat Sample01_S1_R2_001.fastq.gz) \
      <(zcat Sample01_S1_R3_001.fastq.gz) | \
      awk -F '\t' '{ if(NR%4==1||NR%4==3) {print $1} else {print $1 $2} }' | \
      gzip > Sample01_S1_CB_UMI.fastq.gz
```

- long read
```bash
# long read QC
$ zcat xaa.fq.gz | seqtk seq -A -L 10000 - | grep -v "^>" | tr -dc "ACGTNacgtn" | wc -m

# zcat ( concatenates the compressed fastq files into one stream )
# seqtk ( converts to fasta format and drops reads less than 10k )
# grep ( -v excludes lines starting with “>”, i.e. fasta headers )
# tr ( -dc removes any characters not in set “ACGTNacgtn” )
# wc ( -m counts characters )
```

- read按barcode拆分
```bash
$ paste <(zcat sample_1.fq.gz|paste - - - -) <(zcat sample_2.fq.gz|paste - - - - ) | awk -v FS="\t" -v OFS="\n" 'FNR==NR {samples[$2]=$1; next} {barcode = substr($6,0,6); if(samples[barcode]) { print $1,$2,$3,$4>>samples[barcode]"_1.fq"; print $5,$6,$7,$8>>samples[barcode]"_2.fq"}}' samples.txt -

# 1. paste - - - -, 四行打包成一行, 默认tab分隔符
# 2. awk -v FS="\t" 指定输入符号为\t, OFS="\n" 指定输出分割符为"\n"
# 3. FNR==NR 处理第一个文件; awk 记录数为文件的行号(line number)
# 当前文件行号和总行号相等，表示第一个文件.
# 当前记录数: FNR(Number of input Record in current input File)
# 总记录数: NR(total Number of Records seen so far )
# 4. substr, awk 内置函数,截取字符串,substr(string,position,length)
# 5. samples.txt - 标准输入samples.txt 内容到awk,(即处理的第一个文件)
```
- 提取barcode
```bash
$ less xxx.fastq.gz |  cut -c 1-8 > barcode.csv
# -c --characters=LIST 截取字符
```

- 按模式提取read
```bash
$ zcat reads.fq.gz \
| paste - - - - \
| awk -v FS="\t" -v OFS="\n" '$2 ~ "AAGTTGATAACGGACTAGCCTTATTTT" {print $1, $2, $3, $4}' \
| gzip > filtered.fq.gz
```

- 16s测序数据read1和read2合并
```bash
$ # 使用的是qiime的vsearch工具合并
$ vsearch --fastq_mergepairs seq/xxx_1.fq.gz \
       --reverse seq/xxx_2.fq.gz \
       --fastqout temp/xxx.merged.fq \
       --relabel xxx.
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
$ seqkit locate -f primers.fasta reference.fasta # 找回primer在基因组中的位置
```

```bash
$ seqkit stat rdp_16s_v18.fa # 统计fasta序列的相关信息
$ # file            format  type  num_seqs     sum_len  min_len  avg_len  max_len
$ # rdp_16s_v18.fa  FASTA   DNA     21,195  30,743,106      455  1,450.5    1,968
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


- 查找基因组上的motif序列
```bash
$ zcat hairpin.fa.gz \
    | seqkit locate -i -d -p AUGGACUN \
    | head -n 4 \
    | csvtk pretty -t

$ cat z | seqkit locate -p 'GTGA(C|T)(G|A)T' -r | csvtk preety -t 
```

- 查找ORF
```bash
$ A[TU]G # 起始密码子
$ [TU](?:AG|AA|GA) 终止密码子
$ .{3})+? # 三联体
$ zcat hairpin.fa.gz \
    | seqkit locate -i -p "A[TU]G(?:.{3})+?[TU](?:AG|AA|GA)" -r \
    | head -n 4 \
    | csvtk pretty -t
```

- bed to fasta
```bash
$ # Extend the peak locations by 50 basepairs left and right.
$ N=50
$ bedtools slop -g $REF.fai -i $BED -b $N > extended.bed

$ # Extract the sequences corresponding to these locations.
$ REF=refs/saccer3.fa
$ bedtools getfasta -bed extended.bed -fi $REF -fo sequences.fa

$ # Run meme-chip on the sequences.
$ rm -rf results
$ meme-chip -meme-maxw 10 -meme-minw 5 -o results sequences.fa
```

### bam处理
SAMtools是li heng开发的用于比对文件处理的利器[samtools](http://www.htslib.org/)。

- 整体查看bam
```bash 
$ samtools flagstat xxx.bam # 生成bam的flags报告
$ # samtools flagstat 给出报告中的正确比对数目properly align为正确比对(properly mapped)且主要比对(primary alignment reads)
$ bamtools stats -in SRR1972739.bwa.bam # 生成bam的flags报告
$ # samtools flagstat和bamtools stats给出的正确比对有差异
$ samtools idxstats xxx.bam # 生成每个染色体上比对了多少reads的报告

```

- properly_alignment
```bash
$ samtools flags UNMAP,SECONDARY,SUPPLEMENTARY
$ 0x904	2308	UNMAP,SECONDARY,SUPPLEMENTARY
$ samtools view -c -f 2 -F 2308 xxx.bam # 正确比对AND主要比对
$ samtools view -c -f 2 xxx.bam         # 正确比对

```
到底什么是正确比对proper-pair?

对于PROPER_PAIR(正确比对)SAM文件给出的解释是：每个片段根据比对软件能正确的比对。这个定义不清晰。一般来说正确的proper-pairs是两条read都比对上reference，方向一正一反，且两端边界距离在正确的范围之内。

- primary alignment

对于每条read只会有一个primary alignment和其他的secondary、supplementary alignment。
samtools没有primary的flag，我们主要通过取反secondary和supplementary来实现primary比对的统计
```bash
$ samtools flags SUPPLEMENTARY,SECONDARY
$ # 0x900 2304 SECONDARY,SUPPLEMENTARY
$ samtools view -c -F 4 -F 2304 xxx.bam # 统计primary alignment
```

- secondary alignment
```bash
$ samtools flags SECONDARY
$ # 0x100 256 SECONDARY
$ samtools view -c -F 4 -f 256 xxx.bam

```
- supplementary alignment

也称嵌合比对chimeric alignment,是一种部分匹配基因组不同区域而不重叠的比对。
```bash
$ samtools flags SUPPLEMENTARY
$ # 0x800 2048 SUPPLEMENTARY
$ samtools view -c -F 4 -f 2048 xxx.bam
```


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

- 统计bam中没有比对上的数目
```bash
$ samtools view -c -f 4 xxx.bam # 没有比对上的
$ samtools view -c -F 4 xxx.bam # 比对上的
$ samtools view -b -F 4 xxx.bam > aligned.bam
```

- 选出比对上负链的
```bash
$ samtools flags 4
0x4	4	UNMAP

$ samtools flags 16
0x10	16	REVERSE

$ samtools view -F 4 -f 16 xxx.bam  # 比对上负链的
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
$ cat xxx.gff3 | awk '$3=="gene" {print $0}' | grep protein_coding | awk -v OFS="\t" '{if ($7=="+") {print $1, $4, $4+1}} else {print $1, $5-1, $5}' > tss.bed

# promoter 5k upstream from tss
$ cat xxx.gff3 | awk '$3=="gene" {print $0}' | grep protein_coding | awk -v OFS="\t" '{if ($7=="+") {print $1, $4, $4+5000} else {print $1, $5-5000, $5}}' > promoters.bed
```
- gff提取intron awk脚本
```bash
BEGIN{OFS="\t"}
{ start[NR]=$4; 
  end[NR]=$5; 
  strand[NR]=$7; 
  ID[NR]=$11; 
  chr[NR]=$1; 
  ens[NR]=$9; 
  symb[NR]=$10} 

END { 
    for (i=1; i<=NR; i++){
        if(ID[i]==ID[i+1]){ 
            if(strand[i]=="+"){
                intron_start=start[i+1]-1;
                intron_end=end[i]+1;
                print chr[i], 
                      intron_end, 
                      intron_start,
                      "intron",
                      strand[i], 
                      ens[i],
                      symb[i], 
                      ID[i]

            }   
            else{
                intron_start=end[i]+1;
                intron_end=start[i+1]-1;
                print chr[i], 
                      intron_start, 
                      intron_end,
                      "intron",
                      strand[i], 
                      ens[i], 
                      symb[i], ID[i]
            }   
        }   
    }   
}
```

### vcf 
- 按染色体提取vcf
```bash
$ for chr in {1..22}; do awk -v chr=$chr -v threshold=10 '$1 == chr && $6 > threshold' xxx.vcf > xxx_chr$chr.vcf; done
# https://samtools.github.io/hts-specs/VCFv4.2.pdf
# shell中的变量在awk中无法识别，可以使用-v指定变量传入awk
# -v chr=$chr, 将shell中的变量$chr传入awk


$ awk '$1 == 5 && $7 == "PASS"' xxx.vcf > xxx_filter.vcf
# 挑选五号染色体上的SNP
# 过滤值为PASS

$ # 值提取PASS的位点
$ awk -F '\t' '{if($0 ~ /\#/) print; else if($7 == "PASS") print}' xxx,vcf > xxx_filter.vcf
```

- 替换vcftools 过滤后genotype(GT)不正常的位点 原先是.|.，转换后变成了.
```bash
cat raw.mm8mac1.recode.vcf | sed $'s/\t.:/\t.|.:/g' > clean.mm8mac1.recode.vcf
```

- bcftools
```bash
$ # 使用bcftools的query功能，指定位置和显示的信息
$ bcftools query -r '19:400300-400800' -f '%CHROM\t%POS\t%REF\t%ALT\n' subset_hg19.vcf.gz
$ # This will produce:
$ # 19 400410 CA C
$ # 19 400666 G C
$ # 19 400742 C T
```
- 查看VCF头文件
```bash
$ bcftools view -h  newfile.vcf.gz
```

- 排除缺失位点和ref纯合位点
```bash
$ bcftools view -e 'GT= "." | GT="0|0"' subset_hg19.vcf.gz |bcftools query -f '%POS[\t%GT\t]\n' | head -n 3
$ # 402556 0|1 0|1 1|1 1|0 0|1 1|1
$ # 402707 0|1 0|1 1|1 1|0 0|1 1|1
$ # 402723 0|1 0|1 1|1 1|0 0|1 1|1
```

- 提取indel
```bash
$ bcftools view -v indels subset_hg19.vcf.gz | bcftools query -f '%POS\t%TYPE\n' |wc -l
$ # 141
$ bcftools view -i 'TYPE="indel"' subset_hg19.vcf.gz | bcftools query -f '%POS\t%TYPE\n' | wc -l
```
- 打印所有的基因型genotype(GT)
```bash
$ bcftools query -H -f "%CHROM\t%POS[\t%GT\t]\n" subset_hg19.vcf.gz | head -n 3
```
- 排除样本输出
```bash
$ bcftools view -s ^HG00115,HG00118 subset_hg19.vcf.gz | bcftools query -H -f '%POS[\t%GT\t]\n' | head -n 3
```

- 指定阈值
-c, --min-ac helps to set the minimum allele count of sites to be printed.

```bash
$ bcftools view -c 5 subset_hg19.vcf.gz | grep -v "#" | wc -l 
```

- 选择所有的纯合位点输出
```bash
$ bcftools view -g ^het subset_hg19.vcf.gz | bcftools view -g ^miss | bcftools query -f '%POS[\t%GT\t]\n' | head -n 3
$ # -g, --genotype
$ # ^, exclude  
```

- 选择n个以上的样本含有杂合位点的数据
[SnpSift](https://pcingola.github.io/SnpEff/snpsift/introduction/)
```bash
$ # alias snpsift="java -jar ~/bin/SnpSift.jar"
$ snpsift filter "countHet() >= 3" subset_hg19.vcf.gz | snpsift extractFields - "POS" "GEN[0].GL[1]" | head -n 3
```

- 质量和深度过滤
```bash
$ bcftools query -i 'QUAL>50 && DP>5000' -f '%POS%QUAL%DP\n' subset_hg19.vcf.gz | head -n 3
$ # -i --include
$ # 
```
- 合并vcf或者bcf
```bash
$ bcftools merge -l samplelist > multi-sample.vcf
$ cat samplelist
$ # sample1.vcf.gz
$ # sample2.vcf.gz
$ # sample3.vcf.gz
```

- 找多个vcf的交集
```bash
$ bcftools isec -p outdir -n=3 sample1.vcf.gz sample2.vcf.gz sample3.vcf.gz
$ # isec, intersections
```

- 至少2个以上样本出现的位点
```bash
$ bcftools isec -p outdir -n+2 sample1.vcf.gz sample2.vcf.gz sample3.vcf.gz
```
- 选择只在1个样本中出现的位点
```bash
$ # 只在sample1中出现的位点
$ bcftools isec -p outdir -C sample1.vcf.gz sample2.vcf.gz sample3.vcf.gz
$ -C, --complent
```

- SnpEff
snp注释
```bash
$ # 注意对Java版本有要求，java v1.11以上
$ # 查询数据库
$ java -jar snpEff.jar databases | grep -i GRCh38
$ ../jdk-22/bin/java -jar snpEff.jar download GRCh38.86

$ java -jar snpEff.jar databases | grep -i arabidopsis
```

- ebola
```bash
$ java -jar snpEff.jar databases | grep ebola
$ java -jar snpEff.jar download ebola_zaire
$ java -jar snpEff.jar dump ebola_zaire | more

$ curl -O http://data.biostarhandbook.com/variant/find-variants.sh
$ curl http://data.biostarhandbook.com/sra/ebola-runinfo.csv > runinfo.txt
$ cat runinfo.txt | grep "04-14" | cut -f 1 -d ',' | grep SRR | head -5 > samples.txt
$ bash find-variants.sh KJ660346 samples.txt
$ java -jar snpEff.jar ebola_zaire combined.vcf > annotated.vcf
```
- variants types

|Type|What is means|Example|
| :--- | :--- | :--- |
|SNP	|Single-Nucleotide Polymorphism	|Reference = 'A', Sample = 'C'|
|Ins	|Insertion	|Reference = 'A', Sample = 'AGT'|
|Del	|Deletion	|Reference = 'AC', Sample = 'C'|
|MNP	|Multiple-nucleotide polymorphism	|Reference = 'ATA', Sample = 'GTC'|
|MIXED	|Multiple-nucleotide and an InDel	|Reference = 'ATA', Sample = 'GTCAGT'|


### table
```bash
$ cut -f5 --complement # 去除第五列，输出剩余列

$ # csvtk pretty 使文件对齐，便于查看
$ csvtk pretty runinfo.csv

$ csv2xlsx csv/tsv 转变为Excel格式，便于在windows上查看
$ csvtk csv2xlsx -t HICD2_Kleb_blast_anno.csv -o HICD2_Kleb_blast_anno.xlsx
$ # -t 表示输入的csv文件采用tab分隔

### http
```

### http

```bash
$ # 在Linux中经常要给Windows用户传递文件，可以使用http服务，减少文件的拷贝步骤
$ mkdir result # 创建文件夹，用于共享结果；将需要共享的文件考入该目录
$ cd result；
$ python -m http.server 8000 # 使用Python的简易http服务，共享文件，
$ # windows用户 http://ip:8000/即可网页下载
```

### 命令执行
```bash
$ # 执行file中最后的一行shell命令
$ cat file.txt | tail -n 1 | xargs -I {} bash -c '{}'
$ xargs -I 从标准输入中读取作为参数传递给后续命令
$ bash -c 用于执行shell脚本
```

- 每十行作为一个参数
```bash
$ # xargs -n
$ cut -d "," -f 11 SRP342448.runinfo.csv | sed '1d' | xargs -n 10
```

- java
```bash
$ java -X
$ # -Xmx<size>        set maximum Java heap size
```
### 服务器配置信息
- 查看物理CPU个数
```bash
$ cat /proc/cpuinfo| grep "physical id"| sort| uniq| wc -l
```
- 查看每个物理CPU中core的个数(即核数)
```bash
$ cat /proc/cpuinfo| grep "cpu cores"| uniq
```
- 查看逻辑CPU的个数
```bash
$ cat /proc/cpuinfo| grep "processor"| wc -l
```

- 查看CPU信息（型号）
```bash
$ cat /proc/cpuinfo | grep name | cut -f2 -d: | uniq -c
```

- 查看有多少个虚拟CPU
```bash
$ cat /proc/cpuinfo | grep ^processor | sort | uniq |wc -l
```

- 查看CPU是几个核心的
```bash
$ cat /proc/cpuinfo | grep 'cpu cores' | uniq
```

- 查看集群CPU配置信息
```bash
$ cat /proc/cpuinfo
```
- 查看系统内核
```bash
$ uname -a          # 方法1
$ cat /proc/version # 方法2
```

- 查看内存信息
```bash
$ cat /proc/meminfo
```
- 使用free查看内存信息
```bash
$ free -m
```
- 查看各分区使用情况
```bash
$ df -h
```

- 远程关机
```bash
$ # 方法1
$ shutdown -h now # 系统立即关机
$ shutdown -r now # 系统立即重启
$ shutdown -h +10 # 系统十分钟后关机

# 方法2
$ poweroff   # 关机

# 方法3
$ reboot        # 重启
$ reboot --halt # 停止机器
$ reboot -p     # 关闭机器
```

- 查看Linux系统版本的命令
```bash
$ cat /etc/redhat-release # 方法1
$ lsb_release -a          # 方法2
$ 
```
