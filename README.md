## Linux_toolkit
Linux commands and tricks in bioinformatics


### One-liner
```bash
$ history | awk '{a[$2]++} END{for(i in a){print a[i]" "i}}' | sort -rn | head # 列出常用的命令
$ wtach vmstat -sSM     # 实时监控
$ vmstat -sSM           # 监控一次
$ du -h -d 1 | sort -rh # 找出最大文件夹
```

### Conda tips
1.创建环境(create)
```bash
$ conda --version          # 显示conda版本
$ conda create -n env_name python=3.7.2 # 创建环境
$ conda activate env_name  # 激活环境
$ conda deactivate         # 退出环境
$ conda env list           # 显示当前所有环境
$ conda info --env         # 显示当前所有环境
```

2. 列出环境中的安装包(list)
```bash
$ conda list             # 显示当前环境中已经安装的包
$ conda list -n env_name # 列出env_name环境中的安装包
$ conda list --export > package-list.txt # 保存安装包便于后续使用
$ conda create -n new_env_name --file package-list.txt # 参考文件重新安装并创建新环境
```

3. 清理安装包的缓存(clean)
```bash
$ conda clean -h
$ conda clean -a # 快速删除
$ conda clean -p # 从可写包缓存中删除没有使用的包，但是不会检查已经安装的包是否有软连接到其中的缓存包
$ conda clean -t # 一键删除anaconda pkgs下面的压缩包

```
4. 包的安装(install)
```bash
$ conda install scipy # 当前环境下安装软件
$ conda install -n env_name scipy # 指定环境下安装软件
```
5. 配置文件(config)
```bash
$ conda config --show         # 显示已经设置好的配置文件的值
$ conda config --describe     # 显示所有可用的配置文件选项
$ conda config --get channels # 显示已经添加的channel
$ conda config add channels x # 添加镜像
$ # 国内提供conda镜像的大学
$ # 清华大学: https://mirrors.tuna.tsinghua.edu.cn/help/anaconda/
$ # 北京外国语大学: https://mirrors.bfsu.edu.cn/help/anaconda/
$ # 南京大学: http://mirrors.nju.edu.cn/
$ # 上海交通大学: https://mirror.sjtu.edu.cn/
$ # 哈尔滨工业大学: http://mirrors.hit.edu.cn/#/home
$ conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/free/
$ conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main/
$ conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge/
$ conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/bioconda/
$ vi ~/.condarc # 直接添加镜像网址也可以
```
### Data download
生信数据庞杂，如何下载自己想要的数据，除了要找对数据库，也得找到合适的工具。
1. SRA tools 下载测序数据
```bash
# 事先安装好SRA toolkit，安装步骤 https://github.com/ncbi/sra-tools/wiki
$ prefetch SRR11180057
$ # prefetch --option-file SraAccList.txt 一次下载多个，使用文件输入
$ fastq-dump --split-files SRR11180057.sra # 将SRA文件转化成fastq
$ # 注意，你也可以一步同时实现下载和解压
$ fastq-dump --split-files SRR11180057 # 除去结尾的.sra后缀
```
具体可以参考 [NCBI sra](https://www.ncbi.nlm.nih.gov/sra/docs/sradownload/)

2. efetch 下载任一格式的数据
Entrez Direct [link](https://www.ncbi.nlm.nih.gov/books/NBK179288/)使用命令行来对NCBI中的数据进行下载，efetch就是其中子程序之一。
```bash
# 以MEDLINE的格式下载文献
$ efetch -db pubmed -id 25359968 -format medline
Retrieve the record for PMID 26287646 in XML format:
# 使用XML下载文献
$ efetch -db pubmed -id 26287646 -format xml
# 以 abstract下载多个文献
$ efetch -db pubmed -id 24102982,21171099,17150207 -format abstract

```
### BAM
```bash
$ samtools view QC.sort.bam | grep "XM:i:0" > noMismatch.sam # 找出没有mismatch的比对
$ samtools view QC.sor.bam | grep "AS:" | grep –v "XS:" > unique_alignments.sam # 从bowtie2筛选唯一比对
$ samtools view reads.bam | grep 'XT:A:U' | samtools view -bS -T referenceSequence.fa - > reads.uniqueMap.bam # 从bwa比对文件中筛选唯一比对
$ (samtools view -H QC.sort.bam; samtools view QC.sort.bam | grep -P "\tNH:i:1\t|\tNH:i:1$" | grep -v "ZS:i" ) | samtools view -bS - > unique.bam # 从hisat2中筛选唯一比对
```
