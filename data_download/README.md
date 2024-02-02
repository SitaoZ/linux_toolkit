## Table of content 
* [SRAtools 下载测序数据](#SRAtools)
* [Entrez Direct](#EntrezDirect)
* [数据检查](#数据检查)

生信数据庞杂，如何下载自己想要的数据，除了要找对数据库，也得找到合适的工具。
### SRAtools 下载测序数据
```bash
# 事先安装好SRA toolkit，安装步骤 https://github.com/ncbi/sra-tools/wiki
$ prefetch SRR11180057
$ # prefetch --option-file SraAccList.txt 一次下载多个，使用文件输入
$ fastq-dump --split-files SRR11180057.sra # 将SRA文件转化成fastq
$ # 注意，你也可以一步同时实现下载和解压
$ fastq-dump --split-files SRR11180057 # 除去结尾的.sra后缀
```
具体可以参考 [NCBI sra](https://www.ncbi.nlm.nih.gov/sra/docs/sradownload/)。

### EntrezDirect 
[Entrez Direct](https://www.ncbi.nlm.nih.gov/books/NBK179288/)使用命令行来对NCBI中的数据进行下载，efetch就是其中子程序之一。

```bash
# efetch 下载任意格式的数据
# 以MEDLINE的格式下载文献
$ efetch -db pubmed -id 25359968 -format medline
# 使用XML下载文献
$ efetch -db pubmed -id 26287646 -format xml
# 以 abstract下载多个文献
$ efetch -db pubmed -id 24102982,21171099,17150207 -format abstract


$ esearch -db sra -query PRJNA830912 | efetch -format summary > summary.xml # 根据项目编号下载
$ # 下载GEO数据库
$ esearch -db gds -query GSE201349 | efetch > GSE201349.txt
$ esearch -db sra -query GSE201349 |efetch -format docsum | xtract -pattern DocumentSummary -element Run@acc

$ # 下载SRA run的信息
$ esearch -db sra -query PRJNA347885 | efetch --format runinfo > GSE87822_runinfo.csv
$ cut -d ',' -f1 GSE87822_runinfo.csv | grep SRR | xargs fastq-dump

$ esearch -db sra -query PRJNA257197 | efetch -format docsum > docsum.xml
$ cat docsum.xml | xtract -pattern DocumentSummary -element Bioproject,Biosample,Run@acc | head -n 5
$ # PRJNA257197	SAMN03253746	SRR1972976
$ # PRJNA257197	SAMN03253745	SRR1972975
$ # PRJNA257197	SAMN03253744	SRR1972974
$ # PRJNA257197	SAMN03254300	SRR1972973
$ # PRJNA257197	SAMN03254299	SRR1972972

$ bash sra-runinfo.sh 2022 10 # Downloading SRA run info for Year=2022 Month=10

$ #运用上述脚本，下载2007-2020的测序数据，结合csvkit (pip install csvkit) 统计不同条件下的频数
$ cat allruns.csv | csvcut -c 29 | sort | uniq -c | sort -rn | head # 物种测的数据频数
$ cat allruns.csv | csvcut -c 19 | sort | uniq -c | sort -rn | head # 测序公司频数
$ cat allruns.csv | csvcut -c 20 | sort | uniq -c | sort -rn | head # 测序仪频数
```
[NBK25497, database name](https://www.ncbi.nlm.nih.gov/books/NBK25497/table/chapter2.T._entrez_unique_identifiers_ui/?report=objectonly)

[NCBI edirect guide](https://www.nlm.nih.gov/dataguide/eutilities/utilities.html)

[NCBI esearch guide](https://www.nlm.nih.gov/dataguide/edirect/esearch.html)

[stowers institute example](https://research.stowers.org/cws/CompGenomics/Tutorial/geo.html)


### 数据检查
文件下载后需要对文件进行检查，确保文件在下载中没有出现错误，导致文件不完整。
```
$ # linux 使用md5sum来对文件进行检查
$ md5sum xxxxx       # 对于单个文件，可以执行该命令两次，看产生的md5值是否一致
$ md5sum -c MD5.txt  # 输入MD5文件检查，可以批量检查很多文件
$ # Mac 使用md5 对文件进行检查
```
