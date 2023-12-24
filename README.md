# Linux_toolkit
Linux commands and tricks in bioinformatics  
## Table of content
* [One-liner](#One-liner)
* [Linux bash strict model](#Linux-bash-strict-model)
* [Linux symbol](#linux-symbol)
* [shell 常用命令](#shell-常用命令)
  * [系统命令](#系统命令)
  * [系统用户登录信息](#系统用户登录信息)
  * [文件夹操作](#文件夹操作)
  * [查找文件](#查找文件)
  * [用户管理](#用户管理)
  * [网络操作](#网络操作)
    * [远程登陆](#远程登陆)
    * [文件下载和拷贝](#文件下载和拷贝)
    * [网络配置](#网络配置)
    * [网络诊断](#网络诊断)
    * [网络连接](#网络连接)
    * [流量统计](#流量统计)
   * [文件处理](#文件处理)
   * [文件权限](#文件权限)
   * [进程管理](#进程管理)
* [Linux三剑客](#Linux三剑客)
   * [awk](#awk)
     * [awk 内置变量](#awk-内置变量)
     * [awk 运算与判断](#awk-运算与判断)
     * [awk 正则运算](#awk-正则运算)
     * [awk 读取](#awk-读取)
     * [awk 流程控制](#awk-流程控制)
     * [awk 数组应用](#awk-数组应用)
     * [awk 内置函数](#awk-内置函数)
   * [grep](#grep)
     * [pattern syntax](#pattern-syntax)
     * [match control](#match-control)
     * [general output control](#general-output-control)
     * [output line prefix control](#output-line-prefix-control)
     * [context line control](#context-line-control)
     * [grep 正则表达式](#grep-正则表达式)
   * [sed](#sed)
     * [sed 参数](#sed-参数)
     * [动作说明](#动作说明)
     * [print 打印命令](#print-打印命令)
     * [deletion 删除命令](#deletion-删除命令)
     * [substitute 替换标记](#substitute-替换标记)
     * [transform 转换](#transform-转换)
     * [quit 退出](#quit-退出)
* [Vim](#Vim)
* [库文件](#库文件)
* [Conda tips](#Conda-tips)
* [Data download](#Data-download)

## One-liner
### 基本文件处理
```bash
$ history | awk '{a[$2]++} END{for(i in a){print a[i]" "i}}' | sort -rn | head # 列出常用的命令

$ watch vmstat -sSM     # 实时监控
$ vmstat -sSM           # 监控一次

$ du -h -d 1 | sort -rh # 找出最大文件夹

$ cat file.txt | sort | uniq -c | sort -k1nr | head # 排序找出现最多的

$ echo $PATH | tr ":" "\n" | nl # 打印全部路径按行排列

$ sed '/^$/d' file.txt # 去除空白行
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

## Linux bash strict model 
Linux bash strict model非官方模式，和perl `use strict;`类似。
```bash
$ set -e # error exit，当一个未处理的错误出现时立刻跳出程序，不会继续执行。
$ set -u # 当调用没有设置的变量时，报错。
$ set -o pipefail # 设置pipefail时，如果多级管道报错，会导致整个管道程序报错。
```
IFS
```bash
$ IFS # interval Field Separator，分隔符，会影响split函数
$ IFS=$' '
$ items='a b c'
$ for x in $items;
$ do 
$	 echo $x;
$ done 
```

## Linux symbol
linux中一些符号具有特定的含义，需要注意
```bash
$ # # 井号 (comments)

$ # ~ 家目录 (home)

$ # > 重定向符号
$ ls -s > result.txt # 覆盖文件

$ # >> 重定向追加符号
$ date >> date_time.txt

$ # 2> 错误重定向符号
$ ls -l 2>stderr.csv

$ # 2>> 错误追加符号
$ ls -lth 2>>stderr.csv

$ # < 重定向输入符号
$ head -n 10 < file.fa


$ # << 追加重定向输入符号

$ # | 管道符号
$ cat xxx.csv | head -n 10

$ # & 后台进程符号
$ sh work.sh & # 放到后台执行命令

$ # "" 双引号(double quote) 所包含的字符作为普通字符，但是\和$除外

$ # '' 单引号(single quote) 所包含的字符视为普通字符， 无例外

$ # `` 倒引号(backticks) 执行它所包含的命令
$ `date`

$ # \ 倒斜线 1.作为转义字符，将符号特定的含义去掉，使其变成普通标点或者 取消alias的含义。2.放在命令的最末端表示接下一行

$ # - 短横线 表示标准输出，或从标准输出中获得输入。一般用于程序需要多个输入的时候。
$ paste <(zcat LJ0612_1.fq.gz|paste - - - -) <(zcat LJ0612_2.fq.gz|paste - - - - ) | awk -v FS="\t" -v OFS="\n" 'FNR==NR {samples[$2]=$1; next} {barcode = substr($2,0,6); if(samples[barcode]) { print $1,$2,$3,$4>>samples[barcode]"_1.fq"; print $5,$6,$7,$8>>samples[barcode]"_2.fq"}}' samples.txt -
$ # 命令解析
$ # awk 命令，samples.txt是第一个输入文件，FNR: 当前文件的行号，NR: 总文件的行号，当有多个文件的是时候，FNR==NR表示第一个文件的内容
$ # (zcat | paste - - - -)  将四行文件转成单行, 括弧表示整体执行，执行后作为标准输入传给 paste
$ # 因此，paste 将两个标准输入合并成一个，read1 和对应的read2 八行变成一行
$ # awk直接输入的samles.txt 作为第一个输入, - 表示将将管道符来源的输入传给当前命令
$ # 所以，管道符号接入的文件作为第二个输入，表示上面FNR==NR表示第一个文件的内容是正确的

$ # samtools 使用 - 表示输入的文件来至于标准输出
$ hisat2 --phred64 --sensitive --no-discordant \
$        --no-mixed -I 1 -X 1000 \
$        -x hg19/GenomeHisat2Index/chrALL \
$        -1 HBRR1_1.fq.gz -2 HBRR1_2.fq.gz 2>HBRR1.Map2GenomeStat.xls | samtools view -b -S -o HBRR1.bam -

$ # $ 变量调用符号

$ # ; 命令分割符号(command separator)，一行中使用;链接的命令顺序执行

$ # ;; 连续分号(terminator)，专门在case语句中承担终止的角色

$ # () 整体执行(command group)
$ cat <(head -n 1 xx.table) <(cat xx.table | grep -w "^Chr1") > result
$ # 保留表头添加到筛选出的文件中，然后重定向新的文件

$ # {} 变量分离

$ # [] 中括号流程控制中使用


$ # 组合字符
$ $? 状态值 status variable
$ $# 脚本中参数的个数
$ $* 获取所有对应参数的值
$ $0 脚本名
$ $n 第n个参数
$ $@ 获取所有对应的参数值

```

## shell 常用命令

### 系统命令
```bash
$ uname -a      # 显示系统和内核
$ hostname      # 显示当前系统的主机名
$ du            # 显示磁盘使用空间
$ du -sh dir    # 显示文件的总的占用空间
$ du -h -d 1 ./ # 显示当前文件中深度为1 的文件的大小
$ df            # 报告文件系统磁盘空间的使用情况
$ df -h         # 易读取查看磁盘
$ mount         # 挂载文件系统
$ umount        # 卸载文件系统
$ fsck          # 检查并修复文件系统
$ fdisk         # 磁盘分区命令，适用于2TB以下的磁盘分区
$ parted        # 磁盘分区命令，没有磁盘大小限制，常用与2TB一下的磁盘分区
$ swapon        # 启用交换分区
$ swapoff       # 关闭交换分区
```

```bash
$ ldd # 查看文件动态库依赖 print shared object dependencies
$ ldd /bin/ls # 查看文件的动态依赖
```
```bash
$ ulimit    # 控制shell程序资源
$ ulimit -a # 显示目前资源的限定
$ ulimit -n # 指定同一时间最多可开启的文件数
$ ulimit -H # 设定资源的硬性限制，也就是管理员所设下的限制

$ launchctl limit maxfiles                   # macOS 最大文件数目
$ sudo launchctl limit maxfiles 65536 200000 # 将文件句柄数设置到最大 
```

- crontab
命令可以在固定的间隔时间执行指定的系统指令或 shell script脚本。
时间间隔的单位可以是分钟、小时、日、月、周及以上的任意组合。这个命令非常适合周期性的日志分析或数据备份等工作。
```bash
$ cat /etc/crontab
SHELL=/bin/bash
PATH=/sbin:/bin:/usr/sbin:/usr/bin
MAILTO=root

# For details see man 4 crontabs

# Example of job definition:
# .---------------- minute (0 - 59)
# |  .------------- hour (0 - 23)
# |  |  .---------- day of month (1 - 31)
# |  |  |  .------- month (1 - 12) OR jan,feb,mar,apr ...
# |  |  |  |  .---- day of week (0 - 6) (Sunday=0 or 7) OR sun,mon,tue,wed,thu,fri,sat
# |  |  |  |  |
# *  *  *  *  * user-name  command to be executed

$ # 编辑文件，按照上面的额格式添加任务
$ vi /etc/crontab 
$ # 启动服务
$ /bin/systemctl start crond.service

$ # 非root用户想使用crontab
$ # 首先创建一个cron.allow文件,将用户名写入
$ touch /etc/cron.allow

```
[crontab](https://linuxtools-rst.readthedocs.io/zh-cn/latest/tool/crontab.html)
### 系统用户登录信息
```bash
$ whoami  # 显示当前用户的名称
$ who     # 显示目前系统的用户信息
$ w       # 显示已经登录的用户列表，并显示用户正在执行的指令
$ last    # 显示登录系统的用户
$ lastlog # 显示当前系统中所有用户最近一次的登陆信息
$ users   # 显示当前登录系统的所有用户
```

### 文件夹操作
```bash
$ pwd       # 显示当前文件夹 print working dir
$ mkdir dir # 创建目录 make directory 
$ cd dir    # 进入dir目录 change directory to dir 
$ cd ..     # 进入上一个层级目录 go up a directory 
$ ls        # 列出文件 list files
```
### 查找文件
```bash
$ find  dir -name *fasta # 在dir 目录下查找后缀为fasta的文件
$ whereis                # 给命令找到二进制，源代码和手册文件
$ which                  # 查找二进制命令，按环境变量PATH的路径找
$ locate                 # 从数据库查找命令
```

### 用户管理
```bash
$ useradd  # 添加用户, 注意Ubuntu命令不一致, adduser
$ userdel  # 删除用户
$ groupadd # 添加用户组
$ passwd   # 修改用户密码
$ change   # 修改用户密码的有效期限
$ su       # 切换用户
$ sudo     # 以另外一个身份执行sudoers文件中允许的命令
$ 
```

### 网络操作

#### 远程登陆 
SSH(secure shell protocol)安全外壳协议，一种加密的网络传输协定
```bash
$ ssh    # 使用SHH加密协议远程登录
$ man ssh
$ ssh -p 22 username@ip # 远程指定端口登陆
$
$ telnet # 使用TELNET协议远程登录
```

#### 文件下载和拷贝
```bash 
$ scp    # secure copy 用于不同主机之间复制文件
$ wget   # 下载文件
$ wget http://data.biostarhandbook.com/data/sequencing-adapters.fa # sequencing adapter
$ wget http://www.example.com/filename.txt -o /path/filename.txt
$ wget -c http://example.com/samplefile.tar.gz # -c 恢复下载文件
$ wget -nc -P ~/refs http://geneontology.org/ontology/go-basic.obo
$ # -nc/--no-clobber 文件只下载一次，如果文件已经下载在~/refs目录，则不需要重复下载
$ # -P/--prefix-directory 目录前缀， 默认是当前目录

$ wget -r -np -nH -R index.html  http://data.biostarhandbook.com/circos/config/ # 下载目录和其中全部文件
$ # -r/--recursive 遍历所有子目录
$ # -np/--no-parent 不到父目录去
$ # -nH/--no-host-directories 不要将文件保存到主机名文件夹
$ # -R/--reject index.html 不下载index.html文件

$ curl   # (client URL) transfer a URL 用于数据传输，支持各种传输协议
$ man curl
$ curl -# # --progress-bar 显示进度条
$ curl -O # --remote-name 使用远程文件的文件名作为写入文件名
$ curl -O http://data.biostarhandbook.com/rnaseq/mouse-gene-expression.txt
```
#### 网络配置
```bash
$ ifconfig # 查看、配置、重启或者禁用网络接口
$ ip 
```

#### 网络诊断
```bash
$ ping   # 测试主机之间网络的联通性
$ mtr    # mytraceroute，它集成了 ping、 traceroute、 nslookup 的功能，诊断网络问题非常方便
$ mtr -n www.baidu.com
```
#### 网络连接
```bash 
$ netstat       # 查看网络状态
$ netstat -ntpl # 查询TCP类型端口
$       # -n/--numeric 显示数字地址
$       # -t/--tcp     显示tcp, transmission control protocol 传输控制协议
$       # -p/--program 显示占用的程序
$       # -l/--listen  正在监听
$ netstat -nupl # 查询UDP 端口类型
$ ss            # 查看网络状态
```

[TCP/UDP/PORT](https://zhuanlan.zhihu.com/p/57987304)

#### 流量统计
```bash
$ man ifstat
$ man sar
$ 
```

### 文件处理
```bash
$ touch filea     # 创建文件 create file 
$ cat file1 file2 # 合并文件1和2 concatenate files and output
$ less file1      # 查看文件
$ file file1      # 查看文件类型
$ cp file1 file2  # 复制文件1到文件2 copy file1 to file2
$ mv file1 file2  # 重命名文件 move file1 to file2
$ rm file1        # 删除文件
$ head file1      # 前十行
$ tail file1      # 后十行
$ tail -F file1   # 实时查看文件
$ sort            # 对文件排序
$ uniq            # 对文件去重
$ wc              # 统计文件行数，单词数，字符数
$ diff file1 file2# 比较文件差异
$ paste           # 按行合并文件内容
$ split           # 将分割文件成小部分
$ tr              # 替换或删除字符
$ rev             # 反向输出内容
$ echo agctagtcg | tr a-z A-Z | tr ATCG TAGC | rev # 反向互补DNA序列

```

### 文件权限
```bash
$ chmod 755 file1            # change mode of file
$ chomod -R zhusitao folder  # recursively chmod folder to zhusitao
$ chown zhusitao ath.csv     # 将文件的所属者改成zhusitao
$ chown user:group file1     # 将文件的所属者改成user, 将文件的所属组改成group
$ chown :gruop file1         # 只改变所属组, 不改变所有者
$ chgrp -v group_id file1    # 更改文件用户组
$ 
```


### 进程管理

```bash
$ ps         # 查看进程 show snapshoot of processes
$ top        # 查看实时进程 show real time processes
$ kii pid    # 删除进程 kill process with id pid
$ pkill name # 使用程序名称删除进程 kill process with name 
```

## Linux三剑客
### awk
文本处理的工具之一,[awk](https://wangchujiang.com/linux-command/c/awk.html)
```bash
$ man awk
$ # -f/--file 从文件中读取awk脚本执行,是命令行执行wak的一种补充方式，可以使用多个-f
$ # -F/--field-separator fs，指定输入文件分隔符号
$ # -v/--assign 指定变量和其对应的值，

$ awk执行方式
$ awk [options] 'script' file(s)
$ awk [options] -f scriptfile files(s)

$ awk 'BEGIN{ print "start" } pattern{ commands } END{ print "end" }' file # 三个语句块都是可选的

$ awk 模式和操作 # pattern {action statement}
$ 模式(pattern)主要包括: /正则表达式/: 使用通配符;
$            关系表达式: 使用运算符进行操作;
$            模式匹配表达式: ~匹配 和 !~不匹配;
$            BEGIN语句、pattern语句、END语句;
$ 操作(action)主要包括：有一个或者多个命令、函数、表达式组成，用换行符和分号隔开，并位于大括号内。包括变量数组赋值，输出命令，内置函数和控制语句


```

#### awk 内置变量
```bash
$ ARGC 命令行参数数目
$ ARGV 包含命令行参数的数组
$ FILENAME 当前文件名
$ FNR 当前文件的行号，记录数
$ FS 字符分割符号，默认是空格 Field Separator
$ OFS 输出字段分隔符号，默认是一个空格, Output Record Separator
$ RS 记录分割符, 默认是一个换行符
$ ORS 输出记录分隔符，默认是一个换行符
$ NF 字段数 Number of Fields
$ NR 文件行号，记录数
$ awk 传入外部变量
$ awk -v variable=100 '{print $variable}' files
```

#### awk 运算与判断
```bash
$ awk 算术运算符
$ # + - 加减
$ awk '{print $1-1, $1+2}' files
$ * / % 乘 除求余
$ awk "NR%4==2" fastq_files
$  ^**n 求幂
$ awk '{print 2**$1}' files
$ ++ -- 自加和自减
$ awk

$ akw 赋值运算符
$ =, +=, -=, *=, /=, %=, ^=, **=

$ 逻辑运算符
$ || 逻辑或
$ && 逻辑与
$ awk 'BEGIN{a=1;b=2;print (a>5 && b<=2),(a>5 || b<=2);}'
```
#### awk 正则运算
```bash 
$ 正则运算符
$ ~, !~
$ ^ 行首
$ $ 行尾
$ . 除了换行符以外的任意单个字符
$ * 前导字符的零个或者多个
$ .* 所有字符
$ [] 字符组内的任一字符
$ [^] 不匹配字符组内的任一字符
$ ^[^] 非字符组内的字符开头的行
$ [a-z] 小写字母
$ [A-Z] 大写字母
$ [a-Z] 小写和大写字母
$ [0-9] 数字
$ \< 单词头
$ \> 单词尾
$ 正则需要用 /正则/ 包围

$ 关系运算符
$ <, <=, >, >=, !=, ==
```
#### awk 读取
```bash
$ # next 语句
$ awk 'NR%4==3{next}{print NR,$0}' fastqs # 去读fastq文件，判断去除其第三行
$ next 语句，循环逐行匹配，匹配就跳过，然后进行下一行匹配，next一般用于多行合并
$ cat text.txt
$ # >chr1
$ # AAAAAAAAA
$ # >chr2
$ # TTTTTTTTT
$ # >chr3
$ # CCCCCCCCC
$ # >chr4
$ # GGGGGGGGG
$ awk '/^>/{ID=$0;next;}{print ID","$0}' text.txt
$ # >chr1,AAAAAAAAA
$ # >chr2,TTTTTTTTT
$ # >chr3,CCCCCCCCC
$ # >chr4,GGGGGGGGG

$ # getline 从标准输入、管道和正在处理的文件之外的其他输入文件或得输入
$ awk 'BEGIN{ "date" | getline out; print out }' test
$ getline 当其左右没有重定向符|或<时：getline作用于当前文件
$ getline 当其左右有重定向符|或<时：getline作用于输入文件，

$ awk 'BEGIN{ "date" | getline out; split(out,mon); print mon[2] }'
$ seq 5 | awk 'BEGIN{getline; print "The first line "$0};{print $0}'
```

#### awk 流程控制
```bash
$ while, for, do-while; break, continue语句控制流程
$ break 退出循环; continue 中断当前正在执行的循环并跳到下一循环。

$ 条件语句
if(表达式)
  {语句1}
else if(表达式)
  {语句2}
else
  {语句3}

$ while循环语句
while(表达式)
  {语句}
$ awk 'BEGIN{
test=100;total=0;
while(i<=test){
  total += i;
  i++;
}
print total;
}'

$ for循环
$ awk 'BEGIN{
total=0;
for(i=0;i<=100;i++){
  total+=i;
}
print total;
}'
```

#### awk 数组应用
```bash
$ # awk 数组应用十分便利。数组不必提前声明
$ # 数组下标1-based
$ # Array[1] = 'AAAA'
$ # split 分割函数 and length 函数
$ awk 'BEGIN{info="it is a test";lens=split(info,tA," ");print length(tA),lens;}'
```
#### awk 内置函数
```bash
$ 1.1 awk 算术函数
$ # atan2, cos, sin, exp, log, sqrt, int, rand, strand

$ 1.2 awk 字符串函数
$ # gsub, sub, index, length, blength, substr, match, split, tolower, toupper, sprintf

$ 1.2.1 gsub # g 表示 global, target string 可选
$ gsub(/^>/,   "chr",  "Chr")
$ #    \_/      \_/     \_/
$ #     |        |       |
$ #  regular     |    target string  
$ # expression   |
$ #           replace string
$
$ awk 'BEGIN{info="this is a test2010test!";gsub(/[0-9]+/,"!",info);print info}' # gsub 替换

$ 1.2.2 substr
$ substr($3,  1,   2)
$ #      \_/  \_/ \_/
$ #       |    |   |
$ #     string | end position (optional)
$ #            |
$ #         start position
$
$ awk 'BEGIN{info="this is a test2010test!";print substr(info,4,10);}'  #截取字符串， 从第四个字符串开始，，截取10个长度字符串


$ 1.2.3 split 
$ awk '{split($0, array, ":")}'
$ #           \_/  \___/  \_/
$ #            |     |     |
$ #         string   |   delimiter
$ #                  |
$ #               array to store the pieces
$
$ awk 'BEGIN{info="this is a test";split(info,tA," ");print length(tA);for(k in tA){print k,tA[k];}}' # 字符串分割


$ 1.2.4 match
$ match(string, regexp)
$
$ awk 'BEGIN{info="this is a test2010test!";print match(info,/[0-9]+/)?"ok":"no found";}' # match 正则表达式匹配查找

$ 1.2.5 index
$ index($0, "target string")
$
$ awk 'BEGIN{info="this is a test2010test!";print index(info,"test")?"ok":"no found";}' # index 查找字符串
 






$ 1.3 一般函数
$ # close, system, getline,
```

### grep
grep(global search regular expression and print out the line),全面搜索正则表达式并把行打印。
```bash
$ man grep 
```
#### pattern syntax 
```bash
$ -E/--extended-regexp # 拓展的正则匹配
$ -F/--fixed-strings   # 将模式固定，不支持正则表达
$ -G/--basic-regexp    # 一般的正则表达
$ -P/--perl-regexp     # perl正则
```

#### match control 

```bash
$ -e/--regexp # 指定字符串作为模式
$ -f/--file   # 从文件中一行行读取作为模式
$ -i//--ignore-case # 忽略大小写
$ --no-ignore-case  # 不忽略大小写
$ -v/--invert-match # 反向匹配
$ -w/--word-regexp  # 匹配单词
$ -x/--line-regexp  # 匹配整行
```

#### general output control
```bash
$ -c/--count # 打印出匹配的行数， 而不是输出行
$ --color    # 显示颜色
$ -L/--files-without-match # 输出没有匹配的文件名
$ -l/--files-with-matches  # 输出匹配的文件名
$ -m/--max-count     # 指定最大的匹配次数
$ -o/--only-matching # 只打印出匹配的,同时显示文件名
```

#### output line prefix control
```bash
$ -b/--byte-offset   # 输出匹配的位置，0-base开始；可以和-o配合使用
$ -H/--with-filename # 打印出每个匹配的文件名
$ -h/--no-filename   # 不打印文件名
$ --lable            # 转化输入文件成新的文件名，
$ cat text.txt | grep --label=text 'aaa' -H -i
$ -n/--line-number   # 打印行号，1-based
$ -T/--initial-tab   # 在匹配的文件前使用一个tab符号, 结合-H, -n, -b效果更佳
```

#### context line control 

```bash
$ -A/--after-context  # 打印匹配行后面NUM行
$ -B/--before-context # 打印匹配文件前面NUM行
$ -C/--context        # 打印前后NUM行
```

#### grep 正则表达式
```bash
$ ^   # 锚定行的开始 
$ $   # 锚定行的结束
$ .   # 匹配一个非换行字符
$ *   # 匹配零个或者多个先前字符
$ .*  # 任意字符
$ []  # 指定范围内的字符 '[Gg]rep'匹配Grep和grep
$ [^] # 匹配一个不在指定范围内的字符，'[^A-Z]rep'匹配不包含A-Z开头的行
$ \<  # 锚定单词开始
$ \>  # 锚定单词结束
$ x\{m\}     # 重复字符x, 'A\{9\}' 匹配9个A
$ x\{m,\}\   # 重复字符x, 至少m次
$ x\{m, n\}} # 重复字符x, 至少m次，至多n次
$ \w         # 匹配文字或数字字符 [a-zA-Z0-9]
```

### sed
sed (stream editor for filtering and transforming text)流式编辑器
```bash
$ # 命令格式
$ sed [options] 'command' file(s)
$ sed [options] -f scriptfile file(s)
$ # 命令必须由单引号包住；或者使用双引号包住，双引号主要在传入变量时使用
$ # 每个command 由最多两个地址(addresses)和一个动作(action)组成,
$ # 每个地址可以是正则表达式或者行数，动作见4.3.2
```
#### sed 参数

```bash
$ # 常用的sed参数，用来指定执行的方式
$ -e/--expression # 传入脚本到命令行执行
$ -f/--file       # 从文件传入命令执行
$ -l/--line-length# 输出行的字符长度
$ -i/--in-place   # 原地编辑文件
$ -n/--quiet/--silent  # 仅仅显示脚本处理后的结果
$ -E/--regexp-extended # 正则表达式
$ -V/--version # 显示版本信息
$ -h/--help    # 显示帮助信息
```

#### 动作说明

a: 新增 append  
c: 取代 change  
d: 删除 delete  
i: 插入 insert  
p: 打印 print, 通常和-n一起使用，打印特定的行
s: 取代 substitute  
y: 转换 transform  
q: 退出 quit  


#### append 添加新行

```bash
$ sed '/ATCG/a atcg' file.txt # 在匹配的行后面添加新行
```

#### deletion 删除命令
```bash
$ d
$ sed '/^$/d' file # 去除所有的空白行
$ sed '2d'         # 删除第二行
$ sed '$d'         # 删除最后一行
$ sed '2,$d' file  # 删除第二行到末尾所有行
$ sed '1,/^$/d'    # 从第一行到第一个空白行之间的都删除, 这里的/^$/表示锁定第一个空白行
$ sed '/^$/, $d'   # 删除第一个空白行到最后一行的内容
$ sed '/^$/, 10d'  # 删除第一个空白行到第十行之间的内容
$ sed '/^ya*y/, /[0-9]$/d' # 删除以ya*y模式开始的行到第一个以数字结尾的行
$ sed '/^test/d'           # 删除文件中所有开头是test的行
$ echo this is a test line | sed 's/\w\+/[&]/g' # \w+表示匹配每一个单词，使用[&]替换它，& 对应之前所匹配到的单词

```

#### print 打印命令
```bash
$ # Nth line
$ sed -n '2p' txt     # 打印文件第二行
$ sed -n '$p' txt     # 打印文件最后一行
$ sed -n '1,5p' txt   # 打印第一行到第五行
$ sed -n '/chr/p' txt # 打印含有chr的行
$ sed -n '/add/,/sub/p' txt # 打印匹配add开始和匹配sub结束的行
$ sed -n '2,/pat/p' txt     # 打印第二行开始,/pat/行结束的内容
$ sed '/pat/,$d' txt        # 打印/pat/开始到最后一行的内容
$ yes | head -c100 | tr '\n' ' ' | sed -n l | head -n1 | wc -c
```

#### substitute 替换标记
```bash
$ s
$ g # 表示全局替换
$ p # 表示打印行
$ w # 表示把行写入一个文件
$ x # 表示互换模板中的文本和缓冲区的文本
$ y # 表示把第一个字符翻译成另外的字符
$ \1 # 子串匹配标记
$ & # 已匹配字符串标记

$ sed 's/AAAA/TTTT/g' file # 全局替换
$ # 定界符号一致均可
$ sed 's#AAAA#TTTT#g' file 
$ sed 's:AAAA:TTTT:g' file

$ sed '2s/AAAA/TTTT/' file      # 只替换第二行
$ sed '/DDDD/s/AAAA/TTTT/' file # 替换匹配DDDD的行

$ # 正则
$ echo aaa BBB | sed 's/\([a-z]\+\) \([A-Z]\+\)/\2 \1/'
$ echo this is digit 7 in a number | sed 's/digit \([0-9]\)/\1/'

```

#### transform 转换
```bash
$ echo ATCG | sed 'y/ATCG/TAGC/' | rev # DNA反向互补
```

#### quit 退出 
```bash
$ sed '100q' file # 打印前100行，然后退出
```
#### 组合多个表达式
```bash
$ sed '表达式1; 表达式2' # sed '表达式1' | sed '表达式2'
$ sed 's/chr/Chr/g; s/geneid/GeneID/g' xxx.fa # 同时多个替换操作
```

#### 引用

```bash
$ # 传入shell变量时，需要使用双引号
$ gene_id=AT1G79550
$ cat file.txt | sed 's/$gene_id/PGK/g' # 替换基因名
```

#### 多点编辑
```bash
$ nl /etc/passwd | sed -e '3,$d' -e 's/bash/blueshell/' # -e表示多点编辑，第一个编辑命令删除第三行到末尾的数据，第二条命令搜索bash替换为blueshell
```
[sed lecture](https://cs.nyu.edu/~mohri/unix08/lect5.pdf)


## Vim
Vim编辑器

```bash
:%# #\r#g # 空白换成换行符，在vim命令行模式下
:2,19s/aaa/bbb/g # 指定行号之间替换
:g/^$/d   # 把空白行删除，
:1        # 跳到文件首行
:$        # 跳到文件的尾行
:5p       # 跳转到第5行，然后显示第5行内容
:3d       # 跳转到第3行并删除第3行
: set nu  # 显示行号
: noh     # 除去高亮
```

## 库文件
库，其实就是把多个源文件（.c文件），经过一定的翻译，然后打包————到最后只提供给我们一个文件。
Linux库文件是可执行的公共代码，是为了减少程序开发的重复劳动，分为静态库(.a)的动态库(.so)。
静态库文件后缀为.a，程序编译时，将静态文件中的代码复制，拷贝到程序生成的可执行文件中。
静态库的优点就是执行程序时不需要外部的函数库的支持，缺点是如果静态函数库变了，程序必须重新编译。
动态库文件后缀为.so (Shared Object)，动态库在编译是没有被编译进目标代码，而是程序执行到相关代码时才调用库中的函数。


### 1. 动态库 
- .so 结尾
```bash
$ # 查看某个可执行文件的动态库依赖
$ ldd /bin/ls
linux-vdso.so.1 (0x00007fff0e4d5000)
libselinux.so.1 => /usr/lib64/libselinux.so.1 (0x00007fcaf947e000)
libcap.so.2 => /usr/lib64/libcap.so.2 (0x00007fcaf9474000)
libc.so.6 => /usr/lib64/libc.so.6 (0x00007fcaf926b000)
libpcre2-8.so.0 => /usr/lib64/libpcre2-8.so.0 (0x00007fcaf91cf000)
/lib64/ld-linux-x86-64.so.2 (0x00007fcaf94d1000)
```
- so文件的命名规则
```bash
$ libname.so.x.y.z
$  \___/     | | |
$    |       | | release version
$  libname   | minor version
$      major version
$ # major version表示重大升级，不同major version之间的库是不兼容的。
$ # minor version表示增量更新，一般是增加了一些新接口，原来的接口不变。
$ # release version表示库的一些bug修复，性能改进等，不添加任何新的接口，不改变原来的接口。
```

- 动态库查找
```bash
$ # 动态链接库的查找先后顺序为：
$ # (1) LD_LIBRARY_PATH环境变量中的路径
$ # (2) /etc/ld.so.cache缓存文件
$ # (3) /usr/lib和/lib
```

- gcc 是GNU Compiler Collection，原名为Gun C语言编译器，因为它原本只能处理C语言，但gcc很快地扩展，包含很多编译器（C、C++、Objective-C、Ada、Fortran、 Java），可以说gcc是GNU编译器集合；
- g++既可以处理C/C++语言，而gcc只能处理C语言；一般我们使用g++即可；
- gcc/g++就是将包含了代码的文本文件编译（预处理、编译、汇编、链接）成可执行的文件。
- gcc编译选项
使用gcc编译链接时，有两个参数需要注意，一个是-l（小写的L），一个是-L（大写的L）。
Linux有个约定速成的规则，假如库名是name，那么动态链接库文件名就是libname.so。在使用GCC编译链接时，-lname来告诉GCC使用哪个库。
链接时，GCC的链接器ld就会前往LD_LIBRARY_PATH环境变量、/etc/ld.so.cache缓存文件和/usr/lib和/lib目录下去查找libname.so。
我们也可以用-L/path/to/library的方式，让链接器ld去/path/to/library路径下去找库文件。
```bash
gcc -L/path/to/library -lname myfile.c
```
## Conda tips
Conda软件安装十分便利，可以建立不同的环境对软件进行依赖匹配。

1. 创建环境(create)
```bash
$ conda --version                       # 显示conda版本
$ conda create -n env_name python=3.7.2 # 创建环境
$ conda activate env_name               # 激活环境
$ conda deactivate                      # 退出环境
$ conda remove -n env_name --all         # 删除环境

$ conda create --clone env_name1 -n env_name2 #克隆环境到新环境
```

2. 环境(environment) 

```bash
$ conda env list                             # 显示当前所有环境
$ conda info --env                           # 显示当前所有环境
$ conda install -n env_name pkg1 pkg2        # 在指定的环境中安装多个包
$ conda uninstall pkg1 -n env_name           # 在指定的环境下删除包
$ conda env remove -n env_name               # 删除环境
$ conda rename -n env_name env_new_name      # 更改环境名称
$ conda remove -n env_name -c channel_name pkgname # 删除指定环境下，特定来源的包

$ conda env export -n env_name > ENV.yml               # 导出环境信息
$ conda env create -n env_name --file ENV.yml          # 安装指定的环境
```

3. 列出环境中的安装包(list)
```bash
$ conda list                                           # 显示当前环境中已经安装的包
$ conda list --show-channel-urls                       # 列出安装包和来源信息
$ conda list -n env_name                               # 列出env_name环境中的安装包
$ conda list -n env_name --show-channel-urls           # 列出指定环境下安装的全部包和其来源

$ conda list --export > package-list.txt               # 保存安装包便于后续使用
$ conda create -n new_env_name --file package-list.txt # 参考文件重新安装并创建新环境


```

4. 清理安装包的缓存(clean)
```bash
$ conda clean -h
$ conda clean -a # 快速删除
$ conda clean -p # 从可写包缓存中删除没有使用的包，但是不会检查已经安装的包是否有软连接到其中的缓存包
$ conda clean -t # 一键删除anaconda pkgs下面的压缩包

```
5. 包的安装(install)
```bash
$ conda install scipy                         # 当前环境下安装软件
$ conda install -n env_name scipy             # 指定环境下安装软件
$ conda install -c channel_name scipy         # 在指定的环境下安装
$ conda install "pkgname>2.7,<3.5"            # 安装指定的版本
$ conda install "pkgname [version='2.5|3.2']" # 安装指定的版本
$ conda uninstall pkgname
$
$ conda search pkg --info                     # 搜索包的信息
```
6. 配置文件(config)
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
## Data download
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
具体可以参考 [NCBI sra](https://www.ncbi.nlm.nih.gov/sra/docs/sradownload/)。

2. Entrez Direct [link](https://www.ncbi.nlm.nih.gov/books/NBK179288/)使用命令行来对NCBI中的数据进行下载，efetch就是其中子程序之一。

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


3.数据检查
文件下载后需要对文件进行检查，确保文件在下载中没有出现错误，导致文件不完整。
```
$ # linux 使用md5sum来对文件进行检查
$ md5sum xxxxx       # 对于单个文件，可以执行该命令两次，看产生的md5值是否一致
$ md5sum -c MD5.txt  # 输入MD5文件检查，可以批量检查很多文件
$ # Mac 使用md5 对文件进行检查
```
