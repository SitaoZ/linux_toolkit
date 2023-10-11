# Linux_toolkit
Linux commands and tricks in bioinformatics


## One-liner
```bash
$ history | awk '{a[$2]++} END{for(i in a){print a[i]" "i}}' | sort -rn | head # 列出常用的命令
$ wtach vmstat -sSM     # 实时监控
$ vmstat -sSM           # 监控一次
$ du -h -d 1 | sort -rh # 找出最大文件夹
```

## linux 特殊符号
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

$ # 2>> 错误追缴符号
$ ls -lth 2>>stderr.csv

$ # < 重定向输入符号

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

$ # - 短横线 表示标准输出，或从标准输出中获得输入。一般用于程序需要对个输入的时候。
$ paste <(zcat LJ0612_1.fq.gz|paste - - - -) <(zcat LJ0612_2.fq.gz|paste - - - - ) | awk -v FS="\t" -v OFS="\n" 'FNR==NR {samples[$2]=$1; next} {barcode = substr($2,0,6); if(samples[barcode]) { print $1,$2,$3,$4>>samples[barcode]"_1.fq"; print $5,$6,$7,$8>>samples[barcode]"_2.fq"}}' samples.txt -
$ # 命令解析
$ # awk 命令，samples.txt是第一个输入文件，FNR: 当前文件的行号，NR: 总文件的行号，当有多个文件的是时候，FNR==NR表示第一个文件的内容
$ # (zcat | paste - - - -)  将四行文件转成单行, 括弧表示整体执行，执行后作为标准输入传给 paste
$ # 因此，paste 将两个标准输入合并成一个，read1 和对应的read2 八行变成一行
$ # awk直接输入的samles.txt 作为第一个输入, - 表示将将管道符来源的输入传给当前命令
$ # 所以，管道符号接入的文件作为第二个输入，表示上面FNR==NR表示第一个文件的内容是正确的

$ # $ 变量调用符号

$ # ; 命令分割符号(command separator)，一行中使用;链接的命令顺序执行

$ # ;; 连续分号(terminator)，专门在case语句中承担终止的角色

$ # () 整体执行(command group)

$ # {} 变量分离

$ # [] 中括号流程控制中使用

# head -n 10 < file.fa

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
$ uname -a # 显示系统和内核
$ hostname # 显示当前系统的主机名
$ du       # 显示磁盘使用空间
$ du -sh dir    # 显示文件的总的占用空间
$ du -h -d 1 ./ # 显示当前文件中深度为1 的文件的大小
$ df     # 报告文件系统磁盘空间的使用情况
$ df -h  # 友好查看磁盘
$ mount  # 挂载文件系统
$ umount # 卸载文件系统
$ fsck   # 检查并修复文件系统
$ fdisk  # 磁盘分区命令，适用于2TB以下的磁盘分区
$ parted # 磁盘分区命令，没有磁盘大小限制，常用与2TB一下的磁盘分区
$ swapon # 启用交换分区
$ swapoff # 关闭交换分区
```

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
$ whereis # 给命令找到二进制，源代码和手册文件
$ which   # 查找二进制命令，按环境变量PATH的路径找
$ locate  # 从数据库查找命令
```

### 用户管理
```bash
$ useradd # 添加用户
$ userdel # 删除用户
$ groupadd #  添加用户组
$ passwd   # 修改用户密码
$ change   # 修改用户密码的有效期限
$ su       # 切换用户
$ sudo     # 以另外一个身份执行sudoers文件中允许的命令
$ 
```

### 网络操作
```bash
$ ssh    # 使用SHH加密协议远程登录
$ telnet # 使用TELNET协议远程登录
$ scp    # secure copy 用于不同主机之间复制文件
$ wget   # 下载文件
$ ping   # 测试主机之间网络的联通性
$ ifconfig   # 查看、配置、重启或者禁用网络接口
$ netstat    # 查看网络状态
$ ss         # 查看网络状态
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
$ chgrp                      # 更改文件用户组
$ 
```


### 进程管理

```bash
$ ps         # 查看进程 show snapshoot of processes
$ top        # 查看实时进程 show real time processes
$ kii pid    # 删除进程 kill process with id pid
$ pkill name # 使用程序名称删除进程 kill process with name 
```

## linux三剑客
### awk
```bash
$ awk 
```
### grep

### sed

### cut 


## Conda tips
Conda软件安装十分便利，可以建立不同的环境对软件进行依赖匹配。

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
```

3.数据检查
文件下载后需要对文件进行检查，确保文件在下载中没有出现错误，导致文件不完整。
```
$ # linux 使用md5sum来对文件进行检查
$ md5sum xxxxx       # 对于单个文件，可以执行该命令两次，看产生的md5值是否一致
$ md5sum -c MD5.txt  # 输入MD5文件检查，可以批量检查很多文件
$ # Mac 使用md5 对文件进行检查
```
##  Bio format
### FASTA
FASTA文件的处理
1.seqkit 
seqkit工具可以快速处理fasta文件 [seqkit](https://bioinf.shenwei.me/seqkit/usage/)。
```bash
$ # 按照长度过滤,选取长度小于300bp的fasta子集
$ seqkit seq -M 300 refer.fasta -o lt300.fa
$ # 使用seqkit 选取特定的子集，使用grep子命令 
$ seqkit grep -n -f wanted_gene.csv refer.fasta -o wanted.fa
```
### BAM
SAMtools是li heng开发的用于比对文件处理的利器[samtools](http://www.htslib.org/)。
```bash
$ # 不同比对软件的tag有细微差异，注意区分
$ samtools view QC.sort.bam | grep "XM:i:0" > noMismatch.sam # 找出没有mismatch的比对
$ samtools view QC.sor.bam | grep "AS:" | grep –v "XS:" > unique_alignments.sam # 从bowtie2筛选唯一比对
$ samtools view reads.bam | grep 'XT:A:U' | samtools view -bS -T referenceSequence.fa - > reads.uniqueMap.bam # 从bwa比对文件中筛选唯一比对
$ (samtools view -H QC.sort.bam; samtools view QC.sort.bam | grep -P "\tNH:i:1\t|\tNH:i:1$" | grep -v "ZS:i" ) | samtools view -bS - > unique.bam # 从hisat2中筛选唯一比对
```
