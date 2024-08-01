
## Shell Commands
* [帮助信息](#帮助信息)
* [系统命令](#系统命令)
* [系统用户登录信息](#系统用户登录信息)
* [文件夹操作](#文件夹操作)
* [查找文件](#查找文件)
* [用户管理](#用户管理)
* [网络操作](#网络操作)
* [文件处理](#文件处理)
* [文件权限](#文件权限)
* [进程管理](#进程管理)
* [自定义快捷键](#自定义快捷键)


### 帮助信息

```bash
$ man     # 手册(manual) 命令帮助信息查看
$ man man # 查看man的帮助信息
$ man -w man    # --where 显示man文档的位置
$ # /usr/share/man/man1/man.man-db.1.gz
$ man /usr/share/man/man1/man.man-db.1.gz

$ man -w -a man # --all --where
$ man ls        # 查看ls帮助信息
```

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

$ mkdir -p ChIP-seq/{data,results,scripts,software} # 一次性在某个目录下创建多个目录
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

$ wget https://ftp.arb-silva.de/current/SILVA_138.2_LSUParc_tax_silva.fasta.gz --directory-prefix=Exports # 下载到指定的目录

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
$ tail -n+2       # 从文件的第二行开始显示，一般用于去除表头行

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

### 文件压缩与解压

```bash
$ #---- 压缩 ----#
$ # tar 
$ tar -zcvf test.tar.gz ./test/ # 压缩当前文件夹
$ tar -cvf test.tar ./test/     # 打包当前文件夹，但是不压缩

$ tar -jcvf test.tar.bz2 ./test/ # 换一种压缩方式

$ #---- 解压 ----#
$ # tar.gz
$ tar -zxvf xxx.tar.gz

$ # tar.bz2
$ tar -jxvf xxx.tar.bz2

$ # zip
$ unzip xx.zip
```

### 自定义快捷键

```bash
# alias
alias ll='ls -lh'
alias les='less -S'
alias htopz='htop -u zhusitao'
alias utr='cd /home/zhusitao/project/utr'
alias today="date +%F"
alias software="cd /home/zhusitao/software"
alias ldir='ls -d */'
alias lsd='ls -d */ | sed "s#/##g"'
alias q1='awk "NR%4==1"'
alias q2='awk "NR%4==2"'
alias rc='python -c "import sys;from Bio.Seq import Seq;a = [print(Seq(i.strip()).reverse_complement()) for i in sys.stdin];"'
# awk fish
alias c1="awk 'NR==FNR { lines[\$0]=1; next } \$0 in lines'"
alias c4="awk 'FNR==NR {arr[\$4]++; next}{if(\$4 in arr)print \$0}'"
alias fastq_check="sh ~/software/own/fastq_check.sh "
alias wgetd="wget -r -np -nH -R index.html"
```
