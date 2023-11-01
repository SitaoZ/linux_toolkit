# Linux_toolkit
Linux commands and tricks in bioinformatics
目录
===
* [1. One-liner](#1. One-liner)
* [2. linux 特殊符号](#2. linux 特殊符号)


## 1. One-liner
```bash
$ history | awk '{a[$2]++} END{for(i in a){print a[i]" "i}}' | sort -rn | head # 列出常用的命令
$ wtach vmstat -sSM     # 实时监控
$ vmstat -sSM           # 监控一次
$ du -h -d 1 | sort -rh # 找出最大文件夹
```

## 2. linux 特殊符号
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

## 3. shell 常用命令

### 3.1 系统命令
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

### 3.2 系统用户登录信息
```bash
$ whoami  # 显示当前用户的名称
$ who     # 显示目前系统的用户信息
$ w       # 显示已经登录的用户列表，并显示用户正在执行的指令
$ last    # 显示登录系统的用户
$ lastlog # 显示当前系统中所有用户最近一次的登陆信息
$ users   # 显示当前登录系统的所有用户
```

### 3.3 文件夹操作
```bash
$ pwd       # 显示当前文件夹 print working dir
$ mkdir dir # 创建目录 make directory 
$ cd dir    # 进入dir目录 change directory to dir 
$ cd ..     # 进入上一个层级目录 go up a directory 
$ ls        # 列出文件 list files
```
### 3.4 查找文件
```bash
$ find  dir -name *fasta # 在dir 目录下查找后缀为fasta的文件
$ whereis # 给命令找到二进制，源代码和手册文件
$ which   # 查找二进制命令，按环境变量PATH的路径找
$ locate  # 从数据库查找命令
```

### 3.5 用户管理
```bash
$ useradd # 添加用户, 注意Ubuntu命令不一致, adduser
$ userdel # 删除用户
$ groupadd #  添加用户组
$ passwd   # 修改用户密码
$ change   # 修改用户密码的有效期限
$ su       # 切换用户
$ sudo     # 以另外一个身份执行sudoers文件中允许的命令
$ 
```

### 3.6 网络操作

#### 3.6.1 远程登陆 
SSH(secure shell protocol)安全外壳协议，一种加密的网络传输协定
```bash
$ ssh    # 使用SHH加密协议远程登录
$ man ssh
$ ssh -p 22 username@ip # 远程指定端口登陆
$
$ telnet # 使用TELNET协议远程登录
```

#### 3.6.2 文件下载和拷贝
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
#### 3.6.3 网络配置
```bash
$ ifconfig # 查看、配置、重启或者禁用网络接口
$ ip 
```

#### 3.6.4 网络诊断
```bash
$ ping   # 测试主机之间网络的联通性
$ mtr    # mytraceroute，它集成了 ping、 traceroute、 nslookup 的功能，诊断网络问题非常方便
$ mtr -n www.baidu.com
```
#### 3.6.5 网络连接
```bash 
$ netstat    # 查看网络状态
$ netstat -ntpl # 查询TCP类型端口
$ # -n/--numeric 显示数字地址
$ # -t/--tcp     显示tcp, transmission control protocol 传输控制协议
$ # -p/--program 显示占用的程序
$ # -l/--listen  正在监听
$ netstat -nupl # 查询UDP 端口类型
$ ss         # 查看网络状态
```

[TCP/UDP/PORT](https://zhuanlan.zhihu.com/p/57987304)

#### 3.6.6 流量统计
```bash
$ man ifstat
$ man sar
$ 
```

### 3.7 文件处理
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

### 3.8 文件权限
```bash
$ chmod 755 file1            # change mode of file
$ chomod -R zhusitao folder  # recursively chmod folder to zhusitao
$ chown zhusitao ath.csv     # 将文件的所属者改成zhusitao
$ chown user:group file1     # 将文件的所属者改成user, 将文件的所属组改成group
$ chown :gruop file1         # 只改变所属组, 不改变所有者
$ chgrp -v group_id file1    # 更改文件用户组
$ 
```


### 3.9 进程管理

```bash
$ ps         # 查看进程 show snapshoot of processes
$ top        # 查看实时进程 show real time processes
$ kii pid    # 删除进程 kill process with id pid
$ pkill name # 使用程序名称删除进程 kill process with name 
```

## 4. linux三剑客
### 4.1 awk
文本处理的工具之一,[awk](https://wangchujiang.com/linux-command/c/awk.html)
```bash
$ man awk
$ # -f/--file 从文件中读取awk脚本执行,是命令行执行wak的一种补充方式，可以使用多个-f
$ # -F/--field-separator fs，指定输入文件分隔符号
$ # -v/--assign 指定变量和其对应的值，

$ awk执行方式
$ awk [options] 'script' file(s)
$ awk [options] -f scriptfile files(s)

$ awk 模式和操作
$ 模式主要包括: /正则表达式/: 使用通配符;
$            关系表达式: 使用运算符进行操作;
$            模式匹配表达式: ~匹配 和 !~不匹配;
$            BEGIN语句、pattern语句、END语句;
$ 操作主要包括：有一个或者多个命令、函数、表达式组成，用换行符和分号隔开，并位于大括号内
$             包括。变量数组赋值，输出命令，内置函数和控制语句

$ awk 'BEGIN{ print "start" } pattern{ commands } END{ print "end" }' file # 三个语句块都是可选的

$ awk 的内置变量
$ ARGC 命令行参数数目
$ ARGV 包含命令行参数的数组
$ FILENAME 当前文件名
$ FNR 当前文件的行号，记录数
$ FS 字符分割符号，默认是空格
$ OFS 输出字段分隔符号，默认是一个空格
$ RS 记录分割符, 默认是一个换行符
$ ORS 输出记录分隔符，默认是一个换行符
$ NF 字段数
$ NR 文件行号，记录数

$ awk 传入外部变量
$ awk -v variable=100 '{print $variable}' files

$ awk 运算与判断

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

$ 正则运算符
$ ~, !~
$ ^ 行首
$ $ 行尾
$ . 除了换行符以外的任意单个字符
$ * 前导字符的零个或者多个
$ .*所有字符
$ []字符组内的任一字符
$ [^]不匹配字符组内的任一字符
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

$ AWK 读取
$ # next 语句
$ awk 'NR%4==3{next}{print NR,$0}' fastqs # 去读fastq文件，判断去除其第三行
$ next 语句，循环逐行匹配，匹配就跳过，然后进行下一行匹配，next一般用于多行合并
$ cat text.txt
>chr1
AAAAAAAAA
>chr2
TTTTTTTTT
>chr3
CCCCCCCCC
>chr4
GGGGGGGGG
$ awk '/^>/{ID=$0;next;}{print ID","$0}' text.txt
>chr1,AAAAAAAAA
>chr2,TTTTTTTTT
>chr3,CCCCCCCCC
>chr4,GGGGGGGGG

$ # getline 从标准输入、管道和正在处理的文件之外的其他输入文件或得输入
$ awk 'BEGIN{ "date" | getline out; print out }' test
$ getline 当其左右没有重定向符|或<时：getline作用于当前文件
$ getline 当其左右有重定向符|或<时：getline作用于输入文件，

$ awk 'BEGIN{ "date" | getline out; split(out,mon); print mon[2] }'
$ seq 5 | awk 'BEGIN{getline; print "The first line "$0};{print $0}'

$ 流程控制语句
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


$ 数组应用
$ # awk 数组应用十分便利。数组不必提前声明
$ # 数组下标1-based
$ # Array[1] = 'AAAA'
$ # split 分割函数 and length 函数
$ awk 'BEGIN{info="it is a test";lens=split(info,tA," ");print length(tA),lens;}'

$ awk 内置函数

$ 1.1 awk 算术函数
$ # atan2, cos, sin, exp, log, sqrt, int, rand, strand

$ 1.2 awk 字符串函数
$ # gsub, sub, index, length, blength, substr, match, split, tolower, toupper, sprintf
$ awk 'BEGIN{info="this is a test2010test!";gsub(/[0-9]+/,"!",info);print info}' # gsub 替换
$ awk 'BEGIN{info="this is a test2010test!";print index(info,"test")?"ok":"no found";}' # index 查找字符串
$ awk 'BEGIN{info="this is a test2010test!";print match(info,/[0-9]+/)?"ok":"no found";}' # match 正则表达式匹配查找
$ awk 'BEGIN{info="this is a test2010test!";print substr(info,4,10);}'  #截取字符串， 从第四个字符串开始，，截取10个长度字符串
$ awk 'BEGIN{info="this is a test";split(info,tA," ");print length(tA);for(k in tA){print k,tA[k];}}' # 字符串分割

$ 1.3 一般函数
$ # close, system, getline,
```
### 4.2 grep
grep(global search regular expression and print out the line),全面搜索正则表达式并把行打印。
$ man grep 

#### 4.2.1 pattern syntax 
```bash
$ -E/--extended-regexp # 拓展的正则匹配
$ -F/--fixed-strings   # 将模式固定，不支持正则表达
$ -G/--basic-regexp    # 一般的正则表达
$ -P/--perl-regexp     # perl正则
```

#### 4.2.2 match control 

```bash
$ -e/--regexp # 指定字符串作为模式
$ -f/--file   # 从文件中一行行读取作为模式
$ -i//--ignore-case # 忽略大小写
$ --no-ignore-case  # 不忽略大小写
$ -v/--invert-match # 反向匹配
$ -w/--word-regexp  # 匹配单词
$ -x/--line-regexp  # 匹配整行
```

#### 4.2.3 general outputy control
```bash
$ -c/--count # 打印出匹配的行数， 而不是输出行
$ --color    # 显示颜色
$ -L/--files-without-match # 输出没有匹配的文件名
$ -l/--files-with-matches  # 输出匹配的文件名
$ -m/--max-count     # 指定最大的匹配次数
$ -o/--only-matching # 只打印出匹配的,同时显示文件名
```

#### 4.2.4 output line prefix control
```bash
$ -b/--byte-offset   # 输出匹配的位置，0-base开始；可以和-o配合使用
$ -H/--with-filename # 打印出每个匹配的文件名
$ -h/--no-filename   # 不打印文件名
$ --lable            # 转化输入文件成新的文件名，
$ cat text.txt | grep --label=text 'aaa' -H -i
$ -n/--line-number   # 打印行号，1-based
$ -T/--initial-tab   # 在匹配的文件前使用一个tab符号, 结合-H, -n, -b效果更佳
```

#### 4.2.5 context line control 

```bash
$ -A/--after-context  # 打印匹配行后面NUM行
$ -B/--before-context # 打印匹配文件前面NUM行
$ -C/--context        # 打印前后NUM行
```

#### 4.2.6 grep 正则表达式
```bash
$ ^ # 锚定行的开始 
$ $ # 锚定行的结束
$ . # 匹配一个非换行字符
$ * # 匹配零个或者多个先前字符
$ .* # 任意字符
$ [] # 指定范围内的字符 '[Gg]rep'匹配Grep和grep
$ [^] # 匹配一个不在指定范围内的字符，'[^A-Z]rep'匹配不包含A-Z开头的行
$ \<  # 锚定单词开始
$ \>  # 锚定单词结束
$ x\{m\} # 重复字符x, 'A\{9\}' 匹配9个A
$ x\{m,\}\ # 重复字符x, 至少m次
$ x\{m, n\}} # 重复字符x, 至少m次，至多n次
$ \w # 匹配文字或数字字符 [a-zA-Z0-9]
```

### 4.3 sed
sed (stream editor for filtering and transforming text)流式编辑器
```bash
$ # 命令格式
$ sed [options] 'command' file(s)
$ sed [options] -f scriptfile file(s)
$ # 每个command 由最多两个地址addresses和一个动作action组成,
$ # 每个address可以是正则表达式或者行数，动作见4.3.2
```
#### 4.3.1 sed 参数
```bash
$ -e/--expression # 传入脚本到命令行执行
$ -f/--file       # 从文件传入命令执行
$ -i/--in-place   # 原地编辑文件
$ -n/--quiet/--silent # 仅仅显示脚本处理后的结果
$ -E/--regexp-extended # 正则表达式
$ -V/--version # 显示版本信息
$ -h/--help # 显示帮助信息
```

#### 4.3.2 动作说明

a: 新增 append 
c: 取代 change
d: 删除 delete
i: 插入 insert
p: 打印 print
s: 取代 substitute
y: 转换 transform
q: 退出 quit

#### 4.3.3 print 打印命令
```bash
$ sed -n '1,5p' # 打印第一行到第五行
```

#### 4.3.4 deletion 删除命令
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
$ sed '/^test/d' $ # 删除文件中所有开头是test的行
$ echo this is a test line | sed 's/\w\+/[&]/g' #  \w+表示匹配每一个单词，使用[&]替换它，& 对应之前所匹配到的单词

```

#### 4.3.5 substitute 替换标记
```bash
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

$ sed '2s/AAAA/TTTT/' file # 只替换第二行
$ sed '/DDDD/s/AAAA/TTTT/' file # 替换匹配DDDD的行

$ # 正则
$ echo aaa BBB | sed 's/\([a-z]\+\) \([A-Z]\+\)/\2 \1/'
$ echo this is digit 7 in a number | sed 's/digit \([0-9]\)/\1/'

```

#### 4.3.6 transform 转换
```bash
$ echo ATCG | sed 'y/ATCG/TAGC/' | rev # DNA反向互补
```

#### 4.3.7 quit 退出 
```bash
$ sed '100q' file # 打印前100行，然后退出
```
[sed lecture](https://cs.nyu.edu/~mohri/unix08/lect5.pdf)


## 5 Vim
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
Linux库文件是可执行的公共代码，是为了减少程序开发的重复劳动，分为静态库的动态库。
静态库文件后缀为.a，程序编译时，将静态文件中的代码复制，拷贝到程序生成的可执行文件中。
静态库的优点就是执行程序时不需要外部的函数库的支持，缺点是如果静态函数库变了，程序必须重新编译。
动态库文件后缀为.so，动态库在编译是没有被编译进目标代码，而是程序执行到相关代码时才调用库中的函数。


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


$ esearch -db sra -query PRJNA830912 | efetch -format summary > summary.xml # 根据项目编号下载
$ # 下载GEO数据库
$ esearch -db gds -query GSE201349 | efetch > GSE201349.txt

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
