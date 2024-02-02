## Table of content
* [Linux bash strict model](# Linux-bash-strict-model)
* [Linux symbol](#Linux-symbol)
## Linux-bash-strict-model 
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

## Linux-symbol
linux中一些符号具有特定的含义，需要注意
```bash
$ # # 井号 (comments)

$ # ~ 家目录 (home)

$ # . 当前目录
$ # .. 上一级目录

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
# # `command` 等价于 $(command) 
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

$ cd - # 表示回到上一次所在的目录

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
