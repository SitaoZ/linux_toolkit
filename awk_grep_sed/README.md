## Table of content
* [awk](#1.awk)
* [grep](#2.grep)
* [sed](#3.sed)
### 1.awk
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

#### 1.1 awk 内置变量
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
$ awk -v variable=100 '{print variable}' files
$ for i in `seq 10`;do awk -v a="$i" '$1==a' xxx.csv >> result.csv;done
```

```bash
$ # 当处理两个文件时，可以使用FNR 和 NR来进行挑选
$ # 如在第一个文件中，选出第二个文件出现的第一个文件中的内容
$ cat file1.txt
$ # gene1  gene1_function_descriptions
$ # gene2  gene2_function_descriptions
$ # gene3  gene3_function_descriptions
$ cat file2.txt
$ # gene1
$ # gene3

$ awk 'FNR==NR{a[$1]=$0} NR>FNR{print a[$1]}' file1.txt file2.txt > result
$ cat result
$ # gene1  gene1_function_descriptions
$ # gene3  gene3_function_descriptions
```


#### 1.2 awk 运算与判断
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
#### 1.3 awk 正则运算
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
#### 1.4 awk 读取
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


$ awk multi fasta file to one_line fasta
$ cat multi.fa | awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' > one_line.fa
$ awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' multi.fa > one_line.fa
```

#### 1.5 awk 流程控制
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

#### 1.6 awk 数组应用
```bash
$ # awk 数组应用十分便利。数组不必提前声明
$ # 数组下标1-based
$ # Array[1] = 'AAAA'
$ # split 分割函数 and length 函数
$ awk 'BEGIN{info="it is a test";lens=split(info,tA," ");print length(tA),lens;}'

$ # mRNA.bed
$ cat Araport11_GFF3_genes_transposons.Mar92021.gff | awk '$3=="mRNA"' | awk -v OFS="\t" '{split($9,a,";"); print $1,$4,$5,a[1],$6,$7}' | sed 's/ID=//g' > mRNA.bed
```
#### 1.7 awk 内置函数
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

- awk转义引号
```bash
$ awk '{print "\""}'   # 双引号
$ awk '{print "'\''"}' # 单引号
```




### 2.grep
grep(global search regular expression and print out the line),全面搜索正则表达式并把行打印。
```bash
$ man grep 
```
#### 2.1 pattern syntax 
```bash
$ -E/--extended-regexp # 拓展的正则匹配
$ -F/--fixed-strings   # 将模式固定，不支持正则表达
$ -G/--basic-regexp    # 一般的正则表达
$ -P/--perl-regexp     # perl正则
```

#### 2.2 match control 

```bash
$ -e/--regexp # 指定字符串作为模式
$ -f/--file   # 从文件中一行行读取作为模式
$ -i//--ignore-case # 忽略大小写
$ --no-ignore-case  # 不忽略大小写
$ -v/--invert-match # 反向匹配
$ -w/--word-regexp  # 匹配单词
$ -x/--line-regexp  # 匹配整行
```

#### 2.3 general output control
```bash
$ -c/--count # 打印出匹配的行数， 而不是输出行
$ --color    # 显示颜色
$ -L/--files-without-match # 输出没有匹配的文件名
$ -l/--files-with-matches  # 输出匹配的文件名
$ -m/--max-count     # 指定最大的匹配次数
$ -o/--only-matching # 只打印出匹配的,同时显示文件名
```

#### 2.4 output line prefix control
```bash
$ -b/--byte-offset   # 输出匹配的位置，0-base开始；可以和-o配合使用
$ -H/--with-filename # 打印出每个匹配的文件名
$ -h/--no-filename   # 不打印文件名
$ --lable            # 转化输入文件成新的文件名，
$ cat text.txt | grep --label=text 'aaa' -H -i
$ -n/--line-number   # 打印行号，1-based
$ -T/--initial-tab   # 在匹配的文件前使用一个tab符号, 结合-H, -n, -b效果更佳
```

#### 2.5 context line control 

```bash
$ -A/--after-context  # 打印匹配行后面NUM行
$ -B/--before-context # 打印匹配文件前面NUM行
$ -C/--context        # 打印前后NUM行
```

#### 2.6 grep 正则表达式
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

#### 2.7 zgrep 

```bash
$ zgrep -v "^#" GCF_000001405.40_GRCh38.p14_genomic.gff.gz | awk 'BEGIN{FS="\t";OFS"\t"}$2=="RefSeqFE"&&$3!="biological_region"'
```

### 3.sed
sed (stream editor for filtering and transforming text)流式编辑器
```bash
$ # 命令格式
$ sed [options] 'command' file(s)
$ sed [options] -f scriptfile file(s)
$ # 命令必须由单引号包住；或者使用双引号包住，双引号主要在传入变量时使用
$ # 每个command 由最多两个地址(addresses)和一个动作(action)组成,
$ # 每个地址可以是正则表达式或者行数，动作见4.3.2
```
#### 3.1 sed 参数

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

#### 3.2 动作说明

a: 新增 append  
c: 取代 change  
d: 删除 delete  
i: 插入 insert  
p: 打印 print, 通常和-n一起使用，打印特定的行
s: 取代 substitute  
y: 转换 transform  
q: 退出 quit  


#### 3.3 append 添加新行

```bash
$ sed '/ATCG/a atcg' file.txt # 在匹配的行后面添加新行
```

#### 3.4 deletion 删除命令
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

#### 3.5 print 打印命令
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

#### 3.6 substitute 替换标记
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

#### 3.7 transform 转换
```bash
$ echo ATCG | sed 'y/ATCG/TAGC/' | rev # DNA反向互补
```

#### 3.8 quit 退出 
```bash
$ sed '100q' file # 打印前100行，然后退出
```
#### 3.9 组合多个表达式
```bash
$ sed '表达式1; 表达式2' # sed '表达式1' | sed '表达式2'
$ sed 's/chr/Chr/g; s/geneid/GeneID/g' xxx.fa # 同时多个替换操作
```

#### 3.10 引用

```bash
$ # 传入shell变量时，需要使用双引号
$ gene_id=AT1G79550
$ cat file.txt | sed "s/$gene_id/PGK/g" # 替换基因名
```

#### 3.11 多点编辑
```bash
$ nl /etc/passwd | sed -e '3,$d' -e 's/bash/blueshell/' # -e表示多点编辑，第一个编辑命令删除第三行到末尾的数据，第二条命令搜索bash替换为blueshell
```
[sed lecture](https://cs.nyu.edu/~mohri/unix08/lect5.pdf)
