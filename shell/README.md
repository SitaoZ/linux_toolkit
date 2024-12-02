## Shell Language 

### 1. 变量

#### 1.1 变量定义
在 Shell 编程中，变量是用于存储数据值的名称。
定义变量时，变量名不加`$`符号，变量名和等号之间不能有空格。
- 只包含字母、数字和下划线： 变量名可以包含字母（大小写敏感）、数字和下划线 _，不能包含其他特殊字符。  
- 不能以数字开头： 变量名不能以数字开头，但可以包含数字。  
- 避免使用 Shell 关键字： 不要使用Shell的关键字（例如 if、then、else、fi、for、while 等）作为变量名，以免引起混淆。  
- 使用大写字母表示常量： 习惯上，常量的变量名通常使用大写字母，例如 PI=3.14。  
- 避免使用特殊符号： 尽量避免在变量名中使用特殊符号，因为它们可能与 Shell 的语法产生冲突。  
- 避免使用空格： 变量名中不应该包含空格，因为空格通常用于分隔命令和参数。

```bash
$ # 有效变量
$ ath_gff="Araport11.GFF" # 定义
$ LD_LIBRARY_PATH="/bin/" # 添加动态库常用
$ _var=123                #

$ # 无效变量
if="ABC"  # if是关键词
var_$=123 # 有特殊符号
?var=123  # 有特殊符号
variable with space="value" # 避免空格
```

- 使用语句赋值变量
```bash
$ for sample in `ls *gff`; do echo $sample;done 
```

- 变量类型
```bash
$ # 字符串类型，使用`'`单引号和`"`双引号来定义。
$ my_string='Human genome'
$ my_string="Human genome"

$ # 整型
$ declare -i my_integer=12

$ # 数组变量
$ my_array=(1 2 3 4 5)

$ # 关联数组
$ declare -A associative_array
$ associative_array["name"]="Jone"
$ associative_array["age"]=30

$ # 环境变量
$ echo $PATH

$ # 特殊变量
$ $0 脚本名称
$ $1, $2 表示脚本参数
$ $# 表示传递给脚本参数的数量
$ $? 表示上一个命令的退出状态
```

-  只读变量
使用readonly 可以将变量定义为只读变量，值不会被改变
```bash
$ my_variable="www.baidu.com"
$ readonly my_variable
$ my_variable="www.google.com" # 运行脚本会报错
```

- 删除变量
使用unset命令可以删除变量，变量删除后不能再次使用，unset不能删除只读变量。
```bash
$ unset variable_name
```

#### 1.2 变量的使用

变量使用`$`符号调用
```bash
$ ath="Araport11.gff"
$ echo $ath   # 直接调用
$ echo ${ath} # 加上{}是为了帮助解释器识别变量的边界
```

### 2.字符串
- 单引号的任何字符都会原样输出，单引号字符串中的变量是无效的
- 单引号中不能出现单独一个单引号(转义也不行)

- 双引号里可以有变量
- 双引号可以出现转义字符
```bash
$ gff="Araport11.gff"
$ str="We use the "$gff" as reference \n"
$ echo -e $str
```
- 字符串长度
```bash
$ string="abcd"
$ echo ${#string}
$ echo ${#string[0]}
```

- 提取子字符串
shell 字符串索引值为0
```bash
$ string="The human genome was finished in 2003"
$ echo ${string:1:4}
```

- 查找子字符串
查找字符i或o的位置，谁先出现就计算哪个
```bash
$ string="Human is searching outspace"
$ echo `expr index "$string" io` # 输出 7
```

### 2. 条件语句
```bash
$ # 目录不存在则创建，存在则跳过
$ if [ ! -d CD_KO ];then
$    mkdir CD_KO
$    else
$    echo CD_KO exist
$ fi

$ # 文件不存在
$ if [ ! -f "/data/filename" ];then
$  echo "文件不存在"
$  else
$  touch /data/filename
$ fi
```

- shell中文件比较符号
```bash
-e 判断对象是否存在
-d 判断对象是否存在，并且为目录
-f 判断对象是否存在，并且为常规文件
-L 判断对象是否存在，并且为符号链接
-h 判断对象是否存在，并且为软链接
-s 判断对象是否存在，并且长度不为0
-r 判断对象是否存在，并且可读
-w 判断对象是否存在，并且可写
-x 判断对象是否存在，并且可执行
-O 判断对象是否存在，并且属于当前用户
-G 判断对象是否存在，并且属于当前用户组
-nt 判断file1是否比file2新  [ "/data/file1" -nt "/data/file2" ]
-ot 判断file1是否比file2旧  [ "/data/file1" -ot "/data/file2" ]
```

