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
