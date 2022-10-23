awk '/^>/{s=++d".fa"}{print > s}' xxx.fa
# /^>/ 匹配fasta 文件名
# s=++d".fa" 生成文件名，一每一个fasta格式自追加编号
# {print > s} 输出到该文件名

csplit -z -q -n 4 -f sequence_ sequences.fasta /\>/ {*}  
# /[正则表达式]/   #匹配文本样式，比如/SERVER/，从第一行到包含SERVER的匹配行。
# {*}     #表示根据匹配重复执行分割，直到文件尾停止，使用{整数}的形式指定分割执行的次数。
# -z或--elide-empty-files 删除长度为0 Byte文件。
# -q或-s或--quiet或--silent 不显示指令执行过程。
# -n      #指定分割后的文件名后缀的数字个数。比如01、02、03等。
# -f      #指定分割后的文件名前缀。

