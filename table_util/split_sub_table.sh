cat <(less total.csv | head -n 1) <(less total.csv| grep AT3G49430) > result.csv
# < 表示重定向输入符号,usage: command < file 
# () 表示整体执行括弧中的命令
# head -n 1 表示表头
# <() 接 grep 表示另外一个筛选的内容
# cat 合并; > 重定向文件
