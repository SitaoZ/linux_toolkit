#  使用sort和uniq快速求文件的交集、并集、差集
# 主要使用的uniq 
# uniq -d (only print duplicate lines, one for each group)
# uniq -u (only print unique lines)
# 并集
sort a.txt b.txt | uniq | wc

# 交集
sort a.txt b.txt | uniq -d | wc

# 差集 1
# a.txt - b.txt 
sort a.txt b.txt b.txt | uniq -u | wc

# b.txt - a.txt 
sort b.txt a.txt a.txt | uniq -u | wc

