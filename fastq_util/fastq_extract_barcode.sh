paste <(zcat sample_1.fq.gz|paste - - - -) <(zcat sample_2.fq.gz|paste - - - - ) | awk -v FS="\t" -v OFS="\n" 'FNR==NR {samples[$2]=$1; next} {barcode = substr($6,0,6); if(samples[barcode]) { print $1,$2,$3,$4>>samples[barcode]"_1.fq"; print $5,$6,$7,$8>>samples[barcode]"_2.fq"}}' samples.txt -

# 1. paste - - - -, 四行打包成一行, 默认tab分隔符
# 2. awk -v FS="\t" 指定输入符号为\t, OFS="\n" 指定输出分割符为"\n"
# 3. FNR==NR 处理第一个文件; awk 记录数为文件的行号(line number)
# 当前文件行号和总行号相等，表示第一个文件.
# 当前记录数: FNR(Number of input Record in current input File)
# 总记录数: NR(total Number of Records seen so far )
# 4. substr, awk 内置函数,截取字符串,substr(string,position,length)
# 5. samples.txt - 标准输入samples.txt 内容到awk,(即处理的第一个文件)
