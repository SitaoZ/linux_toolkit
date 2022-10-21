for chr in {1..22}; do awk -v chr=$chr -v threshold=10 '$1 == chr && $6 > threshold' xxx.vcf > xxx_chr$chr.vcf; done
# https://samtools.github.io/hts-specs/VCFv4.2.pdf
# shell中的变量在awk中无法识别，可以使用-v指定变量传入awk
# -v chr=$chr, 将shell中的变量$chr传入awk


awk '$1 == 5 && $7 == "PASS"' xxx.vcf > xxx_filter.vcf
# 挑选五号染色体上的SNP
# 过滤值为PASS
