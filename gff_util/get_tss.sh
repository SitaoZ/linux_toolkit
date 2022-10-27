# tss
cat xxx.gff3 | awk '$3=="gene" {print $0}' | grep protein_coding | awk -v OFS="\t" '{if ($7=="+") {print $1, $4, $4+1}} else {print $1, $5-1, $5}' > tss.bed

# promoter 5k upstream from tss
cat xxx.gff3 | awk '$3=="gene" {print $0}' | grep protein_coding | awk -v OFS="\t" '{if ($7=="+") {print $1, $4, $4+5000} else {print $1, $5-5000, $5}}' > promoters.bed
