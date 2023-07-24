# long read QC
zcat xaa.fq.gz | seqtk seq -A -L 10000 - | grep -v "^>" | tr -dc "ACGTNacgtn" | wc -m


# zcat ( concatenates the compressed fastq files into one stream )
# seqtk ( converts to fasta format and drops reads less than 10k )
# grep ( -v excludes lines starting with “>”, i.e. fasta headers )
# tr ( -dc removes any characters not in set “ACGTNacgtn” )
# wc ( -m counts characters )
