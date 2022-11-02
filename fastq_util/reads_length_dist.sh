# compute reads length distribution from a fastq file
zcat xx.fq.gz | awk 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}' > reads.dist

# read species
zcat xx.fq.gz | awk 'NR%4 == 2 {species[$0]++} END {for (s in species) {print s, species[s]}}' > reads.count
