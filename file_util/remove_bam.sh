# wildcard
rm *bam

# find and xargs
find . -name "*.bam" | xargs rm
