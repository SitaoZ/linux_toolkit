echo 'AGTCATGCAGTGCNNNNT' | rev | tr 'ACTG' 'TGAC'

echo 'AGTCATGCAGTGCNNNNT' | python -c "import sys;from Bio.Seq import Seq;a = [print(Seq(i.strip()).reverse_complement()) for i in sys.stdin];"
