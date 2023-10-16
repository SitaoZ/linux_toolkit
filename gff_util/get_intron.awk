BEGIN{OFS="\t"}
{ start[NR]=$4; 
  end[NR]=$5; 
  strand[NR]=$7; 
  ID[NR]=$11; 
  chr[NR]=$1; 
  ens[NR]=$9; 
  symb[NR]=$10} 

END { 
    for (i=1; i<=NR; i++){
        if(ID[i]==ID[i+1]){ 
            if(strand[i]=="+"){
                intron_start=start[i+1]-1;
                intron_end=end[i]+1;
                print chr[i], 
                      intron_end, 
                      intron_start,
                      "intron",
                      strand[i], 
                      ens[i],
                      symb[i], 
                      ID[i]

            }   
            else{
                intron_start=end[i]+1;
                intron_end=start[i+1]-1;
                print chr[i], 
                      intron_start, 
                      intron_end,
                      "intron",
                      strand[i], 
                      ens[i], 
                      symb[i], ID[i]
            }   
        }   
    }   
}
