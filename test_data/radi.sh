
TIMES=1

READSIM=/home/dino/Documents/proj/readsim-1.6/src/readsim.py
READS=450
RATE=2

#stvara readove 






for i in `seq $TIMES`
do
  /home/dino/Documents/proj/amos/bin/toAmos -s reads.$READS.10x.fasta -o reads.$READS.afg
  /home/dino/Documents/proj/amos/bin/minimus reads.$READS

  ## Dump overlaps
  /home/dino/Documents/proj/amos/bin/bank-report -b reads.$READS.bnk/ OVL > overlaps.$READS.afg


  ## Dump unitig layous
  /home/dino/Documents/proj/amos/bin/bank-report -b reads.$READS.bnk/ LAY > layouts.$READS.afg

  READS=$((RATE*READS))  
done

