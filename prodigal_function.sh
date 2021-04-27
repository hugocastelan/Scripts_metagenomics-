



for file in *.fa ;
  do prodigal -i $file -o mygenes$file -a my_protein_$file -p meta ; 
done
