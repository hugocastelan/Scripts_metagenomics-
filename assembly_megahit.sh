#!/bin/bash
for  f1 in *_1.fq 
 do           
 f2=${fi%%_1.fq}"_2.fq"
 /home/sonia/megahit/megahit  -1  $f1  -2 $f2 -o "assembly$fi.out" 
 done

