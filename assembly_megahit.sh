#!/bin/bash

# -*- coding: utf-8 -*- # NAME [assembly_megahit.sh]	Version [1.0]
# AUTHOR  Hugo Castelan Sanchez 
# CREATED (2019)
# USAGE assembly_megahit.sh directory with fastq

# DESCRIPTION
# Assembly multiple fastq 


for  f1 in *_1.fq 
 do           
 f2=${fi%%_1.fq}"_2.fq"
 /home/sonia/megahit/megahit  -1  $f1  -2 $f2 -o "assembly$fi.out" 
 done

