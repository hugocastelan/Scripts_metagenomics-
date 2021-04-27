#!/bin/bash

# -*- coding: utf-8 -*- # NAME [prodigal_function.sh]	Version [1.0]
# AUTHOR  Hugo Castelan Sanchez 
# CREATED (2019)
# USAGE prodigal_function.sh directory with assembly data 

# DESCRIPTION
# Performs structural annotation using prodigal from multiple fasta files


for file in *.fa ;
  do prodigal -i $file -o mygenes$file -a my_protein_$file -p meta ; 
done
