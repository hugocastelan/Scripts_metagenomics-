# -*- coding: utf-8 -*- # NAME [gc_contect.py]	Version [1.0]
# AUTHOR  Hugo Castelan Sanchez 
# CREATED (2019-08)
# USAGE gc_contect.py -f fasta or multifasta from assembly.fa -o output.txt 
# DESCRIPTION
#Calcular el contenido de GC y AT de un Genoma 
#Calcula el total del genoma 
import argparse
import re

#Para parsear el archivo 
parser=argparse.ArgumentParser()
requiredNamed=parser.add_argument_group('required named arguments')
requiredNamed.add_argument('-f', '--file', help='Input fasta or multifasta from assembly')

args=parser.parse_args()
file=open(args.file)

#Iniciado las variables en cero
g=0
c=0
a=0
t=0
seq=""

lineas=file.readlines()#Lee todas la lineas del archivo
#print(lineas)

#Saltar los encabezados 
for line in lineas:
	if line.startswith('>'): #Este es el encabezado
	    encabezado=line
	else:
		seq=line
		seq=seq.rstrip("\n")#Limpia el salto de linea
		seq=seq.lower()#Hace todas las letras a minusculas 
		for caracter in seq:
			if caracter == "g":
				g+=1
			if caracter == "c":
				c+=1
			if caracter == "a":
				a+=1
			if caracter == "t":
				t+=1
#La longitud de toda la secuencia es la suma de todos los caracteres
longitud=(g+c+a+t)
print("length of sequence is= ",longitud,"pb")
total_de_gc=((g+c)/longitud)*100
total_de_at=((a+t)/longitud)*100
print("GC contect is=  " ,str(round(total_de_gc, 2)),"%")
print("AT contect is= " ,str(round(total_de_at,2)),"%")	

file.close()
