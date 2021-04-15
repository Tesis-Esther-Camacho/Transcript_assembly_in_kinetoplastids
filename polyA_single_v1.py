#!/usr/bin/env python
#-*- coding: UTF-8 -*-

from __future__ import division
import sys
import re
import subprocess
from collections import defaultdict
import os

syntax = """
------------------------------------------------------------------------------------------------------------
 Usage: python polyA_trim.pl.py <reads.fastq> <min_polyA_length> <genome.fasta> <genome_index(bt2)> <bow>
------------------------------------------------------------------------------------------------------------
"""
""" Descripción: Buscamos secuencias de polyT* en el extremo 5' de las secuencias, cuya longitud debe ser al menos la definida por
 el usuario (min_polyA_length). Recorremos el fichero fastq buscando este patrón, cuando se localiza una lectura con éste, se
 escribe en un fichero fastq temporal que se alinea contra el genoma. Entonces tomamos las coordenadas del alineamiento y 
 analizamos si el polyA/polyT (depende de la orientación del alineamiento**) eliminado estaba o no en el genoma, de modo que si 
 esta secuencia estaba en el genoma no recortamos la lecturas, y si no estaba si la recortamos. 
 
 
*All reads are in the forward orientation (5'--3')
Usual cases
5-------------------------------------------------------------------AAAAAAAA3 (fragment from a transcript)
						                3---------------------------TTTTTTTT5 (read)

5----------------------------3 												  (read)
3-------------------------------------------------------------------TTTTTTTT5 (fragment from a transcript)


The unusual cases: len(fragment) - len(poly) < len(read), i.e. very short fragments (we do not consider these cases)
5-----------------------------------------AAAAAAAA3 (fragment from a transcript)
			  3---------------------------TTTTTTTT5 (read)

5-----------------------------------------AAA3      (read)
3-----------------------------------------TTTTTTTT5 (fragment from a transcript)
 

**
5-------------------AAAAAAAA---------------------------------------------TTTTTTTT-----------------------------------------3 (genomic polyA/polyT)
3-------------------TTTTTTTT5 (false polyA)          3-------------------AAAAAAAA5 (not a polyA)



               	   5AAAAAAAA-------------------3 (not a polyA)          5TTTTTTTT-------------------3 (false polyA)
3-------------------TTTTTTTT---------------------------------------------AAAAAAAA-----------------------------------------5 (genomic polyA/polyT)


------------------------------------------------------------------------------------------------------------------------------
"""

def refasta(infile):
	""" Reformatea un fichero fasta para que la secuencia esté en una sóla línea """
	file = open (infile, 'r')
	s = ''
	seq = {}
	for line in file:
		line = line.rstrip('\n')
		if line:
			if line[0] == '>':
				if s:
					seq[name] = s
					s = ''
				name = line
			else:
				s = s + line
	#--La última
	seq[name] = s
	file.close()
	
	#--Re-escribimos
	file = open (infile, 'w')
	for s in seq.keys():
		file.write(s + '\n')
		file.write(seq[s] + '\n')
	file.close()
	
#--Parameters control
if len(sys.argv) != 6:
	print syntax
	sys.exit()	

bin_path = os.path.dirname(os.path.realpath(__file__))

#--Parameters definition 
min_polyA_length = int(sys.argv[2])
polya = re.compile('^T{%d,}'%min_polyA_length, re.IGNORECASE) #--Regex definition
genome = sys.argv[3]
gindex = sys.argv[4]
bow = sys.argv[5]
pol = defaultdict(list)

#--Output_files
polya_file = open ('polyA_temp.fastq', 'w')

#--Open fastq file and polyA search
fastq_infile = open (sys.argv[1], 'r')

count = 0
while True:
	name = fastq_infile.readline().rstrip('\n')
	seq = fastq_infile.readline().rstrip('\n')
	coment = fastq_infile.readline().rstrip('\n')
	qual = fastq_infile.readline().rstrip('\n')

	if not name: #--Finish
		break
	
	id = name[1:].split(' ')[0]
	
	#--Buscamos que la lectura no haya sido recortada por poseer sequencia leader
	if id[-1] == 'L':
		continue
	
	match = polya.search(seq)
	
	if match: #--potential leader sequence
		x = match.start() #--index of the first matched base
		y = match.end() #--index of the next base in the sequence
		
		if y < 0.75*(len(seq)): #--if the polyA is longer than 75% of the read there is to less sequence to get a proper alignment, so do not trim these reads
			count = count + 1
			#--Escribimos estos potenciales hits en un fichero temporal y en un diccionario
			record = [name, seq[y:], coment, qual[y:]]
			polya_file.write('\n'.join(record) + '\n')
			pol[id] = [name, seq, coment, qual, y, 0] #--Añadimos al diccionario la posición de corte y una etiqueta que luego nos permitirá descriminar SL genuinos de los falsos positivos

print 'Potentical PolyA hits: ' + str(count)
	
fastq_infile.close()
polya_file.close()


#--Alineamos los potenciales polyA detectados
subprocess.call(bow + ' -x ' + gindex + ' -U polyA_temp.fastq -S temp.sam', shell=True)

#--Recorremos el fichero sam para extraer las secuencias genómicas inmediatamente anteriores al inicio del alineamiento de los potenciales polyA
#--De este modo podemos analizar si la secuencia estaba en el genoma y por lo tanto es un falso positivo
gtf_temp = open('temp.gtf', 'w')
sam = open ('temp.sam', 'r')
chr = {}

for line in sam:
	line = line.rstrip('\n')
	col = line.rstrip('\n').split('\t')
	
	#--skip header lines and take chromosome sizes
	
	if '@' in line[0]: 
		if col[0] == '@SQ':
			chr[col[1].split(':')[1]] = int(col[2].split(':')[1])
			continue
		else:
			continue

	#--The read do not align
	if col[2] == '*':		
		#del sl[col[0]]
		continue
		
	#--The read align
	else: 
		#--Obtenermos la información del bitwise flag
		bitwise = int(col[1])
		bits = bin(bitwise)[2:]
		byte = "{0:0>12}".format(bits)

		#--Obtenemos las coordenadas en función de la orientación del alineamiento
		if int(byte[-5]): #--int(byte[-5])==1 (reverse)
			strand = '-'					
			cigar = col[5]
			let = re.split('\d+', cigar)[1:] #--Dividimos por número	
			val = re.split('\D+', cigar)[:-1] #--Dividimos por letras

			coord = 0
			for l, v in zip(let, val):
				if l == 'M': #--Fragmento de lectura alineado
					coord = coord + int(v)
				elif l == 'D': #--Delección en la lectura
					coord = coord + int(v)
				else: #--Las inserciones no las sumamos
					continue
					
			coord = int(col[3]) + coord - 1
			
			x = coord + 1
			y = coord + pol[col[0]][-2]
			
			if y > chr[col[2]]:
				y = chr[col[2]]
						
			if (y - x) + 1 < min_polyA_length: #--Si la secuencia esta fuera de los límites de cromosoma damos por bueno al polyA, porque no podemos discriminarlo
				pol[col[0]][-1] = 1
			
			else:		
				#--Escribimos las coordenadas del fragmento eliminado en un fichero gtf para extraer la secuencia genómica
				transcript = [col[2], 'CBMSO', 'transcript', str(x), str(y), '.', '-', '.']
				att = ['gene_id "', col[0], '"; ', 'transcript_id "',col[0], '"; '] 	
				gtf_temp.write('\t'.join(transcript) + '\t' + ''.join(att) + '\n')
	
		else: #--int(byte[-5])==0 (forward)
			strand = '+'
			coord = int(col[3])				
			
			y = coord - 1
			x = coord - pol[col[0]][-2]
			
			if x <= 0:
				x = 1
			
			if (y - x) + 1 < min_polyA_length:  #--Si la secuencia esta fuera de los límites de cromosoma damos por bueno al SL, porque no podemos discriminarlo
				pol[col[0]][-1] = 1
				
			else:
				#--Escribimos las coordenadas del fragmento eliminido en un fichero gtf para extraer la secuencia genómica
				transcript = [col[2], 'CBMSO', 'transcript', str(x), str(y), '.', '+', '.']
				att = ['gene_id "', col[0], '"; ', 'transcript_id "',col[0], '"; '] 	
				gtf_temp.write('\t'.join(transcript) + '\t' + ''.join(att) + '\n')
			
sam.close()
gtf_temp.close()

#--Extraemos las secuencias genómicas
subprocess.call (bin_path + '/cufflinks/gffread -w temp.fa -g ' + genome + ' temp.gtf', shell = True)

#--Reformateo del fichero temp.fa por si la secuencia ocupa de más de 1 línea (en el caso de que la secuencia eliminada sea mayor de 60 pb)
refasta('temp.fa')

#--Analizamos las secuencias genómicas en busca de secuencias leader genómicas (mismas condiciones que en las lecturas)
seqs =  open ('temp.fa', 'r')
while True:

	name = seqs.readline().rstrip('\n')[1:]
 	seq = seqs.readline().rstrip('\n')
 	
 	if not name: #--Finish
 		break
 	
 	#--Como extraigo las secuencias en la misma orientación que la lectura alineada, en todos los casos debo buscar polyT
 	#--Como podría darse el caso de que las coordenadas del fragmento eliminado estuvieran fuera de los cromosomas
 	#--tomo la longitud de la secuencia extraída como longitud para comprobar si es sólo un polyT
 	
 	if seq != 'T'*len(seq):
 		pol[name][-1] = 1

seqs.close()	

#--Volvemos a recorrer el fichero con las lecturas originales y ahora si recortamos los SL que sabemos son correctos
fastq_infile = open (sys.argv[1], 'r')
polya_outfile = open('polyA_ids', 'w')
reads_outfile = open(gindex + '_SL_polyA_trimmed.fastq', 'w')

count = 0

while True:
	name = fastq_infile.readline().rstrip('\n')
	seq = fastq_infile.readline().rstrip('\n')
	coment = fastq_infile.readline().rstrip('\n')
	qual = fastq_infile.readline().rstrip('\n')
	
	if not name: #--Finish
		break
	
	if pol[name[1:].split(' ')[0]] and pol[name[1:].split(' ')[0]][-1]:
		count = count + 1
		new_name = name.split(' ')[0] + 'P ' + ' '.join(name.split(' ')[1:])
		#new_name = name + 'L' #--Añado una L para identificar las lecturas después
		record = [new_name, seq[pol[name[1:].split(' ')[0]][-2]:], coment, qual[pol[name[1:].split(' ')[0]][-2]:]] 
		polya_outfile.write(name[1:].split(' ')[0] + '\n') #--Guardar sólo los ids y luego sacarlos del sam general
		reads_outfile.write('\n'.join(record) + '\n')
	
	else:
		record = [name, seq, coment, qual] 
		reads_outfile.write('\n'.join(record) + '\n')	

print 'Valid polyA hits: ' + str(count)

polya_outfile.close()
reads_outfile.close()
fastq_infile.close()
	
os.remove('polyA_temp.fastq')
os.remove('temp.sam')
os.remove('temp.gtf')
os.remove('temp.fa')
