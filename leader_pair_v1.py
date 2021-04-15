#!/usr/bin/env python
#-*- coding: UTF-8 -*-

from __future__ import division
import sys
import re
import subprocess
from collections import defaultdict
import os


syntax = """
-------------------------------------------------------------------------------------------------------------------------------------------------------------------
Usage: python leader_trim.pl.py <leader.fa> <reads_1.fastq> <reads_2.fastq> <min_leader_length> 
		<max_mm> <genome.fasta> <genome_index(bt2)> <bow>

Parameters description:

<leader.fa>			Fasta file containing the sequence to be searched (leader)
<unaligned.fastq>			Fastq file to search
<min_leader_length>	Initial portion to be searched (minimum length of the leader)
<max_mm>			Maximum number of missmatches in the whole leader
<genome.fasta>			Genome sequence
<genome_index(bt2)> 	Bowtie2 genome index 
-------------------------------------------------------------------------------------------------------------------------------------------------------------------
"""

"""Descripción: 
Buscamos secuencias leader en el 5' de las secuencias*. Para ello localizamos primero una
semilla en las lecturas, cuyo tamaño viene definido por el usuario. Si se encuentra esta
semilla evaluamos el resto de la secuencia, contando los posibles errores de secuencia,
tras lo que evaluamos si supera el número máximo de errores definido por el usuario,
teniendo en cuenta que el fragmento localizado puede ser inferior a la longitud de la
secuencia leader, es decir, calculamos el porcentaje máximo de errores admitidos en
función de la longitud total de la secuencia leader y usamos eso como valor umbral para
compararlo con el porcentaje de errores encontrados en el fragmento de secuencia leader
localizada en la lectura.

*Theoretically, as leader sequences (kinetoplastos) are located at 5' end of the
transcripts, it is not possible to find the leader sequence in the reverse-complement
orientation, unless the fragment from which the read came from was smaller than the read
length plus leader sequence length, i.e., a very small fragment of the beginning of a
fragment that was read from the end. 

Usual cases
5FFFFFFFFF------------------------------------------------------------------3 (fragment from a transcript)
						                3-----------------------------------5 (read)

5FFFFFFFFF------------------------3 (read)
3RRRRRRRRR------------------------------------------------------------------5 (fragment from a transcript)


The unusual cases: len(fragment) - len(leader) < len(read), i.e. very short fragments (we do not consider these cases)
5FFFFFFFFF-------------------------------3 (fragment)
	 3RRRR-------------------------------5 (read)

5FFFFFFFFF------------------------3 (read)
3RRRRRRRRR--------------------------------5 (fragment)
"""

def align(query = [], subject = []):

	""" función que compara dos secuencias suministradas como listas, calculando el número
	de missmatches como porcentaje de la longitud del query"""

	mm = 0
	for q, s in zip(query, subject):
		if q != s:
			mm = mm + 1 

	return mm/len(query)

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
if len(sys.argv) != 9:
	print syntax
	sys.exit()	

bin_path = os.path.dirname(os.path.realpath(__file__))

#--leader: take the leader sequence from a fasta file (only the first sequence)
leader_file = open (sys.argv[1], 'r')
leader_name = leader_file.readline()[1:].rstrip('\n')
leader = leader_file.readline().rstrip('\n').upper()
leader_file.close()

#--Parameters definition 
min_leader_length = int(sys.argv[4])
max_mm = float(sys.argv[5])
max_error = max_mm / len(leader) #--Maximum error rate for full leader match
genome = sys.argv[6]
gindex = sys.argv[7]
bow = sys.argv[8]

#--Open fastq file and leader search
fastq_infile1 = open (sys.argv[2], 'r')
fastq_infile2 = open (sys.argv[3], 'r')
leader_outfile1 = open('seqleader_temp_1.fastq', 'w')
leader_outfile2 = open('seqleader_temp_2.fastq', 'w')

count = 0
sl = defaultdict(list)

while True:
	name1 = fastq_infile1.readline().rstrip('\n')
	seq1 = fastq_infile1.readline().rstrip('\n')
	coment1 = fastq_infile1.readline().rstrip('\n')
	qual1 = fastq_infile1.readline().rstrip('\n')

	name2 = fastq_infile2.readline().rstrip('\n')
	seq2 = fastq_infile2.readline().rstrip('\n')
	coment2 = fastq_infile2.readline().rstrip('\n')
	qual2 = fastq_infile2.readline().rstrip('\n')	
	
	rl = len(seq1) #--Lo necesito más adelante
	
	if not name1: #--Finish
		break
	
	trim = 0 #--Buscamos primero en el par 1, y si no encontramos, buscamos en el 2
	for i in range (0, len(leader) - min_leader_length + 1):

		query = leader[i:]
		subject = seq1[:len(leader) - i]

		#--Si el ratio de errores que devuelve la función align no supera el máximo lo consideramos como un potencial hit
		if align(list(query), list(subject)) <= max_error:
			count = count + 1
			
			#--Escribimos estos potenciales hits en un fichero temporal y en un diccionario para procesarlo más adelante
			record1 = [name1, seq1[len(query):], coment1, qual1[len(query):]]
			record2 = [name2, seq2, coment2, qual2] 
			leader_outfile1.write('\n'.join(record1) + '\n')
			leader_outfile2.write('\n'.join(record2) + '\n')	
			sl[name1.split(' ')[0][1:]] = [name1, seq1, coment1, qual1, 1, len(query), 0] #--Añadimos al diccionario el par que es, la posición de corte y una etiqueta que luego nos permitirá descriminar SL genuinos de los falsos positivos
			trim = 1
			break
			
	if not trim: #--Si no encontramos la secuencia leader en la lectura 1 la buscamos en la lectura 2
	
		for i in range (0, len(leader) - min_leader_length + 1):

			query = leader[i:]
			subject = seq2[:len(leader) - i]

			#--Si el ratio de errores que devuelve la función align no supera el máximo lo consideramos como un potencial hit
			if align(list(query), list(subject)) <= max_error:
				count = count + 1
			
				#--Escribimos estos potenciales hits en un fichero temporal y en un diccionario para procesarlo más adelante
				record1 = [name1, seq1, coment1, qual1] 
				record2 = [name2, seq2[len(query):], coment2, qual2[len(query):]]
				leader_outfile1.write('\n'.join(record1) + '\n')
				leader_outfile2.write('\n'.join(record2) + '\n')	
				sl[name2.split(' ')[0][1:]] = [name2, seq2, coment2, qual2, 2, len(query), 0] #--Añadimos al diccionario el par que es, la posición de corte y una etiqueta que luego nos permitirá descriminar SL genuinos de los falsos positivos
				trim = 1
				break

		
print 'Potential Leader hits: ' + str(count)

leader_outfile1.close()
leader_outfile2.close()
fastq_infile1.close()
fastq_infile2.close()

#--Alineamos los potenciales SL detectados
subprocess.call(bow + ' -x ' + gindex + ' -1 seqleader_temp_1.fastq -2 seqleader_temp_2.fastq -S temp.sam', shell=True)

#--Recorremos el fichero sam para extraer las secuencias genómicas inmediatamente anteriores al inicio del alineamiento de los potenciales SL
#--De este modo podemos analizar si la secuencia leader estaba en el genoma y por lo tanto es un falso positivo
gtf_temp = open('temp.gtf', 'w')
sam = open ('temp.sam', 'r')
chr = {}

#--leemos la cabecera
while True:
	line = sam.readline()
	line = line.rstrip('\n')
	col = line.rstrip('\n').split('\t')	
	
	if col[0] == '@SQ':
		chr[col[1].split(':')[1]] = int(col[2].split(':')[1])
		continue
	
	if col[0] == '@PG': #--Para leer sólo la cabecera
		break

#--Leemos la parte del alineamiento
while True:	
	line1 = sam.readline().rstrip('\n')
	line2 = sam.readline().rstrip('\n')
	
	if not line1: #--Final del fichero
		break
	
	read_name = line1.split('\t')[0]

	#--Analizamos el bitwise flag de la primera lectura que aparece y determinamos cual de las dos es el pair 1 y 2
	col = line1.split('\t')
	bitwise = int(col[1])
	bits = bin(bitwise)[2:]
	byte = "{0:0>12}".format(bits)
		
	if int(byte[-7]) == 1: #--Par 1
		pair1 = line1
		pair2 = line2	
	else:
		pair1 = line2
		pair2 = line1				

	#--Tomamos los datos del diccionario y analizamos la lectura en la que se detectó el SL
	if sl[read_name][-3] == 1: #--SL detectado en el par 1
		col = pair1.split('\t')
		bitwise= int(col[1])
		bits = bin(bitwise)[2:]
		byte = "{0:0>12}".format(bits)
	
	else: #--SL detectado en el par 2
		col = pair2.split('\t')
		bitwise= int(col[1])
		bits = bin(bitwise)[2:]
		byte = "{0:0>12}".format(bits)
	
 	#--The read do not align, not trim
 	if int(byte[-3]) == 1:		
 		del sl[col[0]]
 		continue	
	
	#--The read align
	else:
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
			y = coord + sl[col[0]][-2]
			
			if y > chr[col[2]]:
				y = chr[col[2]]
						
			if (y - x) + 1 < min_leader_length: #--Si la secuencia esta fuera de los límites de cromosoma damos por bueno al SL, porque no podemos discriminarlo
				sl[col[0]][-1] = 1
			
			else:	
				#--Escribimos las coordenadas del fragmento eliminado en un fichero gtf para extraer la secuencia genómica
				transcript = [col[2], 'CBMSO', 'transcript', str(x), str(y), '.', '-', '.']
				att = ['gene_id "', col[0], '"; ', 'transcript_id "',col[0], '"; '] 	
				gtf_temp.write('\t'.join(transcript) + '\t' + ''.join(att) + '\n')
	
		else: #--int(byte[-5])==0 (forward)
			strand = '+'
			coord = int(col[3])				
			
			y = coord - 1
			x = coord - sl[col[0]][-2]
			
			if x <= 0:
				x = 1
			
			if (y - x) + 1 < min_leader_length:  #--Si la secuencia esta fuera de los límites de cromosoma damos por bueno al SL, porque no podemos discriminarlo
				sl[col[0]][-1] = 1
				
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
 	
 	query = leader[-(sl[name][-2]):]
 	
 	if len(seq) < len(query): #--para las secuencias de los extremos en los que falta secuencia genómica, aunque siempre hay más que el mínimo (filtrado arriba)
 		query = leader[-(len(seq)):]
	 	if align(list(query), list(seq)) > max_error:
			sl[name][-1] = 1 #--Cambiamos la etiqueta del diccionaria a 1, es decir, que es un SL genuino y debe ser recortado	
 	else:
 		if align(list(query), list(seq)) > max_error:
 			sl[name][-1] = 1 #--Cambiamos la etiqueta del diccionaria a 1, es decir, que es un SL genuino y debe ser recortado

seqs.close()	

#--Volvemos a recorrer el fichero con las lecturas originales y ahora si recortamos los SL que sabemos son correctos
fastq_infile1 = open (sys.argv[2], 'r')
fastq_infile2 = open (sys.argv[3], 'r')
leader_outfile = open('seqleader_ids', 'w')
reads_outfile1 = open(gindex + '_SL_trimmed_1.fastq', 'w')
reads_outfile2 = open(gindex + '_SL_trimmed_2.fastq', 'w')

count = 0

while True:
	name1 = fastq_infile1.readline().rstrip('\n')
	seq1 = fastq_infile1.readline().rstrip('\n')
	coment1 = fastq_infile1.readline().rstrip('\n')
	qual1 = fastq_infile1.readline().rstrip('\n')

	name2 = fastq_infile2.readline().rstrip('\n')
	seq2 = fastq_infile2.readline().rstrip('\n')
	coment2 = fastq_infile2.readline().rstrip('\n')
	qual2 = fastq_infile2.readline().rstrip('\n')
	
	if not name1: #--Finish
		break
	
	name = name1.split(' ')[0][1:]
	
	if sl[name] and sl[name][-1]:
		
		count = count + 1
		new_name = name1.split(' ')[0] + 'L ' + ' '.join(name1.split(' ')[1:]) #--Añado una L al final del nombre de la lectura para identificarlas después
		
		if sl[name][-3] == 1: #--recortar el par 1		
			record1 = [new_name, seq1[sl[name][-2]:], coment1, qual1[sl[name][-2]:]]
			record2 = [new_name, seq2, coment2, qual2] 
			leader_outfile.write(name + '\t' + str(sl[name][-3]) + '\n') #--Guardar sólo los ids y luego sacarlos del sam general
			reads_outfile1.write('\n'.join(record1) + '\n')
			reads_outfile2.write('\n'.join(record2) + '\n')
			
		else:
			record1 = [new_name, seq1, coment1, qual1]
			record2 = [new_name, seq2[sl[name][-2]:], coment2, qual2[sl[name][-2]:]] 
			leader_outfile.write(name + '\t' + str(sl[name][-3]) + '\n') #--Guardar sólo los ids y luego sacarlos del sam general
			reads_outfile1.write('\n'.join(record1) + '\n')
			reads_outfile2.write('\n'.join(record2) + '\n')
			
	else:
		record1 = [name1, seq1, coment1, qual1]
		record2 = [name2, seq2, coment2, qual2]
		reads_outfile1.write('\n'.join(record1) + '\n')
		reads_outfile2.write('\n'.join(record2) + '\n')

print 'Valid leader hits: ' + str(count)

leader_outfile.close()
reads_outfile1.close()
reads_outfile2.close()
fastq_infile1.close()
fastq_infile2.close()
	
# os.remove('seqleader_temp_1.fastq')
# os.remove('seqleader_temp_2.fastq')
# os.remove('temp.sam')
# os.remove('temp.gtf')
# os.remove('temp.fa')

