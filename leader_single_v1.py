#!/usr/bin/env python
#-*- coding: UTF-8 -*-

from __future__ import division
import sys
import re
import subprocess
from collections import defaultdict
import os


syntax = """
------------------------------------------------------------------------------------------------------------------------------
Usage: python leader_trim.pl.py <leader.fa> <reads.fastq> <min_leader_length> <max_mm> <genome.fasta> <genome_index(bt2)> <bow>

Parameters description:

<leader.fa>			Fasta file containing the sequence to be searched (leader)
<unaligned.fastq>			Fastq file to search
<min_leader_length>	Initial portion to be searched (minimum length of the leader)
<max_mm>			Maximum number of missmatches in the whole leader
<genome.fasta>		
<genome_index(bt2)> 
------------------------------------------------------------------------------------------------------------------------------
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
if len(sys.argv) != 8:
	print syntax
	sys.exit()	

bin_path = os.path.dirname(os.path.realpath(__file__))

#--leader: take the leader sequence from a fasta file (only the first sequence)
leader_file = open (sys.argv[1], 'r')
leader_name = leader_file.readline()[1:].rstrip('\n')
leader = leader_file.readline().rstrip('\n').upper()
leader_file.close()

#--Parameters definition 
min_leader_length = int(sys.argv[3])
max_mm = float(sys.argv[4])
max_error = max_mm / len(leader) #--Maximum error rate for full leader match
genome = sys.argv[5]
gindex = sys.argv[6]
bow = sys.argv[7]


#--Open fastq file and leader search
fastq_infile = open (sys.argv[2], 'r')
leader_outfile = open('seqleader_temp.fastq', 'w')
count = 0
sl = defaultdict(list)

while True:
	name = fastq_infile.readline().rstrip('\n')
	seq = fastq_infile.readline().rstrip('\n')
	coment = fastq_infile.readline().rstrip('\n')
	qual = fastq_infile.readline().rstrip('\n')
	
	rl = len(seq) #--Lo necesito más adelante
	
	if not name: #--Finish
		break
	
	for i in range (0, len(leader) - min_leader_length + 1):

		query = leader[i:]
		subject = seq[:len(leader) - i]

		#--Si el ratio de errores que devuelve la función align no supera el máximo lo consideramos como un potencial hit
		if align(list(query), list(subject)) <= max_error:
			count = count + 1
			
			#--Escribimos estos potenciales hits en un fichero temporal y en un diccionario para procesarlo más adelante
			record = [name, seq[len(query):], coment, qual[len(query):]] 
			leader_outfile.write('\n'.join(record) + '\n')	
			sl[name[1:].split(' ')[0]] = [name, seq, coment, qual, len(query), 0] #--Añadimos al diccionario la posición de corte y una etiqueta que luego nos permitirá descriminar SL genuinos de los falsos positivos
			break
		
print 'Potential Leader hits: ' + str(count)

leader_outfile.close()
fastq_infile.close()

#--Alineamos los potenciales SL detectados
subprocess.call(bow + ' -x ' + gindex + ' -U seqleader_temp.fastq -S temp.sam', shell=True)

#--Recorremos el fichero sam para extraer las secuencias genómicas inmediatamente anteriores al inicio del alineamiento de los potenciales SL
#--De este modo podemos analizar si la secuencia leader estaba en el genoma y por lo tanto es un falso positivo
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
		del sl[col[0]]
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
				#--Escribimos las coordenadas del fragmento eliminado en un fichero gtf para extraer la secuencia genómica
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
fastq_infile = open (sys.argv[2], 'r')
leader_outfile = open('seqleader_ids', 'w')
reads_outfile = open(gindex + '_SL_trimmed.fastq', 'w')

count = 0

while True:
	name = fastq_infile.readline().rstrip('\n')
	seq = fastq_infile.readline().rstrip('\n')
	coment = fastq_infile.readline().rstrip('\n')
	qual = fastq_infile.readline().rstrip('\n')
	
	if not name: #--Finish
		break
	
	if sl[name[1:].split(' ')[0]] and sl[name[1:].split(' ')[0]][-1]:
		count = count + 1
		new_name = name.split(' ')[0] + 'L ' + ' '.join(name.split(' ')[1:])
		#new_name = name + 'L' #--Añado una L para identificar las lecturas después
		record = [new_name, seq[sl[name[1:].split(' ')[0]][-2]:], coment, qual[sl[name[1:].split(' ')[0]][-2]:]] 
		leader_outfile.write(name[1:].split(' ')[0] + '\n') #--Guardar sólo los ids y luego sacarlos del sam general
		reads_outfile.write('\n'.join(record) + '\n')
	
	else:
		record = [name, seq, coment, qual] 
		reads_outfile.write('\n'.join(record) + '\n')	

print 'Valid leader hits: ' + str(count)

leader_outfile.close()
reads_outfile.close()
fastq_infile.close()
	
os.remove('seqleader_temp.fastq')
os.remove('temp.sam')
os.remove('temp.gtf')
os.remove('temp.fa')

