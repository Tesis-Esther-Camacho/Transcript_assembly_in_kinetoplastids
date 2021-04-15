#!/usr/bin/env python
#-*- coding: UTF-8 -*-

from __future__ import division
from collections import defaultdict
import sys
import re
import os
import operator

def diccionarioAnidado():
	return defaultdict(diccionarioAnidado)

syntax = """----------------------------------------------------------------------------------------------------
Usage:	polycistronic_transcripts.py file.gff transcriptome_mod_6.gtf
----------------------------------------------------------------------------------------------------"""
##########################################################################################
if len(sys.argv) != 3:
	print syntax
	sys.exit()
##########################################################################################	
#--obtenemos la anotación del fichero gff del tritrypdb
gff = open (sys.argv[1], 'r')

#--genes[chr][strand][gene_id] = [x, y, att]
genes = diccionarioAnidado()
for line in gff:
	line = line.rstrip('\n')
	col = line.split('\t')[:-1]
	
	if line[:7] == '##FASTA': #--final del fichero (para no leer el fasta asociado)
		break #--Finish
	elif line.startswith('#'): #--Saltamos la cabecera
		continue
	elif col[2] == 'chromosome': #--Saltamos la cabecera
		continue
	elif col[2] == 'gene':
		att = [x.split('=')[1] for x in line.split('\t')[-1].split(';')]
		genes[col[0]][col[6]][att[0]] = [col[3], col[4]]
	else:
		continue

gff.close()
##########################################################################################	
#--obtenemos la información de los transcritos ensamblados (gtf)
gtf = open (sys.argv[2], 'r')

#--trans[chr][strand][trans_id] = [x, y, att]
trans = diccionarioAnidado()

for line in gtf:
	line = line.rstrip('\n')
	col = line.split('\t')[:-1]
	att = line.split('\t')[-1].split('"')
	trans[col[0]][col[6]][att[1]] = [col[3], col[4], att]
	
gtf.close()
##########################################################################################
#--Analizamos cada transcrito para determinar si existen 2 o más genes anotados en sus coordenadas
out = open(sys.argv[2].split('.')[0][:-1] + '7.gtf', 'w')
corregidos = []
for chr in trans.keys():
	for strand in trans[chr].keys():
		for tid in trans[chr][strand].keys():

			tx =  int(trans[chr][strand][tid][0])
			ty =  int(trans[chr][strand][tid][1])

			gids = []

			for gen in genes[chr][strand].keys():
				#--Genes contenidos enteramente dentro de las coordenadas del transcrito
				if tx <= int(genes[chr][strand][gen][0]) <= ty and tx <= int(genes[chr][strand][gen][1]) <= ty:
					gids.append(gen)

			if len(gids) > 1: #--más de un gen anotado en el transcrito						
				corregidos.append(tid)
				transcript = [chr, 'CBMSO', 'transcript', str(tx), str(ty), '.', strand, '.']
			
				new_att = trans[chr][strand][tid][2]
				new_att[13] = '|'.join(gids)
				new_att[15] = 'polycistronic'
				
				out.write('\t'.join(transcript) + '\t' + '"'.join(new_att) + '\n')
				#--Elimino el transcrito del diccionario
				del trans[chr][strand][tid]


					
##########################################################################################
#--Escribimos los mensajeros que no han sido modificados
for chr in trans.keys():
	for strand in trans[chr].keys():
		for tid in trans[chr][strand].keys():

			tx = trans[chr][strand][tid][0]
			ty = trans[chr][strand][tid][1]

			att = trans[chr][strand][tid][2]
			transcript = [chr, 'CBMSO', 'transcript', tx, ty, '.', strand, '.']
			out.write('\t'.join(transcript) + '\t' + '"'.join(att) + '\n')				

out.close()
##########################################################################################
print 'Mensajeros policistrónicos:', len(corregidos)
out_ids = open ('transcripts_mod_6_ids.txt', 'w')
for i in corregidos:
	out_ids.write(i +'\n')
out_ids.close()
##########################################################################################


