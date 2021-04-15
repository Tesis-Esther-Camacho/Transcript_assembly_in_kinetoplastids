#!/usr/bin/env python
#-*- coding: UTF-8 -*-

from __future__ import division
from collections import defaultdict
import sys
import re
import operator

""" Programa que divide mensajeros ensamblados en función de los SLs (secuencias de mini-
exon) que se encuentren dentro de las coordenadas del mensajero, siendo los SL clasifica-
dos en función del número de lecturas y de la proximidad unos de otros """


if len(sys.argv) != 3:
	print 'Syntax: python SL_cut.py SL_unused.gtf transcripts_SL_both.gtf (reviewed)'
	sys.exit()


###--transcripts.gtf--####################################################################	

#--Obtenemos las coordenadas de los mensajeros ensamblados con cufflinks

tfile = open (sys.argv[2], 'r')

#--trans = {chr:{id:[x, y]}
trans = defaultdict(dict)

for line in tfile:
	if line == '\n':
		break
	col = line.rstrip('\n').split('\t')
	att = col[8].split('"')[1::2]
	id = att[0]
	chr = col[0]
	trans[chr][id] = [int(col[3]), int(col[4]), att[4]]

tfile.close()

##########################################################################################

# for c in trans.keys():
# 	for t in trans[c].keys():
# 		print c, t, trans[c][t]


###--SL_unused.gtf--#######################################################################	

#--Obtenemos la información de los SL del fichero SAM.
SL_file = open (sys.argv[1], 'r')
#--starts = {chr:{strand:{coord:reads}}}
starts = defaultdict(lambda : defaultdict(dict))

for line in SL_file:
	line = line.rstrip('\n')
	col = line.split('\t')
	att = col[8].split('"')[1::2]

	if col[6] == '+':
		starts[col[0]][col[6]][int(col[3])] = int(att[2])
	else:
		starts[col[0]][col[6]][int(col[4])] = int(att[2])

SL_file.close()

#print starts
##########################################################################################

###--División por SLs--###################################################################
	
gtf = open ('transcripts_SL_both_cutted.gtf', 'w')	

#--Para cada cromosomas
for chr in sorted(trans.keys()):

	#--Para cada transcrito
	for id in sorted(trans[chr].keys()): #--tomamos cada transcrito	
			
		i = 1 #--Para nombrar los mensajeros que vamos dividiendo
		
		sl = {} #--sl[coord] = reads
		
		orientation = trans[chr][id][2]
		
		if orientation == 'F':	
			strand = '+'
			#--Buscamos SL de cadena positiva
			for coord in sorted(starts[chr][strand].keys()):
				if (trans[chr][id][0] - 50) < coord < (trans[chr][id][1] + 50): #--SL dentro de las coordenadas del transcrito (+-50 pb)
					sl[coord] = starts[chr][strand][coord]
					del starts[chr][strand][coord] 
					
		elif orientation == 'R':
			strand = '-'
			#--Buscamos SL de cadena negativa
			for coord in sorted(starts[chr][strand].keys()):
				if (trans[chr][id][0] - 50) < coord < (trans[chr][id][1] + 50): #--SL dentro de las coordenadas del transcrito (+-50 pb)
					sl[coord] = starts[chr][strand][coord]
					del starts[chr][strand][coord]
					
		elif orientation == 'N': #--No debe ser cortado
			transcript = [chr, 'CBMSO', 'transcript', str(trans[chr][id][0]), str(trans[chr][id][1]), '.', '.', '.']
			att = ['gene_id "', id + '.X', '"; ', 'transcript_id "', id + '.X', '"; ', 'SL "', 'ND', '"; ', 'Alt_SL "', 'ND', '"; ']		
			gtf.write('\t'.join(transcript) + '\t' + ''.join(att) + '\n')
			continue
			
		else:
			print id
			print 'ERROR: not reviewed transcript'
			sys.exit()		


		#--Clasificamos los inicios
		""" Primero buscamos el que más se expresa, luego buscamos en el entorno de éste (+-500 pb), y los SLs
		en dicho entorno son considerados como alternativos. Después buscamos el siguiente que más se expresa, 
		habiendo quitado el anterior y sus alternativos, y repetimos hasta que no queden más SLs dentro de las 
		coordenadas del mensajero """
		
		#--ppales={coord_ppal:{'reads': X, alt: [coord(reads), coord(reads)]}}		
		ppales = defaultdict(dict) #--Almacena los inicios principales asociados a los alternativos

		#print id, sl
		while sl:

			sl_2 = sorted(sl.iteritems(), key=operator.itemgetter(1), reverse = True) #--sl_2 = [(coord, reads), ..] ordenadas por lecturas de más a menos
			ppal = sl_2[0][0] #--tomamos la coordenada del que más se expresa como principal

			#--control de elección del ppal, por si hay algún otro con el mismo número de lecturas.
			for s in sl_2: #--s = (coord, reads)
				if s[1] == sl[ppal]:
					if strand == '+' and s[0] < ppal: #--El más 5'
							ppal = s[0]
					elif strand == '-'  and s[0] > ppal: #--El más 3'
							ppal = s[0]
					else:
						continue		

			#--asigno el ppal en un diccionario
			ppales[ppal]['reads'] = sl[ppal]
			ppales[ppal]['alt'] = []
			del sl[ppal] #--elimino el sl ppal del diccionario sl

			#--busco en el entorno del principal y los defino como alternativos
			for s in sl.keys():
				if ppal - 500 < s < ppal + 500:		
					ppales[ppal]['alt'].append(str(s) + '(' + str(sl[s]) + ')')
					del sl[s]
		
		#--Divimos en función de la hebra
		if strand == '+':
			if sorted(ppales.keys())[0]-trans[chr][id][0] >= 500:
				transcript = [chr, 'CBMSO', 'transcript', str(trans[chr][id][0]), str(sorted(ppales.keys())[0]-1), '.', '.', '.']
				att = ['gene_id "', id + '.' + str(i) + '.X', '"; ', 'transcript_id "', id + '.' + str(i) + '.X', '"; ', 'SL "', 'ND', '"; ', 'Alt_SL "', 'ND', '"; ']		
				gtf.write('\t'.join(transcript) + '\t' + ''.join(att) + '\n')
				i = i + 1
				
			for z in range(0, len(ppales)-1):
				transcript = [chr, 'CBMSO', 'transcript', str(sorted(ppales.keys())[z]), str(sorted(ppales.keys())[z+1]-1), '.', '+', '.']
				att = ['gene_id "', id + '.' + str(i) + '.F', '"; ', 'transcript_id "', id + '.' + str(i) + '.F', '"; ', 'SL "', str(ppales[sorted(ppales.keys())[z]]['reads']), '"; ', 'Alt_SL "', 'ND', '"; ']
				#--Anotamos los alternativos
				if ppales[sorted(ppales.keys())[z]]['alt']:
					att[10] = str('|'.join(ppales[sorted(ppales.keys())[z]]['alt']))	
				else:
					att[10] = 'ND'
											
				gtf.write('\t'.join(transcript) + '\t' + ''.join(att) + '\n')
				i = i + 1
			
			#--Imprimimos el último (conservando la información de polyA si tenía asignado
			transcript = [chr, 'CBMSO', 'transcript', str(sorted(ppales.keys())[-1]), str(trans[chr][id][1]), '.', '+', '.']
			att = ['gene_id "', id + '.' + str(i) + '.F', '"; ', 'transcript_id "', id + '.' + str(i) + '.F', '"; ', 'SL "', str(ppales[sorted(ppales.keys())[-1]]['reads']), '"; ', 'Alt_SL "', 'ND', '"; ']
			#--Anotamos los alternativos
			if ppales[sorted(ppales.keys())[-1]]['alt']:
				att[10] = str('|'.join(ppales[sorted(ppales.keys())[-1]]['alt']))
			else:
				att[10] = 'ND'
			
			gtf.write('\t'.join(transcript) + '\t' + ''.join(att) + '\n')
			i = i + 1
			
			#--elimino el transcrito del diccionario para que sólo queden los no modificados
			del trans[chr][id]
									
		else: #--Para la cadena negativa (el lado oscuro)
			if trans[chr][id][1] - sorted(ppales.keys(), reverse=True)[0] >= 500:
				transcript = [chr, 'CBMSO', 'transcript', str(sorted(ppales.keys(), reverse = True)[0]+1), str(trans[chr][id][1]), '.', '.', '.']
				att = ['gene_id "', id + '.' + str(i) + '.X', '"; ', 'transcript_id "', id + '.' + str(i) + '.X', '"; ', 'SL "', 'ND', '"; ', 'Alt_SL "', 'ND', '"; ']
				gtf.write('\t'.join(transcript) + '\t' + ''.join(att) + '\n')
				i = i + 1
				
			for z in range(0, len(ppales)-1):
				transcript = [chr, 'CBMSO', 'transcript', str(sorted(ppales.keys(), reverse = True)[z+1]+1), str(sorted(ppales.keys(), reverse = True)[z]), '.', '-', '.']
				att = ['gene_id "', id + '.' + str(i) + '.R', '"; ', 'transcript_id "', id + '.' + str(i) + '.R', '"; ', 'SL "', str(ppales[sorted(ppales.keys(), reverse = True)[z]]['reads']), '"; ', 'Alt_SL "', 'ND', '"; ']
			
				#--Anotamos los alternativos
				if ppales[sorted(ppales.keys(), reverse = True)[z]]['alt']:
					att[10] = str('|'.join(ppales[sorted(ppales.keys(), reverse = True)[z]]['alt']))
				else:
					att[10] = 'ND'
						
				gtf.write('\t'.join(transcript) + '\t' + ''.join(att) + '\n')
				i = i + 1
			
			#--Imprimimos el último
			transcript = [chr, 'CBMSO', 'transcript', str(trans[chr][id][0]), str(sorted(ppales.keys(), reverse = True)[-1]), '.', '-', '.']
			att = ['gene_id "', id + '.' + str(i) + '.R', '"; ', 'transcript_id "', id + '.' + str(i) + '.R', '"; ', 'SL "', str(ppales[sorted(ppales.keys(), reverse = True)[-1]]['reads']), '"; ', 'Alt_SL "', 'ND', '"; ']
			
			#--Anotamos los alternativos
			if ppales[sorted(ppales.keys(), reverse = True)[-1]]['alt']:
				att[10] = str('|'.join(ppales[sorted(ppales.keys(), reverse = True)[-1]]['alt']))
				
			else:
				att[10] = 'ND'
			
			gtf.write('\t'.join(transcript) + '\t' + ''.join(att) + '\n')
			i = i + 1
	
			#--elimino el transcrito del diccionario para que sólo queden los no modificados
			del trans[chr][id]


##########################################################################################

###--Imprimimos los Sls no utilizados--###################################################	

out = open ('SL_unused_mod_1.gtf', 'w')
a = 1

for chr in sorted(starts.keys()):
	for strand in sorted(starts[chr].keys()):
		if strand == '+':
			for coord in sorted(starts[chr][strand].keys()):
					transcript = [chr, 'CBMSO', 'transcript', str(coord), str(coord+30), '.', '+', '.']
					att = ['gene_id "', 'SL_' + str(a), '"; ', 'transcript_id "', 'SL_' + str(a), '"; ', 'reads "', str(starts[chr][strand][coord]), '";'] 	
					out.write('\t'.join(transcript) + '\t' + ''.join(att) + '\n')
					del starts[chr][strand][coord]
					a = a + 1
		else:
			for coord in sorted(starts[chr][strand].keys(), reverse=True):
					transcript = [chr, 'CBMSO', 'transcript', str(coord-30), str(coord), '.', '-', '.']
					att = ['gene_id "', 'SL_' + str(a), '"; ', 'transcript_id "', 'SL_' + str(a), '"; ', 'reads "', str(starts[chr][strand][coord]), '";'] 	
					out.write('\t'.join(transcript) + '\t' + ''.join(att) + '\n')
					del starts[chr][strand][coord]
					a = a + 1
out.close()

##########################################################################################


