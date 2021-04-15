#!/usr/bin/env python
#-*- coding: UTF-8 -*-

from __future__ import division
from collections import defaultdict
import sys
import re
import operator

""" Programa que recorta los extremos 3' de los mensajeros en función de las coordenadas
de secuencias con polyA encontradas """

if __name__ == '__main__':

	if len(sys.argv) !=3:
		print 'Syntax: python polyA_recut_v1.py polyA_unused.gtf transcripts_polyA_both.gtf (reviewed)'
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
	trans[chr][id] = [int(col[3]), int(col[4]), col[8]]

tfile.close()

##########################################################################################

###--polyA_unused.gtf--#######################################################################	

#--Obtenemos la información de los polya del fichero SAM.
pfile = open (sys.argv[1], 'r')
#--starts = {chr:{strand:{coord:reads}}}
ends = defaultdict(lambda : defaultdict(dict))

for line in pfile:
	line = line.rstrip('\n')
	col = line.split('\t')
	att = col[8].split('"')[1::2]

	if col[6] == '+':
		ends[col[0]][col[6]][int(col[3])] = int(att[2])
	else:
		ends[col[0]][col[6]][int(col[4])] = int(att[2])

pfile.close()

##########################################################################################

###--Asignación de polyA--################################################################	
 	
gtf = open('transcripts_polyA_both_cutted.gtf', 'w')

for chr in sorted(trans.keys()):
	
	for id in sorted(trans[chr].keys()):		
			
		tlen = trans[chr][id][1] - trans[chr][id][0]
		margen = tlen * 0.2 #--Admitimos buscar en el 20% final del mensajero
		
		#--temp = {coord:reads}
		temp = {} #--vaciamos los SL temporales
		
		orientation = trans[chr][id][2].split('"')[-2]
	
		if orientation == 'F':	#--Busqueda de polyA de la misma hebra que estén en el 3' del mensajero	
			strand = '+' #--hebra del mensajero
			for coord in sorted(ends[chr]['-'].keys()):
				if (trans[chr][id][1] - margen) < coord < (trans[chr][id][1] + 50): #--polyA en el extremo 3' del mensajero
					temp[coord] = ends[chr]['-'][coord] #--temp = {coord:reads}	
					del ends[chr]['-'][coord] #--eliminar coord del diccionario para que no se repita
					
		elif orientation == 'R':	
			strand = '-' #--hebra del mensajero
			for coord in sorted(ends[chr]['+'].keys()):
				if (trans[chr][id][0] - 50) < coord < (trans[chr][id][0] + margen):
					temp[coord] = ends[chr]['+'][coord] #--temp = {coord:reads}	
					del ends[chr]['+'][coord] #--eliminar coord del diccionario para que no se repita

		elif orientation == 'N': #--No debe ser dividido
			transcript = [chr, 'CBMSO', 'transcript', str(trans[chr][id][0]), str(trans[chr][id][1]), '.', '.', '.']
			att = '"'.join(trans[chr][id][2].split('"')[:-3]) + '"; '
			gtf.write('\t'.join(transcript) + '\t' + att + '\n')
			del trans[chr][id]
			continue
			
		else:
			print id
			print 'ERROR: transcript not reviewd'
			sys.exit()

		
		#--control de corte
		if not temp: 
			print id
			print 'ERROR: all transcripts in this step have polyA'
			sys.exit()
		
		#--transcritos con 1 sólo polyA
		elif len(temp) == 1: 
			if strand == '+': #--hebra del mensajero
				ppal = sorted(temp)[0]
				transcript = [chr, 'CBMSO', 'transcript', str(trans[chr][id][0]), str(ppal), '.', '+', '.']		
				att1 = '"'.join(trans[chr][id][2].split('"')[:-7]) + '"; '
				att1 = re.sub('X', 'F', att1)
				att = ['polyA "', str(temp[ppal]), '"; ', 'alt_polyA "', 'ND', '";']
				gtf.write('\t'.join(transcript) + '\t' + att1 + ''.join(att) + '\n')
				
			if strand == '-': #--hebra del mensajero
				ppal = sorted(temp)[0]
				transcript = [chr, 'CBMSO', 'transcript', str(ppal), str(trans[chr][id][1]), '.', '-', '.']
				att1 = '"'.join(trans[chr][id][2].split('"')[:-7]) + '"; '
				att1 = re.sub('X', 'R', att1)
				att = ['polyA "', str(temp[ppal]), '"; ', 'alt_polyA "', 'ND', '";']
				gtf.write('\t'.join(transcript) + '\t' + att1 + ''.join(att) + '\n')
		
		#--Más de un polyA 	
		else: 
			
			#--ppales={coord_ppal:{'reads': X, alt: [coord(reads), coord(reads)]}}		
			ppales = defaultdict(dict) #--Almacena el polya principal y sus alternativos	
			
			#--temp_2 = [(coord, reads), ..] ordenadas por lecturas de más a menos
			temp_2 = sorted(temp.iteritems(), key=operator.itemgetter(1), reverse = True)
			
			#--tomamos como principal el que más se exprese
			ppal = temp_2[0][0]
			
			#--control de elección del ppal, por si hay algún otro con el mismo número de lecturas.
			for t in temp_2:
				if t[1] == temp[ppal]:
					if strand == '+' and t[0] > ppal: #--El más 3'
							ppal = t[0]
					elif strand == '-'  and t[0] < ppal: #--El más 5'
							ppal = t[0]
					else:
						continue		

			#--asigno el ppal en un diccionario
			ppales[ppal]['reads'] = temp[ppal]
			ppales[ppal]['alt'] = []
			del temp[ppal] #--elimino el polyA ppal del diccionario temp

			#--anoto el resto como alternativos
			for polya in temp.keys():
				ppales[ppal]['alt'].append(str(polya) + '(' + str(temp[polya]) + ')')
				del temp[polya]
			
			#--Imprimimos las nuevas coordenadas del mensajero
			if strand == '+':
				transcript = [chr, 'CBMSO', 'transcript', str(trans[chr][id][0]), str(ppal), '.', '+', '.']
				att1 = '"'.join(trans[chr][id][2].split('"')[:-7]) + '"; '
				att1 = re.sub('X', 'F', att1)
				att = ['polyA "', str(ppales[ppal]['reads']), '"; ', 'alt_polyA "', '|'.join(ppales[ppal]['alt']) , '";']
				gtf.write('\t'.join(transcript) + '\t' + att1 + ''.join(att) + '\n')
				
			if strand == '-':
				transcript = [chr, 'CBMSO', 'transcript', str(ppal), str(trans[chr][id][1]), '.', '-', '.']
				att1 = '"'.join(trans[chr][id][2].split('"')[:-7]) + '"; '
				att1 = re.sub('X', 'R', att1)
				att = ['polyA "', str(ppales[ppal]['reads']), '"; ', 'alt_polyA "', '|'.join(ppales[ppal]['alt']), '";']
				gtf.write('\t'.join(transcript) + '\t' + att1 + ''.join(att) + '\n')			
			
		#--elimino el mensajero que ha sido modificado
		del trans[chr][id]

##########################################################################################

###--Imprimimos el resto de polyA no utilizados--########################################

 	poly = open ('polyA_unused_mod_1.gtf', 'w')
 	
 	a = 1
 	
 	for chr in sorted(ends.keys()):
 		for strand in sorted(ends[chr].keys()):
 			if strand == '-':
 				for coord in sorted(ends[chr][strand].keys()):
					transcript = [chr, 'CBMSO', 'transcript', str(coord-30), str(coord), '.', '-', '.']
					att = ['gene_id "', 'polyA_' + str(a), '"; ', 'transcript_id "', 'polyA_' + str(a), '"; ', 'reads "', str(ends[chr][strand][coord]), '";'] 	
					poly.write('\t'.join(transcript) + '\t' + ''.join(att) + '\n')
					a = a + 1
 			else:
 				for coord in sorted(ends[chr][strand].keys(), reverse=True):
					transcript = [chr, 'CBMSO', 'transcript', str(coord), str(coord+30), '.', '+', '.']
					att = ['gene_id "', 'polyA_' + str(a), '"; ', 'transcript_id "', 'polyA_' + str(a), '"; ', 'reads "', str(ends[chr][strand][coord]), '";'] 	
					poly.write('\t'.join(transcript) + '\t' + ''.join(att) + '\n')
					a = a + 1
 						
 	poly.close()
 		
##########################################################################################
