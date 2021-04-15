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
		print 'Syntax: python recort.py file.sam transcripts_SL_cutted.gtf'
		sys.exit()

###--transcripts.gtf--####################################################################	

	#--Obtenemos las coordenadas de los mensajeros ensamblados con cufflinks
	
	tfile = open (sys.argv[2], 'r')
	
	#--trans = {chr:strand:{id:[x, y]}}
	trans = defaultdict(lambda : defaultdict(dict))
	
	for line in tfile:
		if line == '\n':
			break
		col = line.rstrip('\n').split('\t')
		id = col[8].split('"')[1]
		chr = col[0]
		strand = col[6]
		trans[chr][strand][id] = [int(col[3]), int(col[4]), col[8]]
	
	tfile.close()
	
##########################################################################################

###--poliA.sam--##########################################################################	
	
	#--Obtenemos la información de los SL del fichero SAM.
	pfile = open (sys.argv[1], 'r')
	#--ends = {chr:{strand:{coord:reads}}}
	ends = defaultdict(lambda : defaultdict(dict))
	
	for line in pfile:
		line = line.rstrip('\n')
		if '@' in line[0]: #--saltamos la cabecera del sam
			continue
				
		col = line.rstrip('\n').split('\t')
		chr = col[2]
		
		#--extraemos las coordenadas y la orientación
		""" Los secuencias de polyAs que pertenecen a mensajeros en F, alinean en R, y 
		viceversa, de modo que asigno la cadena + a los que alinean en R y la - a los que
		lo hacen en F, es decir, que los que alinean en R pertenecen a mensajeros F y 
		viceversa"""

		#--Obtenermos la información del bitwise flag
		bitwise = int(col[1])
		bits = bin(bitwise)[2:]
		byte = "{0:0>12}".format(bits)
		
		if int(byte[-5]): #--int(byte[-5])==1 (reverse), alineamiento reverse pero pertenece a mensajeros forward
			strand = '+'
			cigar = col[5]
			let = re.split('\d+', cigar)[1:]		
			val = re.split('\D+', cigar)[:-1]

			coord = 0
			for l, v in zip(let, val):
				if l == 'M': #--Fragmento de lectur
					coord = coord + int(v)
				elif l == 'D': #--Delección en la l
					coord = coord + int(v)
				else: #--Las inserciones no las sum
					continue
			coord = int(col[3]) + coord - 1

		else: #--int(byte[-5])==0 (forward), alineamiento forward pero pertenece a mensajeros reverse
			strand = '-'
			coord = int(col[3])
			
		#--Metemos los datos en un diccionario
		try:
			ends[chr][strand][coord] += 1
		except KeyError:
			ends[chr][strand][coord] = 1

	pfile.close()

	#--Imprimimos un fichero con todos los polyA detectados y el número de lecturas
	
	todos = open ('polyA_detected.gtf', 'w')
	
 	a = 1

	for chr in sorted(ends.keys()):
  		for strand in sorted(ends[chr].keys()):
 			for coord in sorted(ends[chr][strand].keys()):
 				if strand == '+':
 					#transcript = [chr, str(coord-1), str(coord+30), 'polyA_' + str(a), str(ends[chr][strand][coord]), strand, str(coord-1), str(coord+30), '255,0,0']
 					#todos.write('\t'.join(transcript) + '\n')
 					
					transcript = [chr, 'CBMSO', 'transcript', str(coord-30), str(coord), '.', '-', '.']
					att = ['gene_id "', 'polyA_' + str(a), '"; ', 'transcript_id "', 'polyA_' + str(a), '"; ', 'reads "', str(ends[chr][strand][coord]), '";'] 	
					todos.write('\t'.join(transcript) + '\t' + ''.join(att) + '\n')
					a = a + 1
 				else:
 					#transcript = [chr, str(coord-31), str(coord), 'polyA_' + str(a), str(ends[chr][strand][coord]), strand, str(coord-31), str(coord), '0,0,255']
 					#todos.write('\t'.join(transcript) + '\n')
 					
 					transcript = [chr, 'CBMSO', 'transcript', str(coord), str(coord+30), '.', '+', '.']
					att = ['gene_id "', 'polyA_' + str(a), '"; ', 'transcript_id "', 'polyA_' + str(a), '"; ', 'reads "', str(ends[chr][strand][coord]), '";'] 	
					todos.write('\t'.join(transcript) + '\t' + ''.join(att) + '\n')
 					a = a + 1
 	 						
 	todos.close()

##########################################################################################

###--Asignación de polyA--################################################################	
	
	gtf = open('transcripts_SL_polyA_cutted.gtf', 'w')

 	for chr in sorted(trans.keys()):
 		
 		for strand in sorted(trans[chr].keys()):
 			
 			if strand == '.': #--Me salto los que no tienen la hebra definida, que proceso más abajo
 				continue
 		
 			for id in sorted(trans[chr][strand].keys()):		
 				
 				tlen = trans[chr][strand][id][1] - trans[chr][strand][id][0]			
 				margen = tlen * 0.2 #--Admitimos buscar en el 20% final del mensajero hasta un máximo de 500 pb
 				margen_out = tlen * 0.05 #--Admitimos buscar fuera de las coordenada un 5% de la longitud mensajero hasta un máximo de 100 pb
 				
  				if margen > 500:
 					margen = 500				
 				
 				if margen_out > 100:
 					margen_out = 100
 								
 				
 				#--print chr, ':', strand, ':', id, ':', str(trans[chr][strand][id])
	
				#--temp = {coord:reads}
				temp = {} #--vaciamos los SL temporales
			
				#--Busqueda de polyA de la misma hebra que estén en el 3' del mensajero	
				for coord in sorted(ends[chr][strand].keys()):
				
					if (strand == '+') and (trans[chr][strand][id][1] - margen) < coord < (trans[chr][strand][id][1] + margen_out): #--polyA en el extremo 3' del mensajero
						temp[coord] = ends[chr][strand][coord] #--temp = {coord:reads}	
						del ends[chr][strand][coord] #--eliminar coord del diccionario para que no se repita
						
					elif (strand == '-') and (trans[chr][strand][id][0] - margen_out) < coord < (trans[chr][strand][id][0] + margen):
						temp[coord] = ends[chr][strand][coord] #--temp = {coord:reads}	
						del ends[chr][strand][coord] #--eliminar coord del diccionario para que no se repita

				#--Imprimo en el fichero de salida los transcritos que no tienen polyA detectado
				if not temp: 
					if strand == '+':
						transcript = [chr, 'CBMSO', 'transcript', str(trans[chr][strand][id][0]), str(trans[chr][strand][id][1]), '.', '+', '.']
						att1 = trans[chr][strand][id][2]
						att = ['polyA "', 'ND', '"; ', 'alt_polyA "', 'ND', '";']
						gtf.write('\t'.join(transcript) + '\t' + att1 + ''.join(att) + '\n')
						
					if strand == '-':
						transcript = [chr, 'CBMSO', 'transcript', str(trans[chr][strand][id][0]), str(trans[chr][strand][id][1]), '.', '-', '.']
						att1 = trans[chr][strand][id][2]
						att = ['polyA "', 'ND', '"; ', 'alt_polyA "', 'ND', '";']
						gtf.write('\t'.join(transcript) + '\t' + att1 + ''.join(att) + '\n')
				
				
				#--transcritos con 1 sólo polyA
				elif len(temp) == 1: 
					if strand == '+':
						ppal = sorted(temp)[0]
						transcript = [chr, 'CBMSO', 'transcript', str(trans[chr][strand][id][0]), str(ppal), '.', '+', '.']
						att1 = trans[chr][strand][id][2]
						att = ['polyA "', str(temp[ppal]), '"; ', 'alt_polyA "', 'ND', '";']
						gtf.write('\t'.join(transcript) + '\t' + att1 + ''.join(att) + '\n')
						
					if strand == '-':
						ppal = sorted(temp)[0]
						transcript = [chr, 'CBMSO', 'transcript', str(ppal), str(trans[chr][strand][id][1]), '.', '-', '.']
						att1 = trans[chr][strand][id][2]
						att = ['polyA "', str(temp[ppal]), '"; ', 'alt_polyA "', 'ND', '";']
						gtf.write('\t'.join(transcript) + '\t' + att1 + ''.join(att) + '\n')
				
				#--Más de un polyA 	
				else: 
					
					#--ppales={coord_ppal:{'reads': X, alt: [coord(reads), coord(reads)]}}		
					ppales = defaultdict(dict) #--Almacena el polya prinpipal y sus alternativos	
					
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
						transcript = [chr, 'CBMSO', 'transcript', str(trans[chr][strand][id][0]), str(ppal), '.', '+', '.']
						att1 = trans[chr][strand][id][2]
						att = ['polyA "', str(ppales[ppal]['reads']), '"; ', 'alt_polyA "', '|'.join(ppales[ppal]['alt']) , '";']
						gtf.write('\t'.join(transcript) + '\t' + att1 + ''.join(att) + '\n')
						
					if strand == '-':
						transcript = [chr, 'CBMSO', 'transcript', str(ppal), str(trans[chr][strand][id][1]), '.', '-', '.']
						att1 = trans[chr][strand][id][2]
						att = ['polyA "', str(ppales[ppal]['reads']), '"; ', 'alt_polyA "', '|'.join(ppales[ppal]['alt']), '";']
						gtf.write('\t'.join(transcript) + '\t' + att1 + ''.join(att) + '\n')			
					
				#--elimino el mensajero que ha sido modificado
				del trans[chr][strand][id]

##########################################################################################

###--Asignación de polyA (mensajeros sin hebra)--#########################################

	gtf_both = open('transcripts_polyA_both.gtf', 'w')
 	for chr in sorted(trans.keys()):

		for id in sorted(trans[chr]['.'].keys()): #--de la hebra + y - ya no quedan ids que analizar
 				
			tlen = trans[chr][strand][id][1] - trans[chr][strand][id][0]			
			margen = tlen * 0.2 #--Admitimos buscar en el 20% final del mensajero hasta un máximo de 500 pb
			margen_out = tlen * 0.05 #--Admitimos buscar fuera de las coordenada un 5% de la longitud mensajero hasta un máximo de 100 pb
			
			if margen > 500:
				margen = 500				
			
			if margen_out > 100:
				margen_out = 100
 											
			#--dict = {coord:reads}
			pplus = {}
			pminus = {}
			
			#--Buscamos polyA en ambas direcciones						
			for coord in sorted(ends[chr]['+'].keys()): #--la hebra se refiera a la de los polya			
				if (trans[chr]['.'][id][1] - margen) < coord < (trans[chr]['.'][id][1] + margen_out): #--polyA en el extremo 3' del mensajero
					pplus[coord] = ends[chr]['+'][coord] #--temp = {coord:reads}	
			
			for coord in sorted(ends[chr]['-'].keys()): #--la hebra se refiera a la de los polya	
				if (trans[chr]['.'][id][0] - margen_out) < coord < (trans[chr]['.'][id][0] + margen):
					pminus[coord] = ends[chr]['-'][coord] #--temp = {coord:reads}		
 				
 			if pplus and pminus: #--No recortamos
 			
				transcript = [chr, 'CBMSO', 'transcript', str(trans[chr]['.'][id][0]), str(trans[chr]['.'][id][1]), '.', '.', '.']
				att = trans[chr]['.'][id][2] + ' polyA "ND"; alt_polyA "ND"; cut "X";'			
				gtf_both.write('\t'.join(transcript) + '\t' + att + '\n')
				
				print 'Warning: ' + id + ' ==> PolyAs detected in both orientation (not trimmed)'
				del trans[chr]['.'][id]
 				continue
 				
 			elif pplus:
 			
 				#--Eliminamos los polyA del diccionario general
 				for p in pplus.keys():
 					del ends[chr]['+'][p]
 						
 				if len(pplus) == 1: #--Un sólo polyA
					ppal = sorted(pplus)[0]
					transcript = [chr, 'CBMSO', 'transcript', str(trans[chr]['.'][id][0]), str(ppal), '.', '+', '.']
					att1 = re.sub('X', 'F', trans[chr]['.'][id][2])
					att = ['polyA "', str(pplus[ppal]), '"; ', 'alt_polyA "', 'ND', '";']
					gtf.write('\t'.join(transcript) + '\t' + att1 + ''.join(att) + '\n')	
									
 				else: #--Más de uno
					#--ppales={coord_ppal:{'reads': X, alt: [coord(reads), coord(reads)]}}		
					ppales = defaultdict(dict) #--Almacena el polya prinpipal y sus alternativos	
					
					#--temp = [(coord, reads), ..] ordenadas por lecturas de más a menos
					temp = sorted(pplus.iteritems(), key=operator.itemgetter(1), reverse = True)
					
					#--tomamos como principal el que más se exprese
					ppal = temp[0][0]
					
					#--control de elección del ppal, por si hay algún otro con el mismo número de lecturas.
					for t in temp:
						if t[1] == pplus[ppal] and t[0] > ppal: #--El más 3'
							ppal = t[0]

					#--asigno el ppal en un diccionario
					ppales[ppal]['reads'] = pplus[ppal]
					ppales[ppal]['alt'] = []
					del pplus[ppal] #--elimino el polyA ppal del diccionario temp

					#--anoto el resto como alternativos
					for p in pplus.keys():
						ppales[ppal]['alt'].append(str(p) + '(' + str(pplus[p]) + ')')
						del pplus[p]
					
					#--escribo el transcrito en el fichero de salida
					transcript = [chr, 'CBMSO', 'transcript', str(trans[chr]['.'][id][0]), str(ppal), '.', '+', '.']
					att1 = re.sub('X', 'F', trans[chr]['.'][id][2])
					att = ['polyA "', str(ppales[ppal]['reads']), '"; ', 'alt_polyA "', '|'.join(ppales[ppal]['alt']) , '";']
					gtf.write('\t'.join(transcript) + '\t' + att1 + ''.join(att) + '\n')

				#--elimino el mensajero que ha sido modificado y rompo el bucle
				del trans[chr]['.'][id]
				
 			elif pminus:
 			
 				#--Eliminamos los polyA del diccionario general
 				for p in pminus.keys():
 					del ends[chr]['-'][p]
 						
 				if len(pminus) == 1: #--Un sólo polyA
					ppal = sorted(pminus)[0]
					transcript = [chr, 'CBMSO', 'transcript', str(ppal), str(trans[chr]['.'][id][1]), '.', '-', '.']
					att1 = re.sub('X', 'R', trans[chr]['.'][id][2])
					att = ['polyA "', str(pminus[ppal]), '"; ', 'alt_polyA "', 'ND', '";']
					gtf.write('\t'.join(transcript) + '\t' + att1 + ''.join(att) + '\n')
									
 				else: #--Más de uno
					#--ppales={coord_ppal:{'reads': X, alt: [coord(reads), coord(reads)]}}		
					ppales = defaultdict(dict) #--Almacena el polya prinpipal y sus alternativos	
					
					#--temp = [(coord, reads), ..] ordenadas por lecturas de más a menos
					temp = sorted(pminus.iteritems(), key=operator.itemgetter(1), reverse = True)
					
					#--tomamos como principal el que más se exprese
					ppal = temp[0][0]
					
					#--control de elección del ppal, por si hay algún otro con el mismo número de lecturas.
					for t in temp:
						if t[1] == pminus[ppal] and t[0] < ppal: #--El más 5'
							ppal = t[0]

					#--asigno el ppal en un diccionario
					ppales[ppal]['reads'] = pminus[ppal]
					ppales[ppal]['alt'] = []
					del pminus[ppal] #--elimino el polyA ppal del diccionario temp

					#--anoto el resto como alternativos
					for p in pminus.keys():
						ppales[ppal]['alt'].append(str(p) + '(' + str(pminus[p]) + ')')
						del pminus[p]
					
					#--escribo el transcrito en el fichero de salida
					transcript = [chr, 'CBMSO', 'transcript', str(ppal), str(trans[chr]['.'][id][1]), '.', '-', '.']
					att1 = re.sub('X', 'R', trans[chr]['.'][id][2])
					att = ['polyA "', str(ppales[ppal]['reads']), '"; ', 'alt_polyA "', '|'.join(ppales[ppal]['alt']), '";']
					gtf.write('\t'.join(transcript) + '\t' + att1 + ''.join(att) + '\n')			

				#--elimino el mensajero que ha sido modificado y rompo el bucle
				del trans[chr]['.'][id] 
							
 			else:
 				continue
				
	gtf_both.close()

##########################################################################################

###--Imprimimos el resto de mensajeros que no tienen polyA--##############################

 	for chr in sorted(trans.keys()):

		for id in sorted(trans[chr]['.'].keys()): #--de la hebra + y - ya no quedan ids que analizar
				
			transcript = [chr, 'CBMSO', 'transcript', str(trans[chr]['.'][id][0]), str(trans[chr]['.'][id][1]), '.', '.', '.']
			att1 = trans[chr]['.'][id][2]
			att = ['polyA "', 'ND', '"; ', 'alt_polyA "', 'ND' , '";']
			gtf.write('\t'.join(transcript) + '\t' + att1 + ''.join(att) + '\n')
			
			del trans[chr]['.'][id]
			
	gtf.close()

##########################################################################################

###--Imprimimos el resto de polyA no utilizados--########################################

 	poly = open ('polyA_unused.gtf', 'w')
 	
 	a = 1
 	
 	for chr in sorted(ends.keys()):
 		for strand in sorted(ends[chr].keys()):
 			if strand == '+':
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
