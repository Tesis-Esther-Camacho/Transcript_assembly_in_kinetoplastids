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


if __name__ == '__main__':

	if len(sys.argv) !=3:
		print 'Syntax: python SL_cut.py file.sam transcripts.gtf'
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
		id = col[8].split('"')[1]
		chr = col[0]
		trans[chr][id] = [int(col[3]), int(col[4])]
	
	tfile.close()
	
##########################################################################################

###--Seqlider.sam--#######################################################################	
	
	#--Obtenemos la información de los SL del fichero SAM.
	SL_file = open (sys.argv[1], 'r')
	#--starts = {chr:{strand:{coord:reads}}}
	starts = defaultdict(lambda : defaultdict(dict))

	for line in SL_file:
		line = line.rstrip('\n')

		if '@' in line[0]: #--saltamos la cabecera del sam
			continue		

		col = line.split('\t')
		chr = col[2]

		#--Obtenermos la información del bitwise flag
		bitwise = int(col[1])
		bits = bin(bitwise)[2:]
		byte = "{0:0>12}".format(bits)

		#--extraemos las coordenadas y la orientación
		if int(byte[-5]): #--int(byte[-5])==1 (reverse)
			strand = '-'
			cigar = col[5]
			let = re.split('\d+', cigar)[1:]		
			val = re.split('\D+', cigar)[:-1]
			
			coord = 0
			for l, v in zip(let, val):
				if l == 'M': #--Fragmento de lectura que hibrida.
					coord = coord + int(v)
				elif l == 'D': #--Delección en la lectura que desplaza el final de la lectura.
					coord = coord + int(v)
				else: #--Las inserciones no las sumanos, porque no modificar la coordenada final.
					continue
			coord = int(col[3]) + coord - 1

		else: #--int(byte[-5]) = 0 (forward)
			strand = '+'
			coord = int(col[3])

		#--Metemos los datos en un diccionario
		try:
			starts[chr][strand][coord] += 1
		except KeyError:
			starts[chr][strand][coord] = 1

	SL_file.close()
	
	#--Imprimimos un fichero con todos los SL detectados y el número de lecturas
	
	todos = open ('SL_detected.gtf', 'w')
	#todos.write('track name="SL" itemRgb="On"' + '\n')
	a = 1
	#chrom	chromStart	chromEnd	name	score	strand
	for chr in sorted(starts.keys()):
  		for strand in sorted(starts[chr].keys()):
 			for coord in sorted(starts[chr][strand].keys()):
 				if strand == '+':
 					#transcript = [chr, str(coord-1), str(coord+30), 'SL_' + str(a), str(starts[chr][strand][coord]), strand, str(coord-1), str(coord+30), '255,0,0']
 					#todos.write('\t'.join(transcript) + '\n')
 					
					transcript = [chr, 'CBMSO', 'transcript', str(coord), str(coord+30), '.', strand, '.']
					att = ['gene_id "', 'SL_' + str(a), '"; ', 'transcript_id "', 'SL_' + str(a), '"; ', 'reads "', str(starts[chr][strand][coord]), '";'] 	
					todos.write('\t'.join(transcript) + '\t' + ''.join(att) + '\n')
					a = a + 1
 				else:
 					#transcript = [chr, str(coord-31), str(coord), 'SL_' + str(a), str(starts[chr][strand][coord]), strand, str(coord-31), str(coord), '0,0,255']
 					#todos.write('\t'.join(transcript) + '\n')
 					
 					transcript = [chr, 'CBMSO', 'transcript', str(coord-30), str(coord), '.', strand, '.']
					att = ['gene_id "', 'SL_' + str(a), '"; ', 'transcript_id "', 'SL_' + str(a), '"; ', 'reads "', str(starts[chr][strand][coord]), '";'] 	
					todos.write('\t'.join(transcript) + '\t' + ''.join(att) + '\n')
 					a = a + 1
	todos.close()
	
 	#--filtrado de los SL de una lectura y las guardamos en un fichero
 	sl = open ('SL_descartadas.gtf', 'w')
 	
 	a = 1
 	
 	for chr in sorted(starts.keys()):
  		for strand in sorted(starts[chr].keys()):
 			for coord in sorted(starts[chr][strand].keys()):
 				if starts[chr][strand][coord] == 1:
 					if strand == '+':
 						transcript = [chr, 'CBMSO', 'transcript', str(coord), str(coord+30), '.', strand, '.']
 						att = ['gene_id "', 'SL_' + str(a), '"; ', 'transcript_id "', 'SL_' + str(a), '"; ', 'reads "', '1', '";'] 	
 						sl.write('\t'.join(transcript) + '\t' + ''.join(att) + '\n')
 					else:
 						transcript = [chr, 'CBMSO', 'transcript', str(coord-30), str(coord), '.', strand, '.']
 						att = ['gene_id "', 'SL_' + str(a), '"; ', 'transcript_id "', 'SL_' + str(a), '"; ', 'reads "', '1', '";'] 	
 						sl.write('\t'.join(transcript) + '\t' + ''.join(att) + '\n')
					del starts[chr][strand][coord]
					a = a + 1
 						
 	sl.close()
 	
##########################################################################################

###--División por SLs--###################################################################
	
	gtf = open ('transcripts_SL_cutted.gtf', 'w')	
	gtf_both = open ('transcripts_SL_both.gtf', 'w')
	
	#--Para cada cromosomas
	for chr in sorted(trans.keys()):
		
		#--Para cada transcrito
		for id in sorted(trans[chr].keys()): #--tomamos cada transcrito	
		
			i = 1 #--Para nombrar los mensajeros que vamos dividiendo
			
			#--Buscamos los SL dentro de las coordenadas del transcrito
			sl_plus = {} #--sl_plus[coord] = reads
			sl_minus = {} #--sl_minus[coord] = reads
			sl = {} #--sl[coord] = reads
			strand = ''

			#--Buscamos SL de cadena positiva
			for coord in sorted(starts[chr]['+'].keys()):
				if (trans[chr][id][0] - 50) < coord < (trans[chr][id][1] + 50): #--SL dentro de las coordenadas del transcrito (+-50 pb)
						sl_plus[coord] = starts[chr]['+'][coord] #--sl_plus[coord] = reads

			#--Buscamos SL de cadena negativa
			for coord in sorted(starts[chr]['-'].keys()):
				if (trans[chr][id][0] - 50) < coord < (trans[chr][id][1] + 50): #--SL dentro de las coordenadas del transcrito (+-50 pb)
						sl_minus[coord] = starts[chr]['-'][coord] #--sl_minus[coord] = reads

			
			#--Control de la dirección de corte		
			if sl_plus and sl_minus:
			
				""" Detectamos SL en ambas direcciones, lo que en teoría no debería ocurrir, de modo que los más probable es
				que se trate de errores del ensamblaje primario (se han unido dos transcritos por su proximidad) o por la 
				identificación de SL errores. Por lo tanto, no dividimos estos transcritos para que sean revisamos manualmente """
				
				transcript = [chr, 'CBMSO', 'transcript', str(trans[chr][id][0]), str(trans[chr][id][1]), '.', '.', '.']
				att = ['gene_id "', id, '"; ', 'transcript_id "', id , '"; ', 'SL "', 'ND', '"; ', 'Alt_SL "', 'ND', '"; ', 'cut "', 'X', '"; '] 			
				gtf_both.write('\t'.join(transcript) + '\t' + ''.join(att) + '\n')
				
				print 'Warning: ' + id + ' ==> SLs detected in both orientation (not divided)'	
				del trans[chr][id]
				continue
						
			elif sl_plus:
				sl = sl_plus
				strand = '+'
				#--Eliminar las coordenadas de los sl del diccionario starts para sólo dejar los no usados
				for s in sl_plus:
					del starts[chr][strand][s] 		
				
			elif sl_minus:
				sl = sl_minus
				strand = '-'
				#--Eliminar las coordenadas de los sl del diccionario starts para sólo dejar los no usados
				for s in sl_minus:
					del starts[chr][strand][s] 

			else: #--No sl detected, el transcrito permanece en el diccionario trans para luego ser escrito en el fichero de salida
				continue

			#--Clasificamos los inicios
			""" Primero buscamos el que más se expresa, luego buscamos en el entorno de éste (+-500 pb), y los SLs
			en dicho entorno son considerados como alternativos. Después buscamos el siguiente que más se expresa, 
			habiendo quitado el anterior y sus alternativos, y repetimos hasta que no queden más SLs dentro de las 
			coordenadas del mensajero """
			
			#--ppales={coord_ppal:{'reads': X, alt: [coord(reads), coord(reads)]}}		
			ppales = defaultdict(dict) #--Almacena los inicios principales asociados a los alternativos
	
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
				
				#--Imprimimos el último
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

	gtf_both.close()
	
##########################################################################################

###--Imprimimos los transcritos no modificados--##########################################	
	
	for chr in sorted(trans.keys()):
		for id in trans[chr].keys():		
			#--columnas y atributos gtf
			transcript = [chr, 'CBMSO', 'transcript', str(trans[chr][id][0]), str(trans[chr][id][1]), '.', '.', '.']
			att = ['gene_id "', id + '.X', '"; ', 'transcript_id "', id + '.X', '"; ', 'SL "', 'ND', '"; ', 'Alt_SL "', 'ND', '"; ']
			gtf.write('\t'.join(transcript) + '\t' + ''.join(att) + '\n')
	
	gtf.close()
	
##########################################################################################

###--Imprimimos los Sls no utilizados--###################################################	

	out = open ('SL_unused.gtf', 'w')
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


