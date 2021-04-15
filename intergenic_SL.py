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
Usage:	intergenic_SL.py file.gff kineto_v0.gtf SL_discarded.gtf
----------------------------------------------------------------------------------------------------"""
##########################################################################################
if len(sys.argv) != 4:
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
	elif line.startswith("#"): #--Saltamos la cabecera
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

#--obtenemos la información de los SL descartados por tener sólo una lectura (gtf)
sl_file = open (sys.argv[3], 'r')

#--sl[chr][strand][[SL_id] = coord

sl = diccionarioAnidado()
for line in sl_file:
	line = line.rstrip('\n')
	col = line.split('\t')[:-1]
	att = line.split('\t')[-1].split('"')[1::2]
	if col[6] == '+':
		sl[col[0]]['+'][att[0]] = col[3]
	else:
		sl[col[0]]['-'][att[0]] = col[4]
sl_file.close()

##########################################################################################
#--Analizamos cada transcrito para determinar si existen 2 o más genes anotados en sus coordenadas
out = open(sys.argv[2].split('.')[0] + '_mod_1.gtf', 'w')
corregidos = []
for chr in trans.keys():
	for strand in trans[chr].keys():
		for tid in trans[chr][strand].keys():

			tx =  int(trans[chr][strand][tid][0])
			ty =  int(trans[chr][strand][tid][1])

			gids = []

			for gen in genes[chr][strand].keys():
				#--Genes contenidos enteramente dentro de las coordenadas del transcrito
				if tx < int(genes[chr][strand][gen][0]) < ty and tx < int(genes[chr][strand][gen][1]) < ty:
					gids.append(gen)

			if len(gids) > 1: #--más de un gen anotado en el transcrito						

##########################################################################################

				if strand == '+':

					temp = {} #--{gen:X}
					for g in gids:
						temp[g] = genes[chr][strand][g][0] #--Sólo coordenada de inicio para ordenar los genes			

					#--Ordenados por la coordenadas de inicio
					ordenados = sorted(temp.iteritems(), key=operator.itemgetter(1))
					ord = [o[0] for o in ordenados] #--sólo nombres ordenadas

					#--Definimos las regiones entre los genes detectados para buscar posibles SL
					regions = []
					for i in range(0,len(ord)-1):
						regions.append(genes[chr][strand][ord[i]][1] + '-' + genes[chr][strand][ord[i+1]][0])

					#--Busqueda de potenciales SL
					cuts = defaultdict(list)
					for r in regions:
						sx = int(r.split('-')[0])
						sy = int(r.split('-')[1])	
								
						#--Buscamos SL en cada región
						for s in sl[chr][strand].keys():
							coord = sl[chr][strand][s]
							if sx < int(coord) < sy:
								cuts[r].append(int(coord))
								del sl[chr][strand][s] #--Eliminar SL utilizados 
 					
 					#--Clasificamos los SL detectados en ppales y alternativos en función de las coordenadas
					ppales = defaultdict(list) 						
					
	
					if cuts:	
						corregidos.append(tid)
															
						#--Clasificamos los SL de cada región (1 ppal y el resto alternativos)
						for r in cuts.keys():			
							#--El principal es el más externo de cada región y el resto se consideran alternativos de éste
							ppales[sorted(cuts[r])[0]] = sorted(cuts[r])[1:]
										
						cortes = [p for p in sorted(ppales.keys())]

						#--Generamos el primer fragmento (que debe conservar los attributos iniciales en cuanto a SLs se refiere)
						new_att = ['gene_id ', 'ND', '; transcript_id ', 'ND', '; SL ', 'ND', '; Alt_SL ', 'ND', '; PolyA ', 'ND', '; Alt_PolyA ', 'ND', '; Reference ', 'unknown', '; Remarks ', '-', ';']
						transcript = [chr, 'CBMSO', 'transcript', str(tx), str(cortes[0]-1), '.', strand, '.']
						new_att[1] = new_att[3] = tid + '.1'
						new_att[5] = trans[chr][strand][tid][2][5]
						new_att[7] = trans[chr][strand][tid][2][7]
						out.write('\t'.join(transcript) + '\t' + '"'.join(new_att) + '\n')

						i = 2
						for c in range(0, len(cortes) - 1):
							new_att = ['gene_id ', 'ND', '; transcript_id ', 'ND', '; SL ', 'ND', '; Alt_SL ', 'ND', '; PolyA ', 'ND', '; Alt_PolyA ', 'ND', '; Reference ', 'unknown', '; Remarks ', '-', ';']
							transcript = [chr, 'CBMSO', 'transcript', str(cortes[c]), str(cortes[c+1]-1), '.', strand, '.']

							if ppales[cortes[c]]:
								alt = []
								for s in ppales[cortes[c]]:
									alt.append(str(s) + '(1)')
								alt2 = str('|'.join(alt))
								new_att[7] = alt2 #--Alt_SL

							new_att[1] = new_att[3] = tid + '.' + str(i)
							new_att[5] = '1'

							out.write('\t'.join(transcript) + '\t' + '"'.join(new_att) + '\n')
							i = i + 1

						#--Generamos el último fragmento (que debe conservar los attributos iniciales en cuanto a polya se refiere)
						new_att = ['gene_id ', 'ND', '; transcript_id ', 'ND', '; SL ', 'ND', '; Alt_SL ', 'ND', '; PolyA ', 'ND', '; Alt_PolyA ', 'ND', '; Reference ', 'unknown', '; Remarks ', '-', ';']
						transcript = [chr, 'CBMSO', 'transcript', str(cortes[-1]), str(ty), '.', strand, '.']

						if ppales[cortes[-1]]:
							alt = []
							for s in ppales[cortes[-1]]:
								alt.append(str(s) + '(1)')
							alt2 = str('|'.join(alt))
							new_att[7] = alt2 #--Alt_SL

						new_att[5] = '1'
						new_att[1] = new_att[3] = tid + '.' + str(i)
						new_att[9] = trans[chr][strand][tid][2][9]
						new_att[11] = trans[chr][strand][tid][2][11]
						out.write('\t'.join(transcript) + '\t' + '"'.join(new_att) + '\n')

						#--Elimino el transcrito del diccionario
						del trans[chr][strand][tid]

##########################################################################################

				if strand == '-':
					
					temp = {} #--{gen:X}
					for g in gids:
						temp[g] = genes[chr][strand][g][1] #--Sólo coordenada de final para ordenar los genes
					
					#--Ordenados por la coordenadas de inicio
					ordenados = sorted(temp.iteritems(), key=operator.itemgetter(1), reverse=True)
					ord = [o[0] for o in ordenados] #--sólo nombres ordenadas
						
					#--Definimos las regiones entre los genes detectados para buscar posibles SL
					regions = []
					for i in range(0,len(ord)-1):
						regions.append(genes[chr][strand][ord[i+1]][1] + '-' + genes[chr][strand][ord[i]][0])
					
					#--Busqueda de potenciales SL
					cuts = defaultdict(list)
					for r in regions:
						sx = int(r.split('-')[0])
						sy = int(r.split('-')[1])

						#--Buscamos SL en cada región
						for s in sl[chr][strand].keys():
							coord = sl[chr][strand][s]
							if sx < int(coord) < sy:
								cuts[r].append(int(coord))
								del sl[chr][strand][s] #--Eliminar SL utilizados 

					#--Clasificamos los SL detectados en ppales y alternativos en función de las coordenadas
					ppales = defaultdict(list) 
	
					if cuts:	
						corregidos.append(tid)
															
						#--Clasificamos los SL de cada región (1 ppal y el resto alternativos)
						for r in cuts.keys():			
							#--El principal es el más externo de cada región y el resto se consideran alternativos de éste
							ppales[sorted(cuts[r])[-1]] = sorted(cuts[r])[:-1]

						cortes = [p for p in sorted(ppales.keys())]

						#--Generamos el primer fragmento (que debe conservar los attributos iniciales en cuanto a polyA se refiere)
						new_att = ['gene_id ', 'ND', '; transcript_id ', 'ND', '; SL ', 'ND', '; Alt_SL ', 'ND', '; PolyA ', 'ND', '; Alt_PolyA ', 'ND', '; Reference ', 'unknown', '; Remarks ', '-', ';']
						transcript = [chr, 'CBMSO', 'transcript', str(tx), str(cortes[0]), '.', strand, '.']
						if ppales[cortes[0]]:
							alt = []
							for s in ppales[cortes[0]]:
								alt.append(str(s) + '(1)')
							alt2 = str('|'.join(alt))
							new_att[7] = alt2 #--Alt_SL						

						new_att[1] = new_att[3] = tid + '.1' #--Name
						new_att[5] = '1' #--SL
						new_att[9] = trans[chr][strand][tid][2][9] #--PolyA
						new_att[11] = trans[chr][strand][tid][2][11] #--ALt_polyA					
						out.write('\t'.join(transcript) + '\t' + '"'.join(new_att) + '\n')			

						i = 2
						for c in range(0, len(cortes) - 1):
							new_att = ['gene_id ', 'ND', '; transcript_id ', 'ND', '; SL ', 'ND', '; Alt_SL ', 'ND', '; PolyA ', 'ND', '; Alt_PolyA ', 'ND', '; Reference ', 'unknown', '; Remarks ', '-', ';']
							transcript = [chr, 'CBMSO', 'transcript', str(cortes[c]+1), str(cortes[c+1]), '.', strand, '.']
							if ppales[cortes[c+1]]:
								alt = []
								for s in ppales[cortes[c+1]]:
									alt.append(str(s) + '(1)')
								alt2 = str('|'.join(alt))
								new_att[7] = alt2 #--Alt_SL
							new_att[1] = new_att[3] = tid + '.' + str(i) #--Name
							new_att[5] = '1' #--SL

							out.write('\t'.join(transcript) + '\t' + '"'.join(new_att) + '\n')	
							i = i + 1

						#--Generamos el último fragmento (que debe conservar los attributos iniciales en cuanto a polya se refiere)
						new_att = ['gene_id ', 'ND', '; transcript_id ', 'ND', '; SL ', 'ND', '; Alt_SL ', 'ND', '; PolyA ', 'ND', '; Alt_PolyA ', 'ND', '; Reference ', 'unknown', '; Remarks ', '-', ';']
						transcript = [chr, 'CBMSO', 'transcript', str(cortes[-1]+1), str(ty), '.', strand, '.']
						new_att[1] = new_att[3] = tid + '.' + str(i)
						new_att[5] = trans[chr][strand][tid][2][5]
						new_att[7] = trans[chr][strand][tid][2][7]
						out.write('\t'.join(transcript) + '\t' + '"'.join(new_att) + '\n')

						#--Elimino el transcrito del diccionario
						del trans[chr][strand][tid]

##########################################################################################
						
				else: #--Strand '.'
					continue
					
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
#--Escribimos los SL que no han sido utilizados
sl_out = open (sys.argv[3].split('.')[0] + '_mod_1.gtf', 'w')
i = 0
for chr in sl.keys():
	for strand in sl[chr].keys():
		for s in sl[chr][strand].keys():
			if strand == '+':
				transcript = [chr, 'CBMSO', 'transcript', sl[chr][strand][s], str(int(sl[chr][strand][s])+30), '.', strand, '.']
				att = ['gene_id ', 'SL_' + str(i), '; transcript_id ', 'SL_' + str(i), '; ']
				sl_out.write('\t'.join(transcript) + '\t' + '"'.join(att) + '\n')
				i = i + 1
			else:
				transcript = [chr, 'CBMSO', 'transcript', str(int(sl[chr][strand][s])-30), sl[chr][strand][s], '.', strand, '.']
				att = ['gene_id ', 'SL_' + str(i), '; transcript_id ', 'SL_' + str(i), '; ']
				sl_out.write('\t'.join(transcript) + '\t' + '"'.join(att) + '\n')
				i = i + 1
sl_out.close()
##########################################################################################
print 'Intergenic SL:', len(corregidos)
out_ids = open ('transcripts_mod_1_ids.txt', 'w')
for i in corregidos:
	out_ids.write(i +'\n')
out_ids.close()
##########################################################################################


