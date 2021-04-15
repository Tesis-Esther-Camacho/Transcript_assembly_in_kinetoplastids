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
Usage:	SL_assign.py transcriptome_mod_2.gtf SL_discarded_mod_1.gtf
----------------------------------------------------------------------------------------------------"""
##########################################################################################
if len(sys.argv) != 3:
	print syntax
	sys.exit()
##########################################################################################	
#--obtenemos la información de los transcritos ensamblados (gtf)
gtf = open (sys.argv[1], 'r')

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
sl_file = open (sys.argv[2], 'r')

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

#--Analizamos cada transcrito para determinar si hay SL de una lectura que podamos asignar
out = open(sys.argv[1].split('.')[0][:-1] + '3.gtf', 'w')
corregidos = []
for chr in trans.keys():
	for strand in trans[chr].keys():
		for tid in trans[chr][strand].keys():

			tx =  int(trans[chr][strand][tid][0])
			ty =  int(trans[chr][strand][tid][1])
			tlen = ty-tx
			inicio = tlen * 0.2 #--Sólo buscamos al principio de los transcritos (20% de la longitud inicial, max 300)
			if inicio > 300:
				inicio = 300	
				
			margen = tlen * 0.05 #--y fuera del mensajero (5% de la longitud, max 100)
			if margen >100:
				margen = 100

##########################################################################################

			if strand == '+':
				
				starts = []
				for s in sl[chr][strand].keys():
					coord = sl[chr][strand][s]
					if tx - margen < int(coord) < tx + inicio:
						starts.append(int(coord))
						del sl[chr][strand][s] #--Eliminar SL utilizados 

				if starts:
					corregidos.append(tid)
					
					if trans[chr][strand][tid][2][5] != 'ND': #--Si ya tiene SL asignado

						att = trans[chr][strand][tid][2]
						transcript = [chr, 'CBMSO', 'transcript', str(tx), str(ty), '.', strand, '.']

						alt = []

						for s in starts:
							alt.append(str(s) + '(1)')

						if att[7] == 'ND':
							att[7] = '|'.join(alt)
						else:
							att[7] = att[7] + '|' + '|'.join(alt)

						out.write('\t'.join(transcript) + '\t' + '"'.join(att) + '\n')

					else: #--Sin SL asignado
						#--Considero que todos los SL encontrados pertenecen al mismo transcritos, es decir, no hago divisiones
						#--con estos SL
						ppal = sorted(starts)[0]
						starts.remove(ppal)

						att = trans[chr][strand][tid][2]
						transcript = [chr, 'CBMSO', 'transcript', str(ppal), str(ty), '.', strand, '.']	
						att[5] = '1'
						

						alt = []
						for s in starts:
							alt.append(str(s) + '(1)')
						
						if alt:
							att[7] = '|'.join(alt)
							
						out.write('\t'.join(transcript) + '\t' + '"'.join(att) + '\n')

					#--Elimino el transcrito del diccionario
					del trans[chr][strand][tid]

##########################################################################################

			elif strand == '-':
				
				starts = []
				for s in sl[chr][strand].keys():
					coord = sl[chr][strand][s]
					if ty - inicio < int(coord) < ty + margen:
						starts.append(int(coord))
						del sl[chr][strand][s] #--Eliminar SL utilizados 

				if starts:
					corregidos.append(tid)
					
					if trans[chr][strand][tid][2][5] != 'ND': #--Si ya tiene SL asignado
						att = trans[chr][strand][tid][2]
						transcript = [chr, 'CBMSO', 'transcript', str(tx), str(ty), '.', strand, '.']

						alt = []

						for s in starts:
							alt.append(str(s) + '(1)')

						if att[7] == 'ND':
							att[7] = '|'.join(alt)
						else:
							att[7] = att[7] + '|' + '|'.join(alt)

						out.write('\t'.join(transcript) + '\t' + '"'.join(att) + '\n')

					else: #--Sin SL asignado
						#--Considero que todos los SL encontrados pertenecen al mismo transcritos, es decir, no hago divisiones
						#--con estos SL
						ppal = sorted(starts)[-1]
						starts.remove(ppal)

						att = trans[chr][strand][tid][2]
						transcript = [chr, 'CBMSO', 'transcript', str(tx), str(ppal), '.', strand, '.']
						att[5] = '1'

						alt = []
						for s in starts:
							alt.append(str(s) + '(1)')
						
						if alt:
							att[7] = '|'.join(alt)
							
						out.write('\t'.join(transcript) + '\t' + '"'.join(att) + '\n')

					#--Elimino el transcrito del diccionario
					del trans[chr][strand][tid]

##########################################################################################

			else: #--Mensajeros sin hebra

				starts = []
				sl_plus = []
				sl_minus = []
				names_plus = []
				names_minus = []
				
				#--Buscamos SL de cadena positiva
				for s in sl[chr]['+'].keys():
					coord = sl[chr]['+'][s]
					if tx - margen < int(coord) < tx + inicio:
						sl_plus.append(int(coord))
						names_plus.append(s)
				
				#--Buscamos SL de cadena negativa
				for s in sl[chr]['-'].keys():
					coord = sl[chr]['-'][s]
					if ty - inicio < int(coord) < ty + margen:
						sl_minus.append(int(coord))
						names_minus.append(s)
				
				#del sl[chr][strand][s] #--Eliminar SL utilizados 

				#--Control de la dirección de corte		
				if sl_plus and sl_minus:
			
					""" Detectamos SL en ambas direcciones, lo que en teoría no debería ocurrir, de modo que los más probable es
					que se trate de errores del ensamblaje primario (se han unido dos transcritos por su proximidad) o por la 
					identificación de SL errores. Por lo tanto, no dividimos estos transcritos para que sean revisamos manualmente """
				
					continue
		
				elif sl_plus:
					starts = sl_plus
					corregidos.append(tid)
					
					#--Eliminar las coordenadas de los sl del diccionario
					for n in names_plus:
						del sl[chr]['+'][n] 		
						
					ppal = sorted(starts)[0]
					starts.remove(ppal)

					att = trans[chr]['.'][tid][2]
					transcript = [chr, 'CBMSO', 'transcript', str(ppal), str(ty), '.', '+', '.']	
					att[1] = att[3] = re.sub('X', 'F', att[1])
					att[5] = '1'
						
					alt = []
					for s in starts:
						alt.append(str(s) + '(1)')
						
					if alt:
						att[7] = '|'.join(alt)
							
					out.write('\t'.join(transcript) + '\t' + '"'.join(att) + '\n')

					#--Elimino el transcrito del diccionario
					del trans[chr]['.'][tid]


				elif sl_minus:
					starts = sl_minus
					corregidos.append(tid)
					
					#--Eliminar las coordenadas de los sl del diccionario starts para sólo dejar los no usados
					for n in names_minus:
						del sl[chr]['-'][n] 
						
					ppal = sorted(starts)[-1]
					starts.remove(ppal)

					att = trans[chr]['.'][tid][2]
					transcript = [chr, 'CBMSO', 'transcript', str(tx), str(ppal), '.', '-', '.']
					att[1] = att[3] = re.sub('X', 'R', att[1])
					att[5] = '1'

					alt = []
					for s in starts:
						alt.append(str(s) + '(1)')
					
					if alt:
						att[7] = '|'.join(alt)
						
					out.write('\t'.join(transcript) + '\t' + '"'.join(att) + '\n')

					#--Elimino el transcrito del diccionario
					del trans[chr]['.'][tid]

				else: #--No sl detected, el transcrito permanece en el diccionario trans para luego ser escrito en el fichero de salida
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
sl_out = open (sys.argv[2].split('.')[0][:-1] + '2.gtf', 'w')
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
print 'SL re-assignment:', len(corregidos)
out_ids = open ('transcripts_mod_3_ids.txt', 'w')
for i in corregidos:
	out_ids.write(i +'\n')
out_ids.close()
##########################################################################################
