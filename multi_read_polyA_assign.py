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
Usage:	polyA_1.py transcriptome_mod_3.gtf polyA_unused_mod_2.gtf
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
#--obtenemos la información de los polyA no usados
trim_file = open (sys.argv[2], 'r')

#--sl[chr][strand][[SL_id] = coord

trim = diccionarioAnidado()
for line in trim_file:
	line = line.rstrip('\n')
	col = line.split('\t')[:-1]
	att = line.split('\t')[-1].split('"')[1::2]
	reads = att[2]
	if col[6] == '+':
		trim[col[0]]['-'][int(col[3])] = reads
	else:
		trim[col[0]]['+'][int(col[4])] = reads
trim_file.close()

##########################################################################################

#--Analizamos cada transcrito para determinar si hay SL de una lectura que podamos asignar
out = open(sys.argv[1].split('.')[0][:-1] + '4.gtf', 'w')
corregidos = []
for chr in trans.keys():
	for strand in trans[chr].keys():
		for tid in trans[chr][strand].keys():

			tx =  int(trans[chr][strand][tid][0])
			ty =  int(trans[chr][strand][tid][1])
			tlen = ty-tx
			ventana = tlen * 0.30 #--Sólo buscamos al principio de los transcritos (20% de la longitud inicial)
			if ventana > 500:
				ventana = 500

#			margen = tlen * 0.05 #--y fuera del mensajero (10% de la longitud)

##########################################################################################

			if strand == '+':

				cuts = []

				for t in trim[chr][strand].keys():
					if (ty - ventana) < t < ty and int(trim[chr][strand][t]) > 1:
						cuts.append(int(t))
						#del trim[chr][strand][t] #--Eliminar elementos utilizados 

				if cuts:
					corregidos.append(tid)

					if trans[chr][strand][tid][2][9] != 'ND': #--Si ya tiene polyA asignado
						
						att = trans[chr][strand][tid][2] #--atributos del mensajero
					
						ppal = sorted(cuts)[-1] #--Más 3'
						alt = []
						
						for c in cuts:
							if int(trim[chr][strand][c]) > int(trim[chr][strand][ppal]):
								ppal = c #--Si hay alguno que se exprese más lo considero principal
						
						#--Si el principal que hemos definido se expresa más que el asignado inicialmente lo cambiamos
						if int(trim[chr][strand][ppal]) > int(att[9]):
							transcript = [chr, 'CBMSO', 'transcript', str(tx), str(ppal), '.', strand, '.']
							alt.append(str(ty) + '(' + att[9] + ')')
							att[9] = trim[chr][strand][ppal]
						
							cuts.remove(ppal)	
							del trim[chr][strand][ppal]
						
							for c in cuts:
								alt.append(str(c) + '(' + trim[chr][strand][c] + ')')
								del trim[chr][strand][c]
			
							if att[11] == 'ND':
								att[11] = '|'.join(alt)
							else:
								att[11] = att[11] + '|' + '|'.join(alt)

							out.write('\t'.join(transcript) + '\t' + '"'.join(att) + '\n')
						
						else:
						
							transcript = [chr, 'CBMSO', 'transcript', str(tx), str(ty), '.', strand, '.']

							for c in cuts:
								alt.append(str(c) + '(' + trim[chr][strand][c] + ')')
								del trim[chr][strand][c]
			
							if att[11] == 'ND':
								att[11] = '|'.join(alt)
							else:
								att[11] = att[11] + '|' + '|'.join(alt)

							out.write('\t'.join(transcript) + '\t' + '"'.join(att) + '\n')


					else: #--Sin polya asignado
					
						#--Consideramos principal al que más se exprese, y el resto como alternativos de este
						ppal = sorted(cuts)[-1] #--Más 3'
						for c in cuts:
							if int(trim[chr][strand][c]) > int(trim[chr][strand][ppal]):
								ppal = c #--Si hay alguno que se exprese más lo considero principal
						
						cuts.remove(ppal)
						
						att = trans[chr][strand][tid][2]
						transcript = [chr, 'CBMSO', 'transcript', str(tx), str(ppal), '.', strand, '.']	
						att[9] = trim[chr][strand][ppal]
						del trim[chr][strand][ppal]
							
						alt = []
						for c in cuts:
							alt.append(str(c) + '(' + trim[chr][strand][c] + ')')
							del trim[chr][strand][c]
							
						if alt:
							att[11] = '|'.join(alt)
							
						out.write('\t'.join(transcript) + '\t' + '"'.join(att) + '\n')

					#--Elimino el transcrito del diccionario
					del trans[chr][strand][tid]

##########################################################################################

			elif strand == '-':

				cuts = []

				for t in trim[chr][strand].keys():
					if (tx) < t < (tx + ventana) and int(trim[chr][strand][t]) > 1:
						cuts.append(int(t))
						#del trim[chr][strand][t] #--Eliminar elementos utilizados 

				if cuts:
					corregidos.append(tid)

					if trans[chr][strand][tid][2][9] != 'ND': #--Si ya tiene polyA asignado
						
						att = trans[chr][strand][tid][2] #--atributos del mensajero
					
						ppal = sorted(cuts)[0] #--Más 5'
						alt = []
						
						for c in cuts:
							if int(trim[chr][strand][c]) > int(trim[chr][strand][ppal]):
								ppal = c #--Si hay alguno que se exprese más lo considero principal
						
						#--Si el principal que hemos definido se expresa más que el asignado inicialmente lo cambiamos
						if int(trim[chr][strand][ppal]) > int(att[9]):
							transcript = [chr, 'CBMSO', 'transcript', str(ppal), str(ty), '.', strand, '.']
							alt.append(str(tx) + '(' + att[9] + ')')
							att[9] = trim[chr][strand][ppal]
						
							cuts.remove(ppal)	
							del trim[chr][strand][ppal]
						
							for c in cuts:
								alt.append(str(c) + '(' + trim[chr][strand][c] + ')')
								del trim[chr][strand][c]
			
							if att[11] == 'ND':
								att[11] = '|'.join(alt)
							else:
								att[11] = att[11] + '|' + '|'.join(alt)

							out.write('\t'.join(transcript) + '\t' + '"'.join(att) + '\n')
						
						else:
						
							transcript = [chr, 'CBMSO', 'transcript', str(tx), str(ty), '.', strand, '.']

							for c in cuts:
								alt.append(str(c) + '(' + trim[chr][strand][c] + ')')
								del trim[chr][strand][c]
			
							if att[11] == 'ND':
								att[11] = '|'.join(alt)
							else:
								att[11] = att[11] + '|' + '|'.join(alt)

							out.write('\t'.join(transcript) + '\t' + '"'.join(att) + '\n')


					else: #--Sin polya asignado
					
						#--Consideramos principal al que más se exprese, y el resto como alternativos de este
						ppal = sorted(cuts)[0] #--Más 5'
						for c in cuts:
							if int(trim[chr][strand][c]) > int(trim[chr][strand][ppal]):
								ppal = c #--Si hay alguno que se exprese más lo considero principal
						
						cuts.remove(ppal)
						
						att = trans[chr][strand][tid][2]
						transcript = [chr, 'CBMSO', 'transcript', str(ppal), str(ty), '.', strand, '.']	
						att[9] = trim[chr][strand][ppal]
						del trim[chr][strand][ppal]
							
						alt = []
						for c in cuts:
							alt.append(str(c) + '(' + trim[chr][strand][c] + ')')
							del trim[chr][strand][c]
							
						if alt:
							att[11] = '|'.join(alt)
							
						out.write('\t'.join(transcript) + '\t' + '"'.join(att) + '\n')

					#--Elimino el transcrito del diccionario
					del trans[chr][strand][tid]
					


##########################################################################################

			else: #--Mensajeros sin hebra
				
				#--Buscamos en forward
				cuts_plus = []
				for t in trim[chr][strand].keys():
					if (ty - ventana) < t < ty and int(trim[chr][strand][t]) > 1:
						cuts_plus.append(int(t))
				
				#--Buscamos en reverse
				cuts_minus = []
				for t in trim[chr][strand].keys():
					if (tx) < t < (tx + ventana) and int(trim[chr][strand][t]) > 1:
						cuts_minus.append(int(t))

				if cuts_plus and cuts_minus:
					continue #--No recortamos
				
				elif cuts_plus:
				
					cuts = cuts_plus
				
					corregidos.append(tid)

					if trans[chr][strand][tid][2][9] == 'ND': #--Si ya tiene polyA asignado
						print tid, 'Imposible'

					else: #--Sin polya asignado
					
						#--Consideramos principal al que más se exprese, y el resto como alternativos de este
						ppal = sorted(cuts)[-1] #--Más 3'
						for c in cuts:
							if int(trim[chr][strand][c]) > int(trim[chr][strand][ppal]):
								ppal = c #--Si hay alguno que se exprese más lo considero principal
						
						cuts.remove(ppal)
						
						att = trans[chr][strand][tid][2]
						transcript = [chr, 'CBMSO', 'transcript', str(tx), str(ppal), '.', '+', '.']	
						att[1] = att[3] = re.sub('X', 'F', att[1])
						att[9] = trim[chr][strand][ppal]
						del trim[chr][strand][ppal]
							
						alt = []
						for c in cuts:
							alt.append(str(c) + '(' + trim[chr][strand][c] + ')')
							del trim[chr][strand][c]
							
						if alt:
							att[11] = '|'.join(alt)
							
						out.write('\t'.join(transcript) + '\t' + '"'.join(att) + '\n')

					#--Elimino el transcrito del diccionario
					del trans[chr][strand][tid]					
					
				elif cuts_minus:
				
					cuts = cuts_minus
					
					corregidos.append(tid)
					
					if trans[chr][strand][tid][2][9] != 'ND': #--Si ya tiene polyA asignado
						print tid, 'Imposible'

					else: #--Sin polya asignado
			
						#--Consideramos principal al que más se exprese, y el resto como alternativos de este
						ppal = sorted(cuts)[0] #--Más 5'
						for c in cuts:
							if int(trim[chr][strand][c]) > int(trim[chr][strand][ppal]):
								ppal = c #--Si hay alguno que se exprese más lo considero principal
				
						cuts.remove(ppal)
				
						att = trans[chr][strand][tid][2]
						transcript = [chr, 'CBMSO', 'transcript', str(ppal), str(ty), '.', '-', '.']
						att[1] = att[3] = re.sub('X', 'R', att[1])
						att[9] = trim[chr][strand][ppal]
						del trim[chr][strand][ppal]
					
						alt = []
						for c in cuts:
							alt.append(str(c) + '(' + trim[chr][strand][c] + ')')
							del trim[chr][strand][c]
					
						if alt:
							att[11] = '|'.join(alt)
					
						out.write('\t'.join(transcript) + '\t' + '"'.join(att) + '\n')

					#--Elimino el transcrito del diccionario
					del trans[chr][strand][tid]		

				else:
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
#--Escribimos los polyA que no han sido utilizados
trim_out = open (sys.argv[2].split('.')[0][:-1] + '3.gtf', 'w')
i = 0
for chr in trim.keys():
	for strand in trim[chr].keys():
		for t in trim[chr][strand].keys():
			if strand == '+':
				transcript = [chr, 'CBMSO', 'transcript', str(t-30), str(t), '.', '-', '.']
				att = ['gene_id ', 'polyA_' + str(i), '; transcript_id ', 'polyA_' + str(i), '; reads ',trim[chr][strand][t], '; ']
				trim_out.write('\t'.join(transcript) + '\t' + '"'.join(att) + '\n')
				i = i + 1
			else:
				transcript = [chr, 'CBMSO', 'transcript', str(t), str(t+30), '.', '+', '.']
				att = ['gene_id ', 'polyA_' + str(i), '; transcript_id ', 'polyA_' + str(i), '; reads ',trim[chr][strand][t], '; ']
				trim_out.write('\t'.join(transcript) + '\t' + '"'.join(att) + '\n')
				i = i + 1
trim_out.close()

##########################################################################################		
print 'Multi-read polyA re-assignment:', len(corregidos)
out_ids = open ('transcripts_mod_4_ids.txt', 'w')
for i in corregidos:
	out_ids.write(i +'\n')
out_ids.close()
##########################################################################################
