#!/usr/bin/env python
#-*- coding: UTF-8 -*-

from __future__ import division
from collections import defaultdict
import sys
import subprocess
import os
import shutil
import datetime
import getopt
import re
import operator

syntax = """----------------------------------------------------------------------------------------------------
Usage:	kineto.py <step> [options] <organism_name>

Input_files: organism_name.fasta, organism_name.fastq, organism_name.gff

Steps: 	-1 (SL-PolyA search, alignment, assembly & SL division)
		-2 (PolyA division)
		-3 (Anotation, refinement & quantification)

Options: [--sp=L/T][--quals=phred33/phred64][--pair][--nt1=int][--nt2=int]
Corrections: [--correct1=file.gtf][--correct2=file.gtf]

Use -h/--help/--man to get more information
----------------------------------------------------------------------------------------------------"""

help = """------------------------------------------------------------------------------------------
Usage: 

	kineto.py <step> [options] <organism_name>

Input files:
	
	organism_name.fastq (phred 33 by default)
		organism_name_1.fastq, organism_name_2.fastq for paired-end reads (--pair)
	organism_name.fasta
	organism_name.gff (from triTrypDB)
		
Transcriptome assembling steps (organism_v0):

	-1	SeqLider search
		Polya search 
		Initial alignment (Bowtie2)
		Primary assembly (Cufflinks)
		SL division 

	###--Manual review--#### transcripts_SL_both_reviewed.gtf
	
	-2	PolyA division

	###--Manual review--#### transcripts_polyA_both_reviewed.gtf
	
	-3	Initial annotation (Blastn)
		Refinement
		Re-annotation (Blastn)
		Quantification (Cufflinks)

	Note: To execute advanced steps, previous ones must be already done

Options (for transcriptome assembly, no for corrections):

	--sp	Kinetoplastid specie (L=leishmania (default), T=Trypanosoma)
	--pair	For paired-end reads (default single end reads are expected)
	--quals	phred33/phred64 (default = phred33)
	--nt1	mininum nucleotides for SL identifications (default = 8)
	--nt2	mininum nucleotides for polyA identifications (default = 6)

Corrections (put input file in gtf format in transcriptome folder):

	The correction steps are incompatible with transcriptome steps.
	--correct1=<input_file.gtf>	quantification
	--correct2=<input_file.gtf>	annotation + quantification

	Note: The transcriptome must be already assembled in a previous run of the program, 
	and in both cases a new version will be created

System requirements:
	python 2.7.X

Included with kineto:
	bowtie-2.2.9, samtools-1.2, cufflinks-2.0.2, blastn 2.2.28+
------------------------------------------------------------------------------------------"""

#--Definimos la carpeta bin
bin_path = os.path.dirname(os.path.realpath(__file__)) + '/bin5/' 

#--Condiciones generales para los alineamientos y ensamblaje
bow = '-p 8 --np 0 --n-ceil L,0,0.02 --rdg 0,6 --rfg 0,6 --mp 6,2 --score-min L,0,-0.24'
cufflinks = 'cufflinks -L transcript -p 8 -u -b'
qual = '--phred33'
nt1 = 8
nt2 = 6
leader_file = 'leishmania.fa'
modfile = ''
max = ''
min = ''
pair = 0

def main():

	global steps
	global org
	global bowtie

	#--control de los ficheros de entrada y definicion de parámetros
	steps, org = options(sys.argv[1:])
	bowtie = 'bowtie2 ' + qual + ' ' + bow

	if pair:  #--Pareadas
		checkin(org + '_1.fastq') #--fichero reads_1.fastq
		checkin(org + '_2.fastq') #--fichero reads.fastq
		
	else: #--No pareadas
		checkin(org + '.fastq') #--fichero reads.fastq

	checkin(org + '.fasta') #--genome.fasta
	checkin(org + '.gff') #--genome.gff

#--Creamos la carpeta de salida
	checkdir ('./transcriptome')
	root = os.getcwd()
	os.chdir (root + '/transcriptome') #--Nos movemos a esta carpeta donde ejecutaremos todo el proceso

#--Registro de la ejecución
	log = open ('run_log.txt', 'a')
	log.write('#' * 90 + '\n' + str(datetime.datetime.now())[:-7] + '\n')
	log.write (' '.join(sys.argv) + ' \n' + '#' * 90 + '\n')
	log.close()

###--Transcriptoma---------------------------------------------------------------------###

###--Step: 1---------------------------------------------------------------------------###
	if '-1' in steps:
		bowindex() #--genome_name
		seqleader() #--Búsqueda de SL
		polya() #--Búsqueda de polyA
		align() #--Alineamiento
		assembly() #--Ensamblaje primario
		sl_division() #--Division por SL

###--Step: 2---------------------------------------------------------------------------###
	if '-2' in steps:
		polya_division()

###--Step: 3---------------------------------------------------------------------------###
	if '-3' in steps:
		refinement()

##########################################################################################

###--Correcciones&Othres---------------------------------------------------------------###

###--correction_1-----quantification---------------------------------------------------###
	if 'b' in steps:
		v = version()
		quantification(modfile, str(v+1))

###--correction_2-----annotation+quantification----------------------------------------###
	if 'B' in steps:
		v = version()
		annot(modfile, str(v+1))
		quantification(org + '_v' + str(v+1) + '.gtf', str(v+1))

##########################################################################################

	#--Registro final
	log = open ('run_log.txt', 'a')
	log.write('-' * 90 + '\n' + '#' * 90 + '\nFIN: ' + str(datetime.datetime.now())[:-7] + '\n')
	log.write ('#' * 90 + '\n')
	log.close()

def options(argv):
	""" control de argumentos y opciones """

	global modfile
	global nt1
	global qual
	global leader_file
	global pair
	global min
	global max
		
	#--Control de la versión de python
	pyV ='.'.join([str(x) for x in sys.version_info[:2]])
	if pyV != '2.7':
		print '-'*90 + '\n' + 'ERROR: python 2.7.X required (use --help/-h)'
		sys.exit()	

	try:
		opts, args = getopt.getopt(argv,'h123a', ['help', 'man', 'nt1=', 'nt2=', 'quals=', 'pair', 'sp=', 'correct1=', 'correct2=',])  
		steps = []
		correct = []

		for opt, val in opts:
			if opt in ('-h', '--help', '--man'):
				print help
				sys.exit()
			elif opt == '--nt1':
				nt1 = int(val)
			elif opt == '--nt2':
				nt2 = int(val)
			elif opt == '--quals':
				qual = '--' + val
			elif opt == '--sp':
				if val == 'L':
					leader_file = 'leishmania.fa'
				elif val == 'T':
					leader_file = 'trypanosoma.fa'
				else:
					print '-'*90 + '\n' + 'ERROR: species not recognize'
					print syntax
					sys.exit()
			elif opt == '--pair':
				# min = val.split('-')[0]
				# max = val.split('-')[1]
				pair = 1
			elif opt == '--correct1':
				modfile = val
				correct.append('b')
			elif opt == '--correct2':
				modfile = val
				correct.append('B')
			else:
				steps.append(opt)

		if len(args) != 1: ##--control de argumentos
			print syntax
			sys.exit()

	except getopt.GetoptError:
		print '-'*90 + '\n' + 'ERROR: erroneous options (use --help/-h)' 
		print syntax
		sys.exit()

	#--Organism name
	org = args[0]

	if len(steps) == 0 and len(correct) == 0:
		print '-'*90 + '\n' + 'ERROR: at least, one step or correction must be specified (use --help/-h)'
		print syntax
		sys.exit()

	elif len(steps) > 0 and len(correct) > 0:
		print '-'*90 + '\n' + 'ERROR: transcriptome and corrections steps cannot run simultaniously (use --help/-h)'
		print syntax
		sys.exit()

	elif len(correct) > 1: #--Sólo una cada vez
		print '-'*90 + '\n' + 'ERROR: corrections steps mus be run one by one (use --help/-h)'
		print syntax
		sys.exit()

	elif len(steps) > 1: #--Sólo se puede de una en una
		print '-'*90 + '\n' + 'ERROR: transcriptome steps mus be run one by one (use --help/-h)'
		print syntax
		sys.exit()

	elif len(steps) == 0:
		return(correct, org)

	else:
		return (steps, org)

def checkin (name): 
	""" comprobación de ficheros de entrada """
	if not os.path.exists(name):
		print "-" * 80
		print "ERROR: el fichero '" + name + "' no existe"
		print "-" * 80
		sys.exit()

def checkdir (name): 
	""" comprobación de la ruta del fichero de salida, y en el 
	caso de no existir la ruta la crea """
	if not os.path.exists(name):
		os.mkdir (name)

def bowindex (): 

	""" comprobación de la existencia del índice del genoma 
	y creación de éste en caso de no existir"""

	log = open ('run_log.txt', 'a')

	if os.path.exists(org + '.1.bt2'):
		log.write('-' * 90 + '\n')
		log.write('El índice para ' + org + ' ya existe\n')

	else:
		fd = open ("index_log", "w")
		subprocess.call (bin_path + "bowtie2/" + "bowtie2-build " + "../" + org + ".fasta " + org, shell=True, stdout = fd)
		fd.close ()

		log.write('-' * 90 + '\n')
		log.write('El índice para ' + org + ' se ha generado con éxito\n')
		
	log.close()

def bam_and_index (file, name):

	""" convierte un fichero .sam en .bam, además de ordenarlo e indexarlo """

	out = open ('temp.bam', "w")
	subprocess.call (bin_path + '/samtools/' + 'samtools ' + 'view ' + '-bS ' + file, shell = True, stdout = out)
	out.close()

	subprocess.call (bin_path + '/samtools/' + 'samtools ' + 'sort ' + 'temp.bam ' + name, shell = True)
	os.remove ('temp.bam')
	subprocess.call (bin_path + '/samtools/' + 'samtools ' + 'index ' + name + ".bam", shell = True)

def gtfclean (file):

	""" programa que elimina de un fichero gtf las líneas redundantes """

	infile = open (file, 'r')
	outfile = open ("temp", "w")

	for line in infile:
		if 'transcript' in line.split('\t')[2]:
			outfile.write(line)
		else:
			continue

	infile.close ()
	outfile.close ()
	os.remove (file)
	os.rename ("temp", file)

def seqleader ():

	#--Copiamos el fichero con la secuencia leader
	shutil.copy (bin_path + leader_file, '.')

	if pair:
		#--Registro del proceso
		log = open ('run_log.txt', 'a')
		log.write('-' * 90 + '\n')
		log.write('Seqleader\n')
		log.write('Cmd: leader_pair_v1.py ' + leader_file + ' ../' + org + '_1.fastq ../' + org + '_2.fastq ' + str(nt1) + ' 4 ../' + org + '.fasta ' + org + ' \n')
		log.write('Aligment: ' + bowtie + '\n')
		log.close()

		#--buscamos SL
		log = open ('run_log.txt', 'a')
		subprocess.call (bin_path + 'leader_pair_v1.py ' + leader_file + ' ../' + org + '_1.fastq ../' + org + '_2.fastq ' + str(nt1) + ' 4 ../' + org + '.fasta ' + org + ' "' + bowtie + '"', shell = True, stdout = log, stderr = log) 
		os.remove(leader_file)
		log.close()

	else:
		#--Registro del proceso
		log = open ('run_log.txt', 'a')
		log.write('-' * 90 + '\n')
		log.write('Seqleader\n')
		log.write('Cmd: leader_single_v1.py ' + leader_file + ' ../' + org + '.fastq ' + str(nt1) + ' 4 ../' + org + '.fasta ' + org + ' \n')
		log.write('Aligment: ' + bowtie + '\n')
		log.close()

		#--buscamos SL
		log = open ('run_log.txt', 'a')
		subprocess.call (bin_path + 'leader_single_v1.py ' + leader_file + ' ../' + org + '.fastq ' + str(nt1) + ' 4 ../' + org + '.fasta ' + org + ' "' + bowtie + '"', shell = True, stdout = log, stderr = log)
		os.remove(leader_file)
		log.close()

def polya ():

	if pair:
		#--Control de la etapa anterior
		checkin(org + '_SL_trimmed_1.fastq')
		checkin(org + '_SL_trimmed_2.fastq')

		#--Registro del proceso
		log = open ('run_log.txt', 'a')
		log.write('-' * 90 + '\n')
		log.write('PolyA\n')
		log.write('Cmd: polyA_pair_v1.py ' + org + '_SL_trimmed_1.fastq ' + org + '_SL_trimmed_2.fastq ' + str(nt2) + ' ../' + org + '.fasta ' + org + ' \n')
		log.write('Alignment: ' + bowtie + ' \n')
		log.close()

		#--buscamos SL
		log = open ('run_log.txt', 'a')
		subprocess.call (bin_path + 'polyA_pair_v1.py ' + org + '_SL_trimmed_1.fastq ' + org + '_SL_trimmed_2.fastq ' + str(nt2) + ' ../' + org + '.fasta ' + org + ' "' + bowtie + '"', shell = True, stdout = log, stderr = log)
		log.close()

	else:
		#--Control de la etapa anterior
		checkin(org + '_SL_trimmed.fastq')

		#--Registro del proceso
		log = open ('run_log.txt', 'a')
		log.write('-' * 90 + '\n')
		log.write('PolyA\n')
		log.write('Cmd: polyA_single_v1.py ' + org + '_SL_trimmed.fastq ' + str(nt2) + ' ../' + org + '.fasta ' + org + ' \n')
		log.write('Alignment: ' + bowtie + ' \n')
		log.close()

		#--buscamos SL
		log = open ('run_log.txt', 'a')
		subprocess.call (bin_path + 'polyA_single_v1.py ' + org + '_SL_trimmed.fastq ' + str(nt2) + ' ../' + org + '.fasta ' + org + ' "' + bowtie + '"', shell = True, stdout = log, stderr = log)
		log.close()

def align ():

	if pair: #--Control de la etapa anterior
		checkin(org + '_SL_polyA_trimmed_1.fastq')
		checkin(org + '_SL_polyA_trimmed_2.fastq')

		#--Registro
		log = open ('run_log.txt', 'a')
		log.write('-' * 90 + '\n')
		log.write('Alignment\n')
		log.write('Cmd: ' + bowtie + " -x " + org + " -1 " + org + "_SL_polyA_trimmed_1.fastq -2 " + org + "_SL_polyA_trimmed_2.fastq -S alignment.sam" + '\n')
		log.close()

		#--Alineamiento
		log = open ('run_log.txt', 'a')
		subprocess.call (bin_path + "bowtie2/" + bowtie + " -x " + org + " -1 " + org + "_SL_polyA_trimmed_1.fastq -2 " + org + "_SL_polyA_trimmed_2.fastq -S alignment.sam", shell = True, stderr = log)
		log.close()

	else: #--Control de la etapa anterior
		checkin(org + '_SL_polyA_trimmed.fastq')

		#--Registro
		log = open ('run_log.txt', 'a')
		log.write('-' * 90 + '\n')
		log.write('Alignment\n')
		log.write('Cmd: ' + bowtie + " -x " + org + " -U " + org +"_SL_polyA_trimmed.fastq -S alignment.sam" + '\n')
		log.close()

		#--Alineamiento
		log = open ('run_log.txt', 'a')	
		subprocess.call (bin_path + "bowtie2/" + bowtie + " -x " + org + " -U " + org +"_SL_polyA_trimmed.fastq -S alignment.sam", shell = True, stderr = log)		
		log.close()

	#--convertimos en bam e indexamos
	bam_and_index ('alignment.sam', 'alignment')

	#--Extraemos los SL y polyA del alineamiento
	if not os.path.exists('seqleader.sam') and not os.path.exists('polyA.sam'):

		if pair:
			sl = defaultdict(int)
			pol = defaultdict(int)

			#--Recuperamos los ids de las lecturas con SL
			sl_ids = open('seqleader_ids', 'r')
			for line in sl_ids:
				line = line.rstrip('\n')
				col = line.split('\t')
				sl[col[0]] = int(col[1]) #--Value = pair (1/2)
			sl_ids.close()

			#--Recuperamos los ids de las lecturas con polyA
			pol_ids = open('polyA_ids', 'r')
			for line in pol_ids:
				line = line.rstrip('\n')
				col = line.split('\t')
				pol[col[0]] = int(col[1]) #--Value = pair (1/2)
			pol_ids.close()

			sam = open('alignment.sam', 'r')
			sl_sam = open('seqleader.sam', 'w')
			pol_sam = open('polyA.sam', 'w')

			#--leemos la cabecera
			while True:
				line = sam.readline()
				line = line.rstrip('\n')
				col = line.rstrip('\n').split('\t')

				#--Imprimimos la cabecera
				sl_sam.write(line + '\n')
				pol_sam.write(line + '\n')

				if col[0] == '@PG': #--Para leer sólo la cabecera
					break

			#--buscamos SL y polyA en las lecturas
			while True:	
				line1 = sam.readline().rstrip('\n')
				line2 = sam.readline().rstrip('\n')

				if not line1: #--Final del fichero
					break

				id = line1.split('\t')[0]

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

				if id[-1] == 'L': #--La lectura contiene un SL

					#--Extraemos la secuencia que contiene el SL
					if sl[id[:-1]] == 1: #--SL detectado en el par 1
						sl_sam.write(pair1 + '\n')
					else: #--SL detectado en el par 2
						sl_sam.write(pair2 + '\n')


				elif id[-1] == 'P': #--La lectura contiene un polyA
					#--Extraemos la secuencia que contiene el SL
					if pol[id[:-1]] == 1: #--polya detectado en el par 1
						pol_sam.write(pair1 + '\n')
					else: #--polya detectado en el par 2
						pol_sam.write(pair2 + '\n')
				else:
					continue

			sl_sam.close()
			pol_sam.close()

			#--Convertimos a bam e indexamos
			bam_and_index ('seqleader.sam', 'seqleader')
			bam_and_index ('polyA.sam', 'polyA')

		else: #--No pareadas
			sl = defaultdict(int)
			pol = defaultdict(int)

			#--Recuperamos los ids de las lecturas con SL
			sl_ids = open('seqleader_ids', 'r')
			for line in sl_ids:
				line = line.rstrip('\n')
				col = line.split('\t')
				sl[col[0]] = 1
			sl_ids.close()

			#--Recuperamos los ids de las lecturas con polyA
			pol_ids = open('polyA_ids', 'r')
			for line in pol_ids:
				line = line.rstrip('\n')
				col = line.split('\t')
				pol[col[0]] = 1
			pol_ids.close()

			sam = open('alignment.sam', 'r')
			sl_sam = open('seqleader.sam', 'w')
			pol_sam = open('polyA.sam', 'w')

			#--leemos la cabecera
			while True:
				line = sam.readline()
				line = line.rstrip('\n')
				col = line.rstrip('\n').split('\t')

				#--Imprimimos la cabecera
				sl_sam.write(line + '\n')
				pol_sam.write(line + '\n')

				if col[0] == '@PG': #--Para leer sólo la cabecera
					break

			#--buscamos SL y polyA en las lecturas
			while True:
				line1 = sam.readline().rstrip('\n')

				if not line1: #--Final del fichero
					break
	
				id = line1.split('\t')[0]

				if id[-1] == 'L': # and sl[id[:-1]]: #--La lectura contiene un SL
					sl_sam.write(line1 + '\n')

				elif id[-1] == 'P': # and pol[id[:-1]]: #--La lectura contiene un polyA
					pol_sam.write(line1 + '\n')
				else:
					continue

			sl_sam.close()
			pol_sam.close()

			#--Convertimos a bam e indexamos
			bam_and_index ('seqleader.sam', 'seqleader')
			bam_and_index ('polyA.sam', 'polyA')

def assembly ():
	
	#--Control de la etapa anterior
	checkin('alignment.bam')
	
	#--Registro	
	log = open ('run_log.txt', 'a')
	log.write('-' * 90 + '\n')
	log.write('Primary_Assembly\n')
	log.write('Cmd: ' + cufflinks + ' ../' + org + '.fasta alignment.bam' + '\n')
	log.close()

	#--Ensamblamos
	log = open ('assembly_log', 'w')
	subprocess.call (bin_path + 'cufflinks/' + cufflinks + ' ../' + org + '.fasta alignment.bam', shell = True, stderr = log)
	log.close()
	
	#--Limpiamos el fichero gtf de salida de cufflinks
	gtfclean('transcripts.gtf')
	
	#--recuento
	result = open ('transcripts.gtf', 'r')
	log = open ('run_log.txt', 'a')
	log.write ('Number of transcripts: ')
	log.write (str(len(result.readlines())) + '\n')
	log.close()
	result.close()
	
	#--Elimino los ficheros que no necesito
	os.remove('isoforms.fpkm_tracking')
	os.remove('genes.fpkm_tracking')
	os.remove('skipped.gtf')
	
	#--Cambio el nombre del ensamblaje inicial de cufflinks para no perder este fichero
	os.rename('transcripts.gtf', 'primary_assembly.gtf')

def sl_division ():

	#--Control de la etapa anterior
	checkin('primary_assembly.gtf')	

	#--Registro
	log = open ('run_log.txt', 'a')
	log.write('-' * 90 + '\n')
	log.write('SL division\n')
	log.write('Cmd: SL_cut_v4.py seqleader.sam primary_assembly.gtf' + '\n')
	log.close()

	#--Ejecutamos el programa que corta los mensajeros por SLs
	log = open ('run_log.txt', 'a')
	subprocess.call (bin_path + 'SL_cut_v4.py seqleader.sam primary_assembly.gtf', shell =  True, stdout = log)
	log.close()

	#--Recuento
	log = open ('run_log.txt', 'a')

	result1 = open ('transcripts_SL_cutted.gtf', 'r')
	log.write ('Number of transcripts: ')
	log.write (str(len(result1.readlines())) + '\n')
	result1.close()

	result2 = open ('transcripts_SL_both.gtf', 'r')
	log.write ('Number of conflicted transcripts: ')
	conflicted = len(result2.readlines())
	log.write (str(conflicted) + '\n')
	result2.close()

	if conflicted == 0:
		log.write ('Not manually revision required\n')
		shutil.copy('transcripts_SL_both.gtf', 'transcripts_SL_both_reviewed.gtf')
		shutil.copy('SL_unused.gtf', 'SL_unused_mod_1.gtf')

	else:
		shutil.copy('transcripts_SL_both.gtf', 'transcripts_SL_both_reviewed.gtf')
		log.write ('Assign orientation in cut attribute (transcripts_SL_both_reviewed.gtf)\n')

	log.close()

def polya_division ():

	#--Control de la etapa anterior
	checkin('transcripts_SL_cutted.gtf')
	checkin('transcripts_SL_both_reviewed.gtf')

	#--Registro
	log = open ('run_log.txt', 'a')
	log.write('-' * 90 + '\n')
	log.write('PolyA division\n')
	log.write('Cmd1: SL_recut_v1.py SL_unused.gtf transcripts_SL_both_reviewed.gtf' + '\n')
	log.write('Cmd2: polyA_cut_v4.py polyA.sam  all_transcripts_SL_cutted.gtf' + '\n')
	log.close()

	#--Ejecutamos el programa que corta los mensajeros por SL de nuevo
	#log = open ('run_log.txt', 'a')
	subprocess.call (bin_path + 'SL_recut_v1.py SL_unused.gtf transcripts_SL_both_reviewed.gtf', shell =  True) #, stdout = log)
	subprocess.call ('cat transcripts_SL_cutted.gtf transcripts_SL_both_cutted.gtf > all_transcripts_SL_cutted.gtf', shell =  True)#, stdout = log)
	#log.close()

	#--Ejecutamos el programa que corta los mensajeros por polyAs
	log = open ('run_log.txt', 'a')
	subprocess.call (bin_path + 'polyA_cut_v4.py polyA.sam  all_transcripts_SL_cutted.gtf', shell =  True, stdout = log)
	log.close()


	#--Recuento
	log = open ('run_log.txt', 'a')

	result1 = open ('transcripts_SL_polyA_cutted.gtf', 'r')
	log.write ('Number of transcripts: ')
	log.write (str(len(result1.readlines())) + '\n')
	result1.close()

	result2 = open ('transcripts_polyA_both.gtf', 'r')
	log.write ('Number of conflicted transcripts: ')
	conflicted = len(result2.readlines())
	log.write (str(conflicted) + '\n')
	result2.close()

	if conflicted == 0:
		log.write ('Not manually revision required\n')
		shutil.copy('transcripts_polyA_both.gtf', 'transcripts_polyA_both_reviewed.gtf')
		shutil.copy('polyA_unused.gtf', 'polyA_unused_mod_1.gtf')

	else:
		shutil.copy('transcripts_polyA_both.gtf', 'transcripts_polyA_both_reviewed.gtf')
		log.write ('Assign orientation in cut attribute (transcripts_polyA_both_reviewed.gtf)\n')

	log.close()

def refinement ():

	#--Control de la etapa anterior
	checkin('transcripts_SL_polyA_cutted.gtf')
	checkin('transcripts_polyA_both_reviewed.gtf')

	#--Ejecutamos el programa que corta los mensajeros por polyA de nuevo
	log = open ('run_log.txt', 'a')
	subprocess.call (bin_path + 'polyA_recut_v1.py polyA_unused.gtf transcripts_polyA_both_reviewed.gtf', shell = True, stdout = log)
	subprocess.call ('cat transcripts_SL_polyA_cutted.gtf transcripts_polyA_both_cutted.gtf > transcripts_processed.gtf', shell = True, stdout = log)
	log.close()

	#--Anotamos
	annot('transcripts_processed.gtf', '0')

	#--Ejecutamos los programas de refinamiento
	log = open ('run_log.txt', 'a')
	subprocess.call (bin_path + 'intergenic_SL.py ../' + org + '.gff ' + org + '_v0.gtf SL_descartadas.gtf', shell = True, stdout = log)
	subprocess.call (bin_path + 'intergenic_polyA.py ../' + org + '.gff ' + org + '_v0_mod_1.gtf polyA_unused_mod_1.gtf', shell = True, stdout = log)
	subprocess.call (bin_path + 'SL_assign.py ' + org + '_v0_mod_2.gtf SL_descartadas_mod_1.gtf', shell = True, stdout = log)
	subprocess.call (bin_path + 'multi_read_polyA_assign.py ' + org + '_v0_mod_3.gtf polyA_unused_mod_2.gtf', shell = True, stdout = log)
	subprocess.call (bin_path + 'polyA_re-assign.py ' + org + '_v0_mod_4.gtf polyA_unused_mod_3.gtf', shell = True, stdout = log)

	#--re-anotación
	annot(org + '_v0_mod_5.gtf', 'X')
	os.rename(org + '_vX.gtf', org + '_v0_mod_6.gtf')

	#--remarks
	subprocess.call (bin_path + 'polycistronic_transcripts_v2.py ../' + org + '.gff ' + org + '_v0_mod_6.gtf', shell =  True, stdout = log)
	subprocess.call (bin_path + 'truncated_transcripts_v2.py ../' + org + '.gff ' + org + '_v0_mod_7.gtf', shell =  True, stdout = log)

	#--Quantificación
	quantification(org + '_v0_mod_8.gtf', 'X')
	os.rename(org + '_vX.gtf', org + '_v0_mod_9.gtf')

	log.close()

def annot (file, v):

	#--Control de la etapa anterior
	checkin(file)

	#--Registro
	log = open ('run_log.txt', 'a')
	log.write('-' * 90 + '\n')
	log.write('Annotation\n')
	log.close()

	#--creamos la carpeta para la base de datos
	checkdir('db')

	#--creamos la base de datos si no existe ya
	if not os.path.exists('./db/db_log.txt') or 'B' in steps: #--Con B rehacemos la base de datos por si ha habido cambios en el organism.gtf de referencia
		#--Extraemos la secuencia de los transcritos descritos
		log = open ('./db/db_log.txt', 'w')
		subprocess.call(bin_path + 'cufflinks/' + 'gffread -w  temp.fa -g ../' + org + '.fasta ../' + org + '.gff', shell = True)
		raw_input('ahora')
		subprocess.call(bin_path + 'makeblastdb -in temp.fa -dbtype nucl -out ./db/' + org, shell = True, stdout = log)
		os.remove('temp.fa')
		log.close()

	#--Extraemos las secuencias de los mensajeros ensamblados
	subprocess.call(bin_path + 'cufflinks/' + 'gffread -w  temp.fa -g ../' + org + '.fasta ' + file, shell = True)

	#--ejecutamos blast
	subprocess.call(bin_path + 'blastn -query temp.fa -db ./db/' + org + ' -outfmt 6 -out blastout -max_target_seqs 1', shell = True)
	#-- fields: query id, subject id, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
	os.remove('temp.fa')

	#--Recuperamos los datos del blast en un diccionario
	blast =  open ('blastout', 'r')
	blastr = {}
	for i in blast:
		i = i.rstrip('\n')
		blast_col = i.split('\t')
		if float(blast_col[2]) == 100:
			blastr[blast_col[0]] = [blast_col[1], blast_col[8], blast_col[9]] #--subject/subject start/subject end
		else:
			continue
	blast.close()
	os.remove('blastout')

	#--Analizamos el fichero de los transcritos y metemos las anotaciones
	gtf = open (file, 'r')
	gtf_out = open ('temp', 'w')
	
	anotados = 0
	for i in gtf:
		i = i.rstrip('\n')
		col = i.split('\t')[:8]
		att = i.split('\t')[8].split('"')[1::2]
		att_names = i.split('\t')[8].split('"')[0::2]

		if len(att) == 6: 
			att.append('unknown')
			att.append('-')
			att_names = ['gene_id ', '; transcript_id ', '; SL ', '; Alt_SL ', '; PolyA ', '; Alt_PolyA ', '; Reference ', '; Remarks ']

		try:
			if blastr[att[0]]:
				anotados = anotados + 1
				annot = blastr[att[0]]
				att[6] = annot[0]

				if col[6] == '.': #--Mensajeros sin hebra asignada
				# Los mensajeros sin hebra definida son extraídos por gffread como F
				# de modo que si las coordenadas del subject (blastn) son directas, es decir, 
				# la resta de start-end debe ser negativo (ej: 1-999), entonces el mensajero 
				# será F y viceversa 

					if int(annot[1])-int(annot[2]) < 0:
						att[0] = att[1] = re.sub('X', 'F', att[0])
						col[6] = '+'
					else:
						att[0] = att[1] = re.sub('X', 'R', att[0])
						col[6] = '-'

				gtf_out.write('\t'.join(col) + '\t')
				z = [x + '"' + y for x,y in zip(att_names,att)]
				gtf_out.write('"'.join(z) + '";\n')

		except KeyError:
			gtf_out.write('\t'.join(col) + '\t')
			z = [x + '"' + y for x,y in zip(att_names,att)]
			gtf_out.write('"'.join(z) + '";\n')

	gtf.close()
	gtf_out.close()

	log = open ('run_log.txt', 'a')
	log.write ('Number of annotated transcripts: ')
	log.write (str(anotados) + '\n')
	log.close()

	#--Creamos en fichero de salida
	os.rename ('temp', org + '_v' + v + '.gtf')

def quantification (inf, v):
	
	#--Control de la etapa anterior
	checkin(inf)


	#--Entrada del log
	log = open ('run_log.txt', 'a')
	log.write('-' * 90 + '\n')
	log.write('7-Quantification\n')
	log.write('Cmd: ' + cufflinks + ' -G ' + inf +  ' alignment.bam' + '\n')
	log.close()
	
	#--cuantificamos
	log = open ('quantification_log', 'w')
	subprocess.call (bin_path +'cufflinks/' + cufflinks + ' ../' + org + '.fasta' + ' -G ' + inf +  ' alignment.bam', shell =  True, stderr = log)
	log.close()
	
	#--Recorremos el fichero isoforms.fpkm_tracking para obtener la cuantificación
	file = open ('isoforms.fpkm_tracking', 'r')
	cabecera = file.readline() #--la leo y ya no entra en el bucle
		
	expr = {}
	for i in file:
		i = i.rstrip('\n')
		col = i.split('\t')
		expr[col[0]] = [col[8], col[9], str((float(col[11])-float(col[9]))/2)] #--expr[id] = cov, FPKM, SD
	file.close()
	
	#--Eliminamos los ficheros que no necesitamos
	os.remove('transcripts.gtf')
	os.remove('genes.fpkm_tracking')
	os.remove('isoforms.fpkm_tracking')
	os.remove('skipped.gtf')
	
	#--Analizamos el fichero de los transcritos y metemos las datos de la cuantificación
	
	gtf = open (inf, 'r')
	gtf_out = open ('temp.gtf', 'w')
	
	for i in gtf:
		i = i.rstrip('\n')
		col = i.split('\t')[:8]
		att = i.split('\t')[8].split('"')[1::2][0:8] #--Sólo los 8 primeros
		att_names = ['gene_id' , 'transcript_id', 'SL', 'Alt_SL', 'PolyA', 'Alt_PolyA', 'Reference', 'Remarks', 'Coverage', 'FPKM', 'FPKM_SD']
		
		try:
			if expr[att[0]]:
				att = att + expr[att[0]]		
 				
 				gtf_out.write('\t'.join(col) + '\t')
 				for name, at in zip (att_names, att):
 					gtf_out.write(name + ' "' + at + '"; ')
 				gtf_out.write('\n')

		except KeyError:
			print i
			print 'ERROR FATAL'

	gtf.close()
	gtf_out.close()
	os.rename('temp.gtf', org + '_v' + v + '.gtf')

def version():

	""" Control de la última versión en formato gtf"""

	versiones = [0,]
	for (path, dirs, files) in os.walk('.'):
		for file in files:
			if file[-3:] == 'gtf':
				if org in file and not 'mod' in file:
					temp1 = file.split('.')[0][-1]
					temp2 = file.split('.')[0][-2]
					if temp1.isdigit() and temp2 == 'v':
						versiones.append(int(file.split('.')[0][-1]))

	ver = sorted(versiones)[-1] #--versión actual (gtf)
	return (ver)

if __name__ == "__main__":

	main()
