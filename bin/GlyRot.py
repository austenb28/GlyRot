#!/usr/bin/env python3
import subprocess
import argparse
import os
import glob
import shutil
import math
from array import array

def GlyRotHelper(glycur,glylen,glynum,glyend,gro_in,args,tempdir):
	print('current glycan contains residues ' + str(glycur) + ' to ' + str(glycur+glylen-1))

	topfile = '../' + args.top
	tprfile = 'my_tpr.tpr'
	ndxfile = 'my_ind.ndx'
	mdpfile = 'dummy.mdp'
	xtcfile = 'trx_out.xtc'

	spargs = [args.gex, 'make_ndx', '-quiet', '-f', gro_in, '-o', ndxfile]
	p = subprocess.Popen(
		spargs,
		cwd=tempdir, 
		stdin=subprocess.PIPE,
		stdout=subprocess.PIPE, 
		stderr=subprocess.PIPE
	)
	p.stdin.write(b'del 0 - 1000\n')
	buff = 'ri 1 - ' + str(glycur-1) + '\n'
	bytes = str.encode(buff)
	p.stdin.write(bytes)
	buff = 'ri ' + str(glycur) + ' - ' + str(glycur+glylen-1) + '\n'
	bytes = str.encode(buff)
	p.stdin.write(bytes)
	if(glycur != int(glyend) - glylen + 1):
		buff = 'ri ' + str(glycur+glylen) + ' - ' + glyend + '\n'
		bytes = str.encode(buff)
		p.stdin.write(bytes)
	p.stdin.write(b'q\n')
	p.communicate()[0]
	p.stdin.close()
	p.wait()

	nstep = str(int((360 / float(args.dtheta))**3))
	
	g1 = 'r_1-' + str(glycur-1)
	g2 = 'r_' + str(glycur) + '-' + str(glycur+glylen-1)
	print('g1: ' + g1)
	print('g2: ' + g2)
	if(glycur != int(glyend) - glylen + 1):
		g3 = 'r_' + str(glycur+glylen) + '-' + str(glyend)
		print('g3: ' + g3)
		buff = ('cutoff-scheme = Group\n' + 'energygrps = '
		       + g1 + ' ' + g2 + ' ' + g3 + '\n' + 'energygrp-excl = '
		       + g1 + ' ' + g1 + ' '
		       + g2 + ' ' + g2 + ' '
		       + g3 + ' ' + g3 + ' '
		       + g1 + ' ' + g3 + ' '
		       + g2 + ' ' + g3 + '\n')
	else:
		buff = ('cutoff-scheme = Group\n' + 'energygrps = '
		       + g1 + ' ' + g2 + '\n' + 'energygrp-excl = '
		       + g1 + ' ' + g1 + ' '
		       + g2 + ' ' + g2 + '\n')

	buff +=   (  'continuation = yes\n'
		       + 'nstlist = 1\n'
		       + 'nstlog = ' + nstep + '\n'
		       + 'nstcalcenergy = 1\n'
		       + 'nstenergy = 1\n'
		       + 'epsilon-r = 1\n')

	with open(tempdir + '/' + mdpfile, 'w') as mdpwrite:
		mdpwrite.write(buff)

	spargs = [args.gex, 'grompp', '-quiet', '-f', mdpfile, '-c', gro_in,
		'-p', topfile, '-n', ndxfile, '-o', tprfile]
	p = subprocess.Popen(
		spargs, 
		cwd=tempdir, 
		stdout=subprocess.PIPE, 
		stderr=subprocess.PIPE
	)
	p.wait()

	print("Generating glycan rotation trajectory...")

	spargs = ['GlyRotHelper', '-quiet', '-s', tprfile, 
		'-n' , ndxfile, '-to', xtcfile, 
		'-dtheta', args.dtheta]

	p = subprocess.Popen(
		spargs,
		cwd=tempdir, 
		stdout=subprocess.PIPE, 
		stderr=subprocess.PIPE
	)
	p.wait()

	print("Generating energy trajectory...")

	spargs = [args.gex, 'mdrun', '-quiet', '-s', tprfile, '-rerun',
		xtcfile, '-deffnm', 'trx_out', '-nsteps', nstep]
	if args.cex != '':
		spargs = [args.cex] + spargs
	p = subprocess.Popen(
		spargs,
		cwd=tempdir, 
		stdout=subprocess.PIPE, 
		stderr=subprocess.PIPE
	)
	p.wait()

	print("Calculating minimium energy...")

	spargs = ['gmx', 'energy', '-f', 'trx_out.edr', '-o', 'trx_out.xvg']

	p = subprocess.Popen(
		spargs,
		cwd=tempdir,
		stdin=subprocess.PIPE,
		stdout=subprocess.PIPE, 
		stderr=subprocess.PIPE
	)
	p.stdin.write(b'9 0\n')
	p.communicate()[0]
	p.stdin.close()
	p.wait()

	with open(tempdir + '/trx_out.xvg', 'r') as myfile:
		energy=myfile.read()

	i = int(0)
	while(i < len(energy)):
		if(energy[i] == '\n'):
			if(energy[i+1] != '#' and energy[i+1] != '@'):
				i += 1
				break
		i += 1

	while(energy[i] == ' '):
		i += 1

	j = i
	while(energy[j] != ' '):
		j += 1

	theta_min = int(float(energy[i:j]))

	i = j
	while(energy[i] == ' '):
		i += 1

	j = i
	while(energy[j] != '\n'):
		j += 1

	min_pot = float(energy[i:j])

	i = j + 1

	while(i < len(energy)):
		while(energy[i] == ' '):
			i += 1

		j = i
		while(energy[j] != ' '):
			j += 1

		theta = int(float(energy[i:j]))

		i = j
		while(energy[i] == ' '):
			i += 1

		j = i
		while(energy[j] != '\n'):
			j += 1

		if(energy[i:j] == '-nan'):
			pot = 2*min_pot
		else:
			pot = float(energy[i:j])

		if(pot < min_pot):
			min_pot = pot
			theta_min = theta	
	
		i = j + 1

	ntheta = int((360/float(args.dtheta)))
	
	theta3 = float(theta_min % ntheta)*float(args.dtheta)
	theta2 = float((theta_min // ntheta) % ntheta)*float(args.dtheta)
	theta1 = float((theta_min // (ntheta*ntheta)) % ntheta)*float(args.dtheta)

	print('theta_min index: ' + str(theta_min) +
		  ' \ntheta1: ' + str(theta1) + 
	      ' degrees\ntheta2: '    +  str(theta2) + 
	      ' degrees\ntheta3: '    +  str(theta3) + 
	      ' degrees\nmin_pot:'    + str(min_pot/1000) + ' MJ/mol.\n')

	spargs = [args.gex, 'trjconv', '-quiet', '-s', tprfile, '-f', 
	        'trx_out.xtc', '-o', gro_in, '-dump', str(theta_min)]
	p = subprocess.Popen(
		spargs,
		cwd=tempdir,
		stdin=subprocess.PIPE,
		stdout=subprocess.PIPE, 
		stderr=subprocess.PIPE
	)
	p.stdin.write(b'0\n')
	p.communicate()[0]
	p.stdin.close()
	p.wait()

	spargs = [args.gex, 'editconf', '-quiet', '-f', gro_in, '-o', 
		gro_in, '-d', '1.2']
	p = subprocess.Popen(
		spargs,
		cwd=tempdir,
		stdout=subprocess.PIPE, 
		stderr=subprocess.PIPE
	)
	p.wait()

def is_sugar(res_name):
	if(res_name[0] == '0' or
	   res_name[0] == '1' or
	   res_name[0] == '2' or
	   res_name[0] == '3' or
	   res_name[0] == '4' or
	   res_name[0] == '5' or
	   res_name[0] == '6' or
	   res_name[0] == 'P' or
	   res_name[0] == 'Q' or
	   res_name[0] == 'R' or
	   res_name[0] == 'S' or
	   res_name[0] == 'T' or
	   res_name[0] == 'U' or
	   res_name[0] == 'V' or
	   res_name[0] == 'W' or
	   res_name[0] == 'X' or
	   res_name[0] == 'Y' or
	   res_name[0] == 'Z'):
		if(res_name[1] == 'A' or
		   res_name[1] == 'D' or
		   res_name[1] == 'R' or
		   res_name[1] == 'X' or
		   res_name[1] == 'N' or
		   res_name[1] == 'E' or
		   res_name[1] == 'L' or
		   res_name[1] == 'G' or
		   res_name[1] == 'K' or
		   res_name[1] == 'I' or
		   res_name[1] == 'M' or
		   res_name[1] == 'T' or
		   res_name[1] == 'C' or
		   res_name[1] == 'P' or
		   res_name[1] == 'B' or
		   res_name[1] == 'J' or
		   res_name[1] == 'F' or
		   res_name[1] == 'Q' or
		   res_name[1] == 'H' or
		   res_name[1] == 'O' or
		   res_name[1] == 'Z' or
		   res_name[1] == 'U' or
		   res_name[1] == 'V' or
		   res_name[1] == 'Y' or
		   res_name[1] == 'W' or
		   res_name[1] == 'S'):
			if(res_name[1] == 'K'):
				if(res_name[2] == 'A' or
				   res_name[2] == 'B'):
					return True;
				elif(res_name[2] == 'N' or
				        res_name[2] == 'O'):
					if(res_name[3] == 'A' or
					   res_name[3] == 'B'):
						return True;
			elif(res_name[1] == 'S' and res_name[2] == 'G'):
				if(res_name[3] == 'A' or
				   res_name[3] == 'B'):
					return True;
			else:
				if(res_name[2] == 'A' or
			       res_name[2] == 'B'):
					return True;
	return False;

def CountGlycans(grofile):
	glylen = array('i')
	glynum = 0
	with open(grofile, 'r') as myfile:
		buff = myfile.read()
	j=0
	j=skipline(buff,j)
	j=skipline(buff,j)
	k = j+5
	while(buff[k] == ' '):
		k+=1
	while(is_sugar(buff[k:k+3]) == False):
		j=skipline(buff,j)
		k = j+5
		while(buff[k] == ' '):
			k+=1
	glycur = int(buff[j:j+5])
	while(1):
		k = j+5
		if(is_sugar(buff[k:k+3]) == False):
			cur_res_ind -= 1
			break
		cur_res_ind = int(buff[j:j+5])
		init_res_ind = int(cur_res_ind)
		termini_count = 1
		while(1):
			k = j+5
			while(buff[k] == ' '):
				k+=1
			if(buff[k] == 'U' or 
			   buff[k] == 'V' or
			   buff[k] == 'W' or
			   buff[k] == 'X' or
			   buff[k] == 'Y' or
			   buff[k] == 'Z'):
				termini_count+=1
			elif(buff[k] == 'Q' or
				 buff[k] == 'R' or
				 buff[k] == 'S' or
				 buff[k] == 'T'):
				termini_count+=2
			elif(buff[k] == 'P'):
				termini_count+=3
			elif(buff[k] == '0'):
				termini_count-=1
			if(termini_count==0):
				break
			while(int(buff[j:j+5]) == cur_res_ind):
				j=skipline(buff,j)
			cur_res_ind = int(buff[j:j+5])
		glylen.append(cur_res_ind-init_res_ind+1)
		glynum += 1
		try:
			while(int(buff[j:j+5]) == cur_res_ind):
				j=skipline(buff,j)
			if(int(buff[j:j+5]) == cur_res_ind + 1):
				cur_res_ind += 1
		except:
			break
	return(glycur,glylen,glynum,str(cur_res_ind))

def skipline(buff,index):
	while(buff[index] != '\n'):
		index += 1
	index += 1
	return index

def main():
	parser = argparse.ArgumentParser()
	parser.add_argument(
		"-gro", 
		help= 'Input Gromacs coordinate file (.gro).',
		required=True
	)
	parser.add_argument(
		"-top", 
		help= 'Input Gromacs topology file (.top).',
		required=True
	)
	parser.add_argument(
		"-o", 
		help= 'Output Gromacs coordinate file (.gro). Default "Rotated.gro".',
		default='Rotated.gro',
		required=False
	)
	parser.add_argument(
		"-dtheta", 
		help= 'Resolution of dihedral rotations (degrees). Default "15".',
		default='15',
		required=False
	)
	parser.add_argument(
		"-gex", 
		help= 'Gromacs executable.  Default "gmx".',
		default='gmx',
		required=False
	)
	parser.add_argument(
		"-cex", 
		help= 'Cluster executable (ex. srun). Only used for mdrun.',
		default='',
		required=False
	)
	args = parser.parse_args()


	os.environ['GMX_MAXBACKUP'] = '0'
	os.environ['GMX_SUPPRESS_DUMP'] ='1'

	tempdir = '.GlyRot_temp'
	try:
		os.mkdir(tempdir)
	except:
		pass

	gro_in = 'in.gro'

	spargs = [args.gex, 'editconf', '-quiet', '-f', args.gro, '-o', 
		tempdir + '/' + gro_in, '-d', '1.2']

	p = subprocess.Popen(spargs, 
			     stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	p.wait()

	(glycur,glylen,glynum,glyend) = CountGlycans(tempdir + '/' + gro_in)
	print(glycur,glylen,glynum,glyend)
	
	for i in range(glynum):
		GlyRotHelper(glycur,glylen[i],glynum,glyend,gro_in,args,tempdir)
		glycur += glylen[i]

	shutil.copyfile(tempdir + '/' + gro_in, args.o)
	for fl in glob.glob(tempdir + '/*'):
		os.remove(fl)
	os.rmdir(tempdir)

main()




