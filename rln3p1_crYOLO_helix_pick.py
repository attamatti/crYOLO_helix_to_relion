#!/usr/bin/env python

import sys
import glob
import subprocess

yoloboxes = glob.glob('{0}/*.box'.format(sys.argv[1]))
print(yoloboxes)

def read_yolobox(boxfile):
	boxdata = open(boxfile,'r').readlines()
	parts=[]
	filename = '{0}_YOLObox.box'.format(boxfile.split('.')[0])
	outbox = open(filename,'w')
	for i in boxdata:
		if i.split()[0] == '#helix:':
			boxsize = int(i.split(',')[-1])
			print(i.strip('\n'))
		elif '#' not in i:
			outbox.write('{0:0.0f} {1:0.0f} {2} {2}\n'.format(float(i.split()[0]),float(i.split()[1]),boxsize))
	outbox.close()

for i in yoloboxes:
	read_yolobox(i)
print ('touch','/'.join(yoloboxes[0].split('/')[:-1]),'coords_suffix_YOLObox.star')
subprocess.call(['touch','{}/coords_suffix_YOLObox.box'.format('/'.join(yoloboxes[0].split('/')[:-2]))])
