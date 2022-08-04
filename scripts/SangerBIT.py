'''
#-------------------------------------------------
#	   File Name:		xxxx.py
#	   Author:		   Lei Huang
#	   Date:		   2021.10.20
#	   E-mail:		   huanglei192@mails.ucas.ac.cn
#-------------------------------------------------
'''

import os
import time
import logging
import argparse
import textwrap
import multiprocessing
import numpy as np
import pandas as pd
from Bio.Seq import Seq
from pyfaidx import Fasta
#from bokeh.plotting import figure
#from bokeh.models import ColumnDataSource, Range1d
#from bokeh.models.glyphs import Text, Rect
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast.Applications import NcbimakeblastdbCommandline

description = '''
------------------------------------------------------------------------------------------------------------------------
xxxxx
xxxxx
------------------------------------------------------------------------------------------------------------------------
'''

parse = argparse.ArgumentParser(prog='PROG', formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent(description))
parse.add_argument("--csv", help="csvinput",required=True)
parse.add_argument("--ref", help="ref",required=True)
parse.add_argument("--seq", help="seq dir",required=True)

parse.add_argument("--o", help="output dir",default='sangerOut',required=False)

args = parse.parse_args()
time0 =time.time()

csv = pd.read_csv(args.csv)
ref = Fasta(args.ref)
seq = os.path.abspath(args.seq)
outdir = args.o

# makedir outdir
if os.path.exists(outdir):
	pass
else:
	if outdir[-1] == '/':
		outdir = outdir[:-1]
	os.system('mkdir -p '+outdir+'/refSeq')


os.chdir(outdir)

# main.log
logger = logging.getLogger()
logger.setLevel(logging.INFO)
rq = time.strftime('%Y%m%d%H%M', time.localtime(time.time()))
fh = logging.FileHandler('main.log', mode='w')
fh.setLevel(logging.DEBUG)
formatter = logging.Formatter("%(asctime)s - %(levelname)s: %(message)s")
fh.setFormatter(formatter)
logger.addHandler(fh)

# refSeq length
def reflen(x):
	with open('refSeq/'+x+'_ref.fa','w') as f:
		rseq = str(ref[x][:].seq)
		f.write('>'+x+'\n'+rseq+'\n')
	return len(rseq)

csv['RefLen'] = csv['RefID'].apply(reflen)

for i,j in csv.groupby('RefID'):
	for v in range(len(j)):
		vi = j['Sample'].values[v]
		os.system('mkdir -p '+i+'/'+vi)
		with open(i+'/'+vi+'/'+vi,'w') as f:
			# Forward sequence 
			fas_f = Fasta(seq+'/'+j['Forward'].values[v])
			seq_f = str(fas_f[list(fas_f.keys())[0]][:].seq)
			# Reverse sequence
			fas_r = Fasta(seq+'/'+j['Reverse'].values[v])
			seq_r = str(Seq(str(fas_r[list(fas_r.keys())[0]][:].seq)).reverse_complement())
			f.write('>f'+str(v)+'\n'+seq_f+'\n>r'+str(v)+'\n'+seq_r+'\n')

def pool_run(t):
	sam = csv['Sample'].values[t]
	sam = csv['RefID'].values[t]+'/'+sam+'/'+sam
	logger.info('cap3 '+sam+' > '+sam+'.o')
	os.system('cap3 '+sam+' > '+sam+'.o && sync')
	if os.path.getsize(sam+'.cap.contigs') != 0:
		cmd1 = NcbimakeblastdbCommandline(dbtype='nucl',input_file=sam+'.cap.contigs')
		stdout,stderr = cmd1()
		cmd2 = NcbiblastnCommandline(query='refSeq/'+csv['RefID'].values[t]+'_ref.fa',db=sam+'.cap.contigs',outfmt=6,out=sam+'.blast.o')
		stdout,stderr = cmd2()
		if os.path.getsize(sam+'.blast.o') != 0:
			dfx = pd.read_csv(sam+'.blast.o',sep='\t',header=None)
			start = dfx[8].values[0]
			end = dfx[9].values[0]
			
			with open(sam+'_trim.fa','w') as f:
				capFas = Fasta(sam+'.cap.contigs')
				seqtrim =  Seq(str(capFas['Contig1'][end-1:start].seq)).reverse_complement()
				f.write('>'+sam.split('/')[-1]+'\n'+str(seqtrim)+'\n')
			cmd3 = NcbimakeblastdbCommandline(dbtype='nucl',input_file=sam+'_trim.fa')
			stdout,stderr = cmd3()
			cmd4 = NcbiblastnCommandline(query='refSeq/'+csv['RefID'].values[t]+'_ref.fa',db=sam+'_trim.fa',outfmt=3,out=sam+'.align.txt')
			stdout,stderr = cmd4()
			cmd5 = NcbiblastnCommandline(query='refSeq/'+csv['RefID'].values[t]+'_ref.fa',db=sam+'_trim.fa',outfmt=6,out=sam+'.fmt6.o')
			stdout,stderr = cmd5()
			cmd6 = NcbiblastnCommandline(query='refSeq/'+csv['RefID'].values[t]+'_ref.fa',db=sam+'_trim.fa',outfmt=5,out=sam+'.xml')
			stdout,stderr = cmd6()
		else:
			logger.error(sam+'.cap.contigs could not be aligned to '+csv['RefID'].values[t])
	else:
		logger.error('cap3 fail to assemble '+sam)

if __name__ == "__main__":
	process = locals()
	for i in range(len(csv)):
		process['q'+str(i)] = multiprocessing.Process(target=pool_run, args=(i,))
	for i in range(len(csv)):
		process['q'+str(i)].start()
	for i in range(len(csv)):
		process['q'+str(i)].join()	

	def cap3(x,y):
		if os.path.getsize(y+'/'+x+'/'+x+'.cap.contigs') != 0:
			return 1
		else:
			return 0

	csv['CAP3']  = csv.apply(lambda row:cap3(row['Sample'],row['RefID']),axis=1)

	def blast(x,y):
		if os.path.exists(y+'/'+x+'/'+x+'.fmt6.o'):
			return 1
		else:
			return 0

	csv['BLAST'] = csv.apply(lambda row:blast(row['Sample'],row['RefID']),axis=1)


	def blastout(x,y,col):
		if os.path.exists(y+'/'+x+'/'+x+'.fmt6.o'):
			if os.path.getsize(y+'/'+x+'/'+x+'.fmt6.o') != 0:
				tab = pd.read_csv(y+'/'+x+'/'+x+'.fmt6.o',sep='\t',header=None)
				return tab[col].values[0]
			else:
				return 0
		else:
			return 0

	csv['pident'] = csv.apply(lambda row:blastout(row['Sample'],row['RefID'],2),axis=1)
	csv['length'] = csv.apply(lambda row:blastout(row['Sample'],row['RefID'],3),axis=1)
	csv['mismatch'] = csv.apply(lambda row:blastout(row['Sample'],row['RefID'],4),axis=1)
	csv['gapopen'] = csv.apply(lambda row:blastout(row['Sample'],row['RefID'],5),axis=1)
	csv['sstart'] = csv.apply(lambda row:blastout(row['Sample'],row['RefID'],8),axis=1)
	csv['send'] = csv.apply(lambda row:blastout(row['Sample'],row['RefID'],9),axis=1)
	csv['evalue'] = csv.apply(lambda row:blastout(row['Sample'],row['RefID'],10),axis=1)	
	csv['bitscore'] = csv.apply(lambda row:blastout(row['Sample'],row['RefID'],11),axis=1)

	csv = csv.sort_values(by=['RefID','Sample'])
	csv.to_csv('result_summary.tab',sep='\t',index=None)
	csv.to_excel('result_summary.xlsx',index=None)
			
	time1=time.time()
	logger.info('CAP3 & visualization finished! Running time: %s seconds'%(round(time1-time0,2))+'.')
