'''
#-------------------------------------------------
#      File Name:       xxxx.py
#      Author:         Lei Huang
#      Date:           2021.10.20
#      E-mail:         huanglei192@mails.ucas.ac.cn
#-------------------------------------------------
'''

from ast import Sub
import os
import sys
import time
import logging
import argparse
import textwrap
import multiprocessing
import numpy as np
import pandas as pd
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Blast import NCBIXML
from pyfaidx import Fasta
from bokeh.plotting import figure, output_file, save
from bokeh.models import ColumnDataSource, Range1d
from bokeh.models.glyphs import Text, Rect
from bokeh.layouts import gridplot
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
parse.add_argument("--color_type", help="color_type",choices=['rainbow','single'],default='single')
parse.add_argument("--o", help="output dir",default='sangerOut',required=False)

args = parse.parse_args()
time0 =time.time()

csv = pd.read_csv(args.csv)
ref = Fasta(args.ref)
seq = os.path.abspath(args.seq)
outdir = args.o
color_type = args.color_type

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
            cmd5 = NcbiblastnCommandline(query='refSeq/'+csv['RefID'].values[t]+'_ref.fa',db=sam+'_trim.fa',outfmt=5,out=sam+'.xml')
            stdout,stderr = cmd5()
            cmd6 = NcbiblastnCommandline(query='refSeq/'+csv['RefID'].values[t]+'_ref.fa',db=sam+'_trim.fa',outfmt=6,out=sam+'.fmt6.o')
            stdout,stderr = cmd6()
        else:
            logger.error(sam+'.cap.contigs could not be aligned to '+csv['RefID'].values[t])
    else:
        logger.error('cap3 fail to assemble '+sam)

def get_single_colors(seqList):
    seqlen = len(seqList[-1])
    colorBox = list((seqlen*len(seqList))*"C")
    for pos in range(0,len(seqList[0])):
        for seqID in range(0,len(seqList)-1):
            if seqList[seqID][pos] == seqList[-1][pos]:
                colorBox[(len(seqList)-1)*seqlen+pos] = "white"
                colorBox[seqID*seqlen+pos] = "white"
            else:
                if seqList[seqID][pos] == "-" and seqList[-1][pos] != "-": 
                    ### its a deletion
                    colorBox[(len(seqList)-1)*seqlen+pos] = "white"
                    colorBox[seqID*seqlen+pos] = "blue"
                elif seqList[seqID][pos] != "-" and seqList[-1][pos] == "-":
                    ### its a insertion
                    colorBox[(len(seqList)-1)*seqlen+pos] = "white"
                    colorBox[seqID*seqlen+pos] = "green"
                else:
                    colorBox[(len(seqList)-1)*seqlen+pos] = "white"
                    colorBox[seqID*seqlen+pos] = "red"
    return colorBox

def get_rainbow_colors(seqList):
    base = [bases for sequences in list(seqList) for bases in sequences]
    clrs = {"A":"red","T":"green","C":"blue","G":"orange","-":"gray","N":"gray"}
    colorBox = [clrs[i] for i in base]
    return colorBox

def PlottingBOKEHPlot(seqs,ids,inflag,outhtml):
    if inflag == "single":
        colors = get_single_colors(seqs)
    else:
        colors = get_rainbow_colors(seqs)
    seqLen = len(seqs[0])
    seqNum = len(seqs)
    text = [i for s in list(seqs) for i in s]
    x = np.arange(1,seqLen+1)
    y = np.arange(0,seqNum,1)
    xx, yy = np.meshgrid(x,y)
    gx = xx.ravel()
    gy = yy.flatten()
    recty = gy + 0.5
    source = ColumnDataSource(data=dict(x=gx,y=gy,recty=recty,text=text,colors=colors))
    plot_height_text = seqNum*20+50
    plot_height_full = seqNum*10+50
    x_range = Range1d(0, seqLen+1, bounds='auto')
    if seqLen >= 70:
        viewlentext = 70
    else:
        viewlentext = seqLen
    tools = "xpan,xwheel_zoom,reset,save"

    if seqLen >= 300:
        viewlenfull = 300
    else:
        viewlenfull = seqLen

    ### Full alignment view
    view_range_full = (0,viewlenfull)
    p = figure(title="Full view",plot_width=800,plot_height=plot_height_full,x_range=view_range_full,y_range=ids,tools=tools,min_border=0,toolbar_location='below')
    rects = Rect(x="x",y="recty",width=1,height=1,fill_color="colors",line_color=None,fill_alpha=0.5)
    p.add_glyph(source,rects)

    ### Sequence text view
    view_range_text = (0,viewlentext)
    pe = figure(title="Sequence text view", plot_width=800,plot_height=plot_height_text,x_range=view_range_text,y_range=ids,min_border=0,tools="xpan,reset")
    glyph = Text(x="x",y="y",text="text",text_align="center",text_color="black",text_font_size="15px",text_font_style="bold")
    rects = Rect(x="x",y="recty",width=1,height=1,fill_color="colors",line_color=None,fill_alpha=0.5)
    pe.add_glyph(source,rects)
    pe.add_glyph(source,glyph)
    pe.yaxis.visible = True
    pe.xaxis.major_label_text_font_style = "bold"

    ### generate plot
    p = gridplot([[p],[pe]],toolbar_location='below')
    #print(outhtml)
    ### save p as html file 
    output_file(outhtml)
    #print("Success!")
    save(p)

def BuildUniqueIDBox(seqbox,idbox):
    print(idbox)
    UniqueTempDict = {}
    UniqueFinalDict = {}
    RankCounter = 1
    for i in range(0,len(seqbox)):
        if not seqbox[i] in UniqueTempDict:
            UniqueTempDict[seqbox[i]] = [idbox[i]]
        else:
            UniqueTempDict[seqbox[i]].append(idbox[i])
    for keys, values in UniqueTempDict.items():
        ### Set the longest idbox as the first element
        values.sort(key=len,reverse=True)
        print(values[0])
        MismatchesCount = values[0].split("%")[1]
        GapOpenCount = values[0].split("%")[2]
        UniqueFinalDict[keys] = "U"+str(RankCounter)+"_N"+str(len(values))+"_"+MismatchesCount+"_"+GapOpenCount
        RankCounter += 1
    return UniqueFinalDict

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
    logger.info('CAP3 finished! Running time: %s seconds'%(round(time1-time0,2))+'.')
    print("CAP3 finished!")    

    PlottingLocation = os.path.join(os.getcwd(),"plot/")
    if not os.path.exists(PlottingLocation):
        os.makedirs(PlottingLocation)
    logger.info("Creating directory: "+PlottingLocation)

    with open(os.path.join(os.getcwd(),"result_summary.tab"),'r') as ff:
        FileRecords = {}
        for BlastRecords in ff.readlines()[1:]:
            BLASTIdentified = BlastRecords.split('\t')[6]
            if BLASTIdentified != '0':
                RefseqID = BlastRecords.split('\t')[3]
                SampleID = BlastRecords.split('\t')[0]
                MismatchesCount = str(BlastRecords).split('\t')[9]
                GapOpens = str(BlastRecords).split('\t')[10]
                if not RefseqID in FileRecords:
                    FileRecords[RefseqID] = [SampleID+"%M"+MismatchesCount+"%G"+GapOpens]
                else:
                    FileRecords[RefseqID].append(SampleID+"%M"+MismatchesCount+"%G"+GapOpens)
    if len(FileRecords.keys()) == 0:
        logger.error('No BLAST result found!')
        sys.exit(1)
    for References in FileRecords.keys():
        ReferencePlottingSequenceBoxRaw = []
        ReferencePlottingSequenceIDRaw = []
        SubjectPatternSeqBox = []
        SubjectPatternIDBox = []
        ReferenceFASTAFile = os.path.join(os.getcwd(),"refSeq/"+References+"_ref.fa")
        ReferenceSequence = str(SeqIO.read(ReferenceFASTAFile,'fasta').seq)
        ReferenceSequenceLen = len(ReferenceSequence)
        for Samples in FileRecords[References]:
            SampleID = Samples.split('%')[0]
            BlastResultHandle = open(os.path.join(os.getcwd(),References,SampleID,SampleID+".xml"),'r')
            SampleBlastRecords = NCBIXML.read(BlastResultHandle).alignments[0].hsps[0]
            SubjectPatternRaw = str(SampleBlastRecords.sbjct)
            SubjectPatternRawLen = len(SubjectPatternRaw)
            SubjectPatternStartPos = int(SampleBlastRecords.sbjct_start)
            if SubjectPatternRawLen == ReferenceSequenceLen:
                SubjectPattern = SubjectPatternRaw
            else:
                if SubjectPatternStartPos == 1:
                    SubjectPattern = SubjectPatternRaw + "N"*(ReferenceSequenceLen-SubjectPatternRawLen)
                else:
                    SubjectPattern = "N"*(SubjectPatternStartPos-1) + SubjectPatternRaw + "N"*(ReferenceSequenceLen-SubjectPatternRawLen-SubjectPatternStartPos+1)
            SubjectPatternSeqBox.append(SubjectPattern)
            SubjectPatternIDBox.append(Samples)
        UniqueSampleDict = BuildUniqueIDBox(SubjectPatternSeqBox,SubjectPatternIDBox)
        for Ukeys, Uvalues in UniqueSampleDict.items():
            ReferencePlottingSequenceBoxRaw.append(str(Ukeys))
            ReferencePlottingSequenceIDRaw.append(str(Uvalues))
        ReferencePlottingSequenceBoxRaw.append(ReferenceSequence)
        ReferencePlottingSequenceIDRaw.append(References)
        ReferencePlottingSequenceBox = ReferencePlottingSequenceBoxRaw[::-1]
        ReferencePlottingSequenceID = ReferencePlottingSequenceIDRaw[::-1]
        PlottingBOKEHPlot(ReferencePlottingSequenceBox,ReferencePlottingSequenceID,color_type,os.path.join(PlottingLocation,References+".html"))
