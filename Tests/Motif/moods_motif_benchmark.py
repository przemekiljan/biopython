import MOODS
import numpy
import gzip
from  Bio import SeqIO
from Bio import Motif
import time
import StringIO
import urllib

ITER=10

t_before=time.clock()
#read in the sequence from fasta (from the web)
ftp_stream=urllib.urlopen("ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r5.21_FB2009_08/fasta/dmel-2L-chromosome-r5.21.fasta.gz")
str_data=StringIO.StringIO(ftp_stream.read())
uncomp_stream= gzip.GzipFile(fileobj=str_data) 
seq2l=SeqIO.read(uncomp_stream,"fasta")

#read a simpler motif from JASPAR pfm file (from biopython tests on the web)
http_stream=urllib.urlopen("http://github.com/biopython/biopython/raw/master/Tests/Motif/SRF.pfm")
str_data=StringIO.StringIO(http_stream.read())
m1=Motif.read(str_data,"jaspar-pfm")

# read a complex motif from JASPAR pfm file (from JASPAR DB on the web)
http_stream=urllib.urlopen("http://jaspar.genereg.net/html/DOWNLOAD/MatrixDir/JASPAR_CORE_2008/MA0010.pfm")
str_data=StringIO.StringIO(http_stream.read())
m2=Motif.read(str_data,"jaspar-pfm")



#set the alphabet to DNA for Bio.Motif
seq2l.seq.alphabet=m1.alphabet

#convert the sequence to string for MOODS 
str2l=seq2l.seq.tostring()

#extract logodds matrix from Bio.Motif instance in a proper format for MOODS
logodds1=numpy.array([map(lambda x: x[1],sorted(x.items())) for x in m1.log_odds()]).transpose().tolist()
logodds2=numpy.array([map(lambda x: x[1],sorted(x.items())) for x in m2.log_odds()]).transpose().tolist()

print "reading the sequence took ",(time.clock()-t_before)/(1.0*ITER),"seconds"

print "First motif: SRF"
t_before=time.clock()
for i in range(ITER):
    l=MOODS.search(str2l,[logodds1],[0],absolute_threshold=True)
print "MOODS calculation took ",(time.clock()-t_before)/(1.0*ITER),"seconds on average"

t_before=time.clock()
for i in range(ITER):
    l=m1.scanPWM(seq2l.seq)>0.0 # this is a numarray with scores, so we can just use the ">" to select
print "Bio.Motif fast calculation took ",(time.clock()-t_before)/(1.0*ITER),"seconds on average"

print "Second motif: Broad complex II"
t_before=time.clock()
for i in range(ITER):
    l=MOODS.search(str2l,[logodds2],[0],absolute_threshold=True)
print "MOODS calculation took ",(time.clock()-t_before)/(1.0*ITER),"seconds on average"

t_before=time.clock()
for i in range(ITER):
    l=m2.scanPWM(seq2l.seq)>0.0 # this is a numarray with scores, so we can just use the ">" to select
print "Bio.Motif fast calculation took ",(time.clock()-t_before)/(1.0*ITER),"seconds on average"



# In case you have a lot of time, you can uncomment the following
#part to see the results of pure python compared to c code

#t_before=time.clock()
#for i in range(ITER):
#    l=list(m1.search_pwm(seq2l.seq,both=False))
#print "Bio.Motif pure python (SLOW) calculation took ",(time.clock()-t_before)/(1.0*ITER),"seconds on average"

#t_before=time.clock()
#for i in range(ITER):
#    l=list(m2.search_pwm(seq2l.seq,both=False))
#print "Bio.Motif pure python (SLOW) calculation took ",(time.clock()-t_before)/(1.0*ITER),"seconds on average"
