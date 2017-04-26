#!/usr/bin/env python3.6

#finds GC content of features in watermelon mt genome
import sys

#function to clean up sequence
def clean_seq(raw_seq):
	clean=raw_seq.upper()
	clean=clean.replace('N','')
	return clean
	
def nuc_freq(sequence, base, sig_figs=2):
	lngth=len(sequence)
	nuc_count=sequence.count(base)
	frequency=nuc_count/lngth*100
	return (lngth,round(frequency, sig_figs))

#set up usage statement
usage=sys.argv[0]+' genome.fasta features.gff'
if len(sys.argv)<3:
	print (usage)
	sys.exit('\nFasta and GFF files required\n')
#open the files
fsa=open (sys.argv[1])
gff=open (sys.argv[2])

#for line in fsa:
	#print(line)
#remove fasta header using line counter
count=0
dna=''
for line in fsa:
	if count>0:
		dna=dna+line
	count+=1
#print (dna)
#declare variables
exon=''
intron=''
misc_feature=''
rrna=''
repeat=''
trna=''
#get start and stop sites
for line in gff:
#clean up data
	# clean_seq(line)
# 	print (clean)
# 	sys.exit()
	column=line.split('	')
	start=int(column[3])
	stop=int(column[4])
	
	fragment=dna[start-1:stop]
	
	fragment=clean_seq(fragment)
# 	print (column[3])

#sort into strings by feature type and end loop


	
	#get exons
	if column[2]=='CDS':
# 		print (len(dna[start-1:stop]))
# 		print (column)
		exon=exon+fragment
# 		print (exon)
# 	print (len(exon))
	#get introns
	if column[2]=='intron':
		intron=intron+fragment
		#print (intron)
	#get misc_features
	if column[2]=='misc_feature':
		misc_feature=misc_feature+fragment
	if column[2]=='rRNA':
		rrna=rrna+fragment
	if column[2]=='repeat_region':
		repeat=repeat+fragment
	if column[2]=='tRNA':
		trna=trna+fragment
# print(len(exon))
#print(rrna)

#find %of genome occupied by each feature type
ex_pc=len(exon)/len(dna)*100
in_pc=len(intron)/len(dna)*100
misc_pc=len(misc_feature)/len(dna)*100
rrna_pc=len(rrna)/len(dna)*100
rep_pc=len(repeat)/len(dna)*100
trna_pc=len(trna)/len(dna)*100

count=-1
feature_type=['exon', 'intron', 'misc_feature', 'repeat', 'rrna', 'trna']

for type in [exon, intron, misc_feature, repeat, rrna, trna]:
	count=count+1
	for nucleotide in ['A', 'C', 'G', 'T']:
		(seq_len,nt_freq)=nuc_freq(type, nucleotide)
		#print (seq_len)
		print (feature_type[count]+'	'+str(seq_len)+'	'+str(nt_freq))
#calculate GC content of each feature
ex_gc=((exon.count('C')+exon.count('G'))/len(exon))*100
in_gc=((intron.count('C')+intron.count('G'))/len(intron))*100
misc_gc=((misc_feature.count('C')+misc_feature.count('G'))/len(misc_feature))*100
rep_gc=((repeat.count('C')+repeat.count('G'))/len(repeat))*100
rrna_gc=((rrna.count('C')+rrna.count('G'))/len(rrna))*100
trna_gc=((trna.count('C')+trna.count('G'))/len(trna))*100

#print results
# print ('exon		'+str(len(exon))+'	('+str("%.1f"% ex_pc)+'%)	'+str("%.2f"% ex_gc))
# print ('intron		'+str(len(intron))+'	('+str("%.1f"% in_pc)+'%)	'+str("%.2f"% in_gc))
# print ('misc_feature	'+str(len(misc_feature))+'	('+str("%.1f"% misc_pc)+'%)	'+str("%.2f"% misc_gc))
# print ('rRNA		'+str(len(rrna))+'	('+str("%.1f"% rrna_pc)+'%)	'+str("%.2f"% rrna_gc))
# print ('repeat_region	'+str(len(repeat))+'	('+str("%.1f"% rep_pc)+'%)	'+str("%.2f"% rep_gc))
# print ('tRNA		'+str(len(trna))+'	('+str("%.1f"% trna_pc)+'%)	'+str("%.2f"% trna_gc))

fsa.close()
gff.close()