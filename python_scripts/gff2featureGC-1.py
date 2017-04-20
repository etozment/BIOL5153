#!/usr/bin/env python3.6

#finds GC content of features in watermelon mt genome
import sys

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
	column=line.split('	')
	start=int(column[3])
	stop=int(column[4])
# 	print (column[3])

#sort into strings by feature type and end loop
	#get exons
	if column[2]=='CDS':
# 		print (len(dna[start-1:stop]))
# 		print (column)
		exon=exon+dna[start-1:stop]
# 		print (exon)
# 	print (len(exon))
	#get introns
	if column[2]=='intron':
		intron=intron+dna[start-1:stop]
		#print (intron)
	#get misc_features
	if column[2]=='misc_feature':
		misc_feature=misc_feature+dna[start-1:stop]
	if column[2]=='rRNA':
		rrna=rrna+dna[start-1:stop]
	if column[2]=='repeat_region':
		repeat=repeat+dna[start-1:stop]
	if column[2]=='tRNA':
		trna=trna+dna[start-1:stop]
# print(len(exon))
#print(rrna)

#find %of genome occupied by each feature type
ex_pc=len(exon)/len(dna)*100
in_pc=len(intron)/len(dna)*100
misc_pc=len(misc_feature)/len(dna)*100
rrna_pc=len(rrna)/len(dna)*100
rep_pc=len(repeat)/len(dna)*100
trna_pc=len(trna)/len(dna)*100

#calculate GC content of each feature
ex_gc=((exon.count('C')+exon.count('G'))/len(exon))*100
in_gc=((intron.count('C')+intron.count('G'))/len(intron))*100
misc_gc=((misc_feature.count('C')+misc_feature.count('G'))/len(misc_feature))*100
rep_gc=((repeat.count('C')+repeat.count('G'))/len(repeat))*100
rrna_gc=((rrna.count('C')+rrna.count('G'))/len(rrna))*100
trna_gc=((trna.count('C')+trna.count('G'))/len(trna))*100

#print results
print ('exon		'+str(len(exon))+'	('+str("%.1f"% ex_pc)+'%)	'+str("%.2f"% ex_gc))
print ('intron		'+str(len(intron))+'	('+str("%.1f"% in_pc)+'%)	'+str("%.2f"% in_gc))
print ('misc_feature	'+str(len(misc_feature))+'	('+str("%.1f"% misc_pc)+'%)	'+str("%.2f"% misc_gc))
print ('rRNA		'+str(len(rrna))+'	('+str("%.1f"% rrna_pc)+'%)	'+str("%.2f"% rrna_gc))
print ('repeat_region	'+str(len(repeat))+'	('+str("%.1f"% rep_pc)+'%)	'+str("%.2f"% rep_gc))
print ('tRNA		'+str(len(trna))+'	('+str("%.1f"% trna_pc)+'%)	'+str("%.2f"% trna_gc))

fsa.close()
gff.close()