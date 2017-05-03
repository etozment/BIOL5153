#!/usr/bin/env python3.6

#add function to find reverse compliment
#modify nuc_freq to find GC content
#something to print/store/build the CDS for each gene

#finds GC content of features in watermelon mt genome
import sys

#function to clean up sequence
def clean_seq(raw_seq):
	clean=raw_seq.upper()
	clean=clean.replace('N','')
	return clean
	
def gc_cont(sequence, sig_figs=2):
	length=len(sequence)
	gc_content=(sequence.count('C')+sequence.length('G'))/lngth*100
	return (length,round(gc_content, sig_figs))

#reverse compliment
def rev_com (sequence):
	sequence=sequence.upper()
	#compliment the bases
#A to T
	com1=sequence.replace('A','t')
#C to G
	com2=com1.replace('C','g')
#G to C
	com3=com2.replace('G','c')
#T to A
	com4=com3.replace('T','a')
#reverse
	r_com=com4[::-1].upper()
	return r_com
	
#key-gene name, value- another dictionary[key-exon number, value- exon sequence]
#gene sequences[cox1][1]='The sequence for the the first exon of cox 1'
#gene sequences[cox1][2]='The sequence for the the second exon of cox 1'
gene_sequences={}
	
#key=feature type, value=sequences of that type
feat_seq={}

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
	
	if column[2]=='CDS':
	#get gene name
		attributes=column[8].split(' ; ')
		#print (attributes[0])
		
		gene_fields=attributes[0].split(' ')
		gene_name=gene_fields[1]
		#print(gene_name)
		
#	print (column[8])
	#get exon number
		#exon_num=gene_fields[3]
		# if 'exon' in gene_fields:
# 			exon_num=gene_fields[3]
# 			print (gene_name,exon_num)
# 		else:
# 			print(gene_name)
	
	#extract/clean sequence
	fragment=dna[start-1:stop]
	
	fragment=clean_seq(fragment)
	#determine whether + or - strand
	if column[6]=='-':
		fragment=rev_com(fragment)
		
	#determine if there are multiple exons
	if 'exon' in gene_fields:
		exon_num=gene_fields[3]
		gene=gene_name+'_'+exon_num
	else:
		gene=gene_name
		
	#print(gene)
# 	print (column[3])

#sort into strings by feature type and end loop

	# if column[2] in feat_seq:
# 		feat_seq[column[2]]+=fragment
# 	else:
# 		feat_seq[column[2]]=fragment
# 
# for type, sequence in feat_seq.items():
# 	print(type+'	'+str(len(sequence)))

	
#store the sequence in gene_sequences

# 	if 'exon' in gene_fields:
# 		if gene_name in gene_sequences:
# 			gene_sequences[gene_name]+=
	if gene in gene_sequences:
		gene_sequences[gene]+=fragment
	else:
		gene_sequences[gene]=fragment
# 	
# print (gene_sequences['cox1_1'])
#create and open output file
	file_name='watermelon.nt/'+gene_name+'.fasta'
	file=open(file_name,'a')
	#append data into file
	file_name.write('>'+gene+'\n'+gene_sequences[gene])
	file.close()		

#find %of genome occupied by each feature type
# ex_pc=len(exon)/len(dna)*100
# in_pc=len(intron)/len(dna)*100
# misc_pc=len(misc_feature)/len(dna)*100
# rrna_pc=len(rrna)/len(dna)*100
# rep_pc=len(repeat)/len(dna)*100
# trna_pc=len(trna)/len(dna)*100
# 
# count=-1
# feature_type=['exon', 'intron', 'misc_feature', 'repeat', 'rrna', 'trna']
# 
# for type in [exon, intron, misc_feature, repeat, rrna, trna]:
# 	count=count+1
# 	for nucleotide in ['A', 'C', 'G', 'T']:
# 		(seq_len,nt_freq)=nuc_freq(type, nucleotide)
# 		#print (seq_len)
# 		print (feature_type[count]+'	'+str(seq_len)+'	'+str(nt_freq))
# #calculate GC content of each feature
# ex_gc=((exon.count('C')+exon.count('G'))/len(exon))*100
# in_gc=((intron.count('C')+intron.count('G'))/len(intron))*100
# misc_gc=((misc_feature.count('C')+misc_feature.count('G'))/len(misc_feature))*100
# rep_gc=((repeat.count('C')+repeat.count('G'))/len(repeat))*100
# rrna_gc=((rrna.count('C')+rrna.count('G'))/len(rrna))*100
# trna_gc=((trna.count('C')+trna.count('G'))/len(trna))*100

#print results
# print ('exon		'+str(len(exon))+'	('+str("%.1f"% ex_pc)+'%)	'+str("%.2f"% ex_gc))
# print ('intron		'+str(len(intron))+'	('+str("%.1f"% in_pc)+'%)	'+str("%.2f"% in_gc))
# print ('misc_feature	'+str(len(misc_feature))+'	('+str("%.1f"% misc_pc)+'%)	'+str("%.2f"% misc_gc))
# print ('rRNA		'+str(len(rrna))+'	('+str("%.1f"% rrna_pc)+'%)	'+str("%.2f"% rrna_gc))
# print ('repeat_region	'+str(len(repeat))+'	('+str("%.1f"% rep_pc)+'%)	'+str("%.2f"% rep_gc))
# print ('tRNA		'+str(len(trna))+'	('+str("%.1f"% trna_pc)+'%)	'+str("%.2f"% trna_gc))

fsa.close()
gff.close()