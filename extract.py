def extract_splice(ingtf,bed3p,bed5p,outfile):
	yui=open(ingtf,'r')
	new5=open(bed5p,'a')
	new3=open(bed3p,'a')
	new=open(outfile,'a') 
	exon_list=[]
	tt,ts=0,0
	for line in yui:
		split=line.split()
		if len(split)<6:
			continue
		if split[2]=='transcript':
			tt+=1
			if len(exon_list)>1:
				ts+=1
				exon_list.sort()
				if last_strand=='+':
					for i in range(1,len(exon_list),1):
						new5.write(chro+'	'+str(exon_list[i-1][1])+'	'+str(exon_list[i-1][1]+2)+'	'+tid+'_i'+str(i)+'	.	'+last_strand+'\n')
						new3.write(chro+'	'+str(exon_list[i][0]-3)+'	'+str(exon_list[i][0]-1)+'	'+tid+'_i'+str(i)+'	.	'+last_strand+'\n')
				elif last_strand=='-':
					for i in range(1,len(exon_list),1):
						new5.write(chro+'	'+str(exon_list[i][0]-3)+'	'+str(exon_list[i][0]-1)+'	'+tid+'_i'+str(i)+'	.	'+last_strand+'\n')
						new3.write(chro+'	'+str(exon_list[i-1][1])+'	'+str(exon_list[i-1][1]+2)+'	'+tid+'_i'+str(i)+'	.	'+last_strand+'\n')
			exon_list=[]
			for i in range(6,len(split),1):
				if split[i]=='transcript_id':
					tid=split[i+1][1:-2]
		elif split[2]=='exon':
			exon_list.append([int(split[3]),int(split[4])])
		chro=split[0]
		last_strand=split[6]
		
	if len(exon_list)>1:
		ts+=1
		exon_list.sort()
		for i in range(1,len(exon_list),1):
			if last_strand=='+':
				new5.write(chro+'	'+str(exon_list[i-1][1])+'	'+str(exon_list[i-1][1]+2)+'	'+tid+'_i'+str(i)+'	.	'+last_strand+'\n')
				new3.write(chro+'	'+str(exon_list[i][0]-3)+'	'+str(exon_list[i][0]-1)+'	'+tid+'_i'+str(i)+'	.	'+last_strand+'\n')
			elif last_strand=='-':
				new5.write(chro+'	'+str(exon_list[i][0]-3)+'	'+str(exon_list[i][0]-1)+'	'+tid+'_i'+str(i)+'	.	'+last_strand+'\n')
				new3.write(chro+'	'+str(exon_list[i-1][1])+'	'+str(exon_list[i-1][1]+2)+'	'+tid+'_i'+str(i)+'	.	'+last_strand+'\n')
	yui.close()
	print ('Total:',tt)
	print ('	Spliced:',ts)
	new.close()
	new5.close()
	new3.close()


def count_dinucleotide(infile,outfile):
	seq_dict=dict()
	import os
	lc=0
	with open(infile,'r') as current:
		for line in current:
			lc+=1
			split=line.split()
			split[-1]=split[-1].upper()
			if split[-1] not in seq_dict:
				seq_dict[split[-1]]=0
			seq_dict[split[-1]]+=1
	print (lc,'splices with',len(seq_dict),'dinucleotides')
	tseq_list=[[seq_dict[seq],seq] for seq in seq_dict]
	tseq_list.sort(reverse=True)
	new=open(outfile,'a')
	for nuc in tseq_list:
		new.write(nuc[1]+'	'+str(nuc[0])+'	'+str(float(nuc[0])*100/lc)+'\n')
	new.close()


def paired_dinucleotides(splice5p,splice3p,outfile):
	this_dict=dict()
	import os
	with open(splice5p,'r') as current:
		lc=0
		for line in current:
			lc+=1
			split=line.split()
			split[-1]=split[-1].upper()
			tid=split[0].split('::')[0]
			if tid not in this_dict:
				this_dict[tid]={5:'',3:''}
			this_dict[tid][5]=split[-1]
		print ('5 prime splices:',lc)
	with open(splice3p,'r') as current:
		lc=0
		for line in current:
			lc+=1
			split=line.split()
			split[-1]=split[-1].upper()
			tid=split[0].split('::')[0]
			if tid not in this_dict:
				this_dict[tid]={5:'',3:''}
			this_dict[tid][3]=split[-1]
		print ('5 prime splices:',lc)
	count_dict=dict()
	tc,pc=0,0
	for tid in this_dict:
		tc+=1
		if len(this_dict[tid])==2:
			junc=this_dict[tid][5]+'--'+this_dict[tid][3]
			if junc not in count_dict:
				count_dict[junc]=0
			count_dict[junc]+=1
			pc+=1
	print ('Total count:',tc)
	print ('	Paired count:',pc)
	my_list=[[count_dict[junc],junc] for junc in count_dict]
	my_list.sort(reverse=True)
	new=open(outfile,'a')
	for one in my_list:
		new.write(one[1]+'	'+str(one[0])+'	'+str(float(one[0])*100/pc)+'\n')
	new.close()
	


def find_CT_AC(splice5p,splice3p,ingtf,ctac):
	import re
	this_dict1=dict()
	with open(splice5p,'r') as current:
		lc=0
		for line in current:
			lc+=1
			split=line.split(':')
			split[-1]=split[-1].upper()
			tid=split[0].split('::')[0]
			if tid not in this_dict1:
				this_dict1[tid]={5:'',3:''}
			this_dict1[tid][5]=split[-1].strip('\n')

		print ('5 prime splices:',lc)
	with open(splice3p,'r') as current:
		lc=0
		for line in current:
			lc+=1
			split=line.split(':')
			split[-1]=split[-1].upper()
			tid=split[0].split('::')[0]
			if tid not in this_dict1:
				this_dict1[tid]={5:'',3:''}
			this_dict1[tid][3]=split[-1].strip('\n')
		print ('3 prime splices:',lc)

	new=open(ctac,'a')
	gtf =open(ingtf,'r')
	n=-1  #设置参数n，以便在第一次循环时（n=0）可进入，并且用tid1标记此时的编号，以便tid相同时不重复记录transcript行
	for tid in this_dict1:        
		if this_dict1[tid][5][-2:] == 'CT' and this_dict1[tid][3][-2:] == 'AC':
			n+=1        
			if this_dict1[tid][5][-5:-4]=='+':
				a=this_dict1[tid][5].split('-')[0]
				b=str(int(''.join(re.findall(".*-(.*)\(.*",this_dict1[tid][3])))+1)
			if this_dict1[tid][5][-5:-4]=='-':
				a=str(int(''.join(re.findall(".*-(.*)\(.*",this_dict1[tid][5])))+1)
				b=this_dict1[tid][3].split('-')[0] 
			tid=tid[:-3]
			for line in  gtf:
				flag=False             ##设立标志
				transcript_id = line.split()[11].strip('"";')
				if n==0 and line.split()[2]=='transcript' and tid==transcript_id:
					while flag==False:  ##当执行到transcript的最后一个exon跳出while
						next_line1=next(gtf)
						if  line.split()[6] == '+' and next_line1.split()[4] == a:
							new.write(line)
							new.write(next_line1)
							new.write(next(gtf))                            
							flag=True   ##设立标志
						if  line.split()[6] == '-' and next_line1.split()[4] == b:
							new.write(line)
							new.write(next_line1)
							new.write(next(gtf))                            
							flag=True   ##设立标志
				elif n!=0 and line.split()[2]=='transcript' and tid==transcript_id:  ##匹配到编号   
					while flag==False:  ##当执行到transcript的最后一个exon跳出while
						next_line1=next(gtf)
						if  line.split()[6] == '+' and next_line1.split()[4] == a:
							if tid1 != tid:
								new.write(line)
								new.write(next_line1)
								new.write(next(gtf))                            
								flag=True   ##设立标志
							if tid1 == tid:
								new.write(next(gtf))
								flag=True   ##设立标志
						if  line.split()[6] == '-' and next_line1.split()[4] == b:
							if tid1 != tid:
								new.write(line)
								new.write(next_line1)
								new.write(next(gtf))
								flag=True   ##设立标志
							if tid1 == tid:
								new.write(next(gtf))
								flag=True   ##设立标志
				if flag==True:  ##当找到比对的exon后，跳出for line in gtf 循环 
					break
			tid1=tid  #用tid1标记此时的编号
			gtf.seek(0)


	new.close()
	gtf.close()

def simple_statistics(gtf,outfile):
	gtf=open(gtf,'r')
	new=open(outfile,'a')
	transcript_num=0
	exon_num=0
	total_transcript_length=0
	total_exon_length=0
	#统计gtf情况
	for line in gtf:
		if line.split()[2]=='transcript':
			transcript_num+=1
			transcript_length=int(line.split()[4])-int(line.split()[3])
			total_transcript_length=total_transcript_length+transcript_length
		if line.split()[2]=='exon':
			exon_num+=1
			exon_length=int(line.split()[4])-int(line.split()[3])
			total_exon_length=total_exon_length+exon_length
	new.write('Average transcript length:' + str(total_transcript_length/transcript_num)+'\n')
	new.write('CT_AC Average exon number:' + str(exon_num/transcript_num)+'\n')
	new.write('CT_AC Average exon length:' + str(total_exon_length/exon_num)+'\n')

'''def find_full_transcript(a,b,c,d):
	ctac = open(a,'r')
	gtf = open(b,'r')
	new = open(c,'a')
	out = open(d,'a')
	exon_length,total_exon_length,transcript_num=0,0,0
	for line in ctac:
		if line.split()[2]=='transcript':
			transcript_num+=1
		if line.split()[2]=='transcript':
			for line1 in gtf:
				#next_line1=next(gtf)            
				if line1.split()[11]==line.split()[11]:
					new.write(line1)
			gtf.seek(0)
	ctac.close()
	gtf.close()
	new.close()
	out.close()'''
		



from optparse import OptionParser
import os,sys
parser = OptionParser()

parser.add_option("-g", "--gtf_file", dest="gtf_file", help="Absolute path to the input gtf file [required]",type='str',default='')
parser.add_option("-o", "--out_file", dest="out_file", help="Absolute path to the output file [required]",type='str',default='')
parser.add_option("-c", "--ctac_file", dest="ctac_file", help="Absolute path to the CT-AC file [required]",type='str',default='')
###parser.add_option("-t", "--fulltrans_file", dest="fulltrans_file", help="Absolute path to the fulltrans file [required]",type='str',default='')
parser.add_option("-5", "--prime5", dest="prime5", help="Absolute path to the 5 prime base",type='str',default='')
parser.add_option("-3", "--prime3", dest="prime3", help="Absolute path to the 3 prime base",type='str',default='')
(options,args)=parser.parse_args()


if (options.gtf_file=='') or (options.out_file==''):
	print ('Required options not provided; -g and -o and -c are required')
	parser.print_help()
	sys.exit(1)

if options.prime5=='':
	options.prime5=options.out_file

if options.prime3=='':
	options.prime3=options.out_file

print ('Extracting the splice sites from the gtf file...')
extract_splice(options.gtf_file,options.prime3+'_3.bed',options.prime5+'_5.bed',options.out_file)

print ('Getting the splice sequences...')
os.system('bedtools getfasta -fi  "/data3/zenghuanbin/reference/human/GRCh38.p13.genome.fa" -tab -name+ -s -fo '+options.prime3+'_3.nuc -bed '+options.prime3+'_3.bed')
os.system('bedtools getfasta -fi "/data3/zenghuanbin/reference/human/GRCh38.p13.genome.fa"  -tab -name+ -s -fo '+options.prime5+'_5.nuc -bed '+options.prime5+'_5.bed')


print ('Counting the splice dinucleotides per end...')
count_dinucleotide(options.prime3+'_3.nuc',options.prime3+'_3.count')
count_dinucleotide(options.prime5+'_5.nuc',options.prime5+'_5.count')

print ('Counting the paired splice dinucleotides per end, each pair per intron...')
paired_dinucleotides(options.prime5+'_5.nuc',options.prime3+'_3.nuc',options.out_file)

print ('Findind out the CT-AC transcript...')
find_CT_AC(options.prime5+'_5.nuc',options.prime3+'_3.nuc',options.gtf_file,options.ctac_file)

print ('Simple Statistics...')
simple_statistics(options.ctac_file,options.out_file)

print ('full transcript...')
find_full_transcript(options.ctac_file,options.gtf_file,options.fulltrans_file,options.out_file)


sys.exit(0)
