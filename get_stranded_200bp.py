def get_stranded_200bp(infile,outfile):
	import os
	yui=open(infile,'r')
	new=open(outfile,'a')
	length=0
	hold=''
	total=0
	passed=0
	for line in yui:
		split=line.split()
		if line[0]=='#' or len(split)<7 or split[6]=='.':
			continue
		if split[2]=='transcript':
			total+=1
			if length>=200:
				passed+=1
				new.write(hold)
			hold=line
			length=0
		elif split[2]=='exon':
			hold+=line
			length+=abs(int(split[4])-int(split[3]))+1
	yui.close()
	if length>=200:
		passed+=1
		new.write(hold)
##	print os.path.split(infile)[-1],
	print ('Total:',total,)
	print ('Passed:',passed)


from optparse import OptionParser
import os,sys
parser = OptionParser()
parser.add_option('-i','--infile', dest='infile',help='Absolute path to the input gtf file',type='str',default='')
parser.add_option('-o','--outfile', dest='outfile',help='Absolute path to the output gtf file',type='str',default='')
(options,args)=parser.parse_args()
#if (options.infile=='') or (options.outfile==''):
#	print '\nRequired filed(s) not supplied\n# The program extracts stranded GTF files of at least 200bp\n'
#	parser.print_help()
#	sys.exit(1)

get_stranded_200bp(options.infile,options.outfile)
