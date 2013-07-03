#!/usr/bin/python
#I used trim zeros from ends, trim with mean, trim until mean, window method trimming (one side or bothside) after every window trim do trim until mean
#check number of parameters passed
#The feature like of 454hiv is added....select the scores until the score falls below mean quality or skip the homopolyer base with lower score
#Low qual region trimming added
import sys
import os,re

#Popen("clear", shell=True).communicate()[0]
#os.system("clear")
version=1.1

if len(sys.argv)==1 or sys.argv[1]=="-h" or sys.argv[1]=="-help" or sys.argv[1]=="--help":
	sys.stdout.flush()
	print "QTrim: Highly sensitive 454 pyrosequence Quality Trimming tool\n\
QTrim Version: " + str(version) +"\n\n\
Required options:\n\
Input File: -fastq fastqfile OR both -fasta fastafile -qual qualityfile\n\n\
Other options:\n\n\
Output file: -o [Default filename: Outputfile]\n\
Mean quality: -m [INT] Range: 0-40 Default: 20\n\
Minimum read length: -l [INT] Default: 50\n\
Mode: -mode [1,2,3,4] Default: 2 \n\
remove keys: -rk [INT]\n\
Verbose: -verbose\n\
WindowSize: -w [INT] [Default is mininum length]\n\
Output file format: -out_format [Output file format: 1) Fastq file with INT quality score 2) Fastq file with ASCII quality score 3) Sequence and Quality in different Fasta file with Base name provided in output filename]\n\
Sequence statistics in id: -seq_id_stat\n\
Analytical plotting: -plot plot_format (supports these formats:eps, pdf, svg, svgz )\n\n\
Example command:\n\
QTrim_v1_1 -fastq myfastqfile #Runs with all default values\n\
QTrim_v1_1 -verbose -fasta fastafile -qual qualityfile -l 10 -m 30 -o outputfilename -mode 2 -out_format 2\n\
OR\n\
/Fullpath/to/QTrim_v1_1 -fastq fastqfile -l 50 -m 30 -o outputfilename -mode 3 -out_format 3 -seq_id_stat -plot pdf"

	sys.exit(0)

from operator import itemgetter, attrgetter

try:
	from Bio import SeqIO
	from Bio.SeqIO.QualityIO import PairedFastaQualIterator
except ImportError:
	print "Biopython is not installed or is not in the path of python"
	sys.exit(1)




outputfile_dirname="./"
output_file_supplied=False
output_fastq_file='Outputfile';output_fastq,ext=os.path.splitext(output_fastq_file)
qual_input_check=False
qual_file=''
fasta_input_check=False
fasta_file=''
fastq_input_check=False
fastq_file=''
#min_quality=20	#default minimum quality (not used)
mean_quality=20	#default mean quality
min_length=50	#default minimum readlength
user_mode_option=2	#default user choice
user_output_format=2 #default output format
window_size=min_length	#default window size is min length
homopolymer_size=2
format_plot="pdf"


verbose_supplied=False
minqualityusersupply=False
meanqualityusersupply=False
minlengthusersupply=False
user_mode_option_supplied=False
output_format_supplied=False
window_size_supplied=False
homopolymer_size_supplied=False
low_qual_trim_supplied=False
seq_id_stat_supplied=False
plot_supplied=False
remove_key_supplied=False

low_qual_trim_limit_score=20
low_qual_trim_length=6

every_seq_len_original=[]
every_seq_mean_original=[]



output_seq_len_discarded=[]
output_seq_mean_discarded=[]





def get_filetype(input_file):
	try:
		check_seqtype_handle=open(input_file)
	except IOError as Err:
		print "File - ", input_file, " not found"
		sys.exit(1)
	
	first_line=check_seqtype_handle.readline()
	second_line=check_seqtype_handle.readline().strip()
	third_line=check_seqtype_handle.readline().strip()
	forth_line=check_seqtype_handle.readline().strip()
	second_matched=re.search(r'^[A|T|C|G|N]*$', second_line, re.I)
	third_matched=re.search(r'^\+', third_line, re.I)
	
	
	if first_line[0]==">" and second_matched :
		#print "It is definitely fasta"
		readtype='fasta'
		check_seqtype_handle.close()
		return readtype
	elif first_line[0]=="@" and second_matched and third_matched:
		#print "It is definitely fastq"
		readtype='fastq'
		check_seqtype_handle.close()
		return readtype
	else:
		#print "The file", input_file, "is not fasta or fastq format"
		check_seqtype_handle.close()
		print "The file " + input_file + " is not fasta or fastq format"
		sys.exit(1)

def total_good_scores(scores):
	total=0
	for each_score in scores:
		if each_score>20:
			total+=1
	return total

def check_value(value):
	if value.isdigit():
		return True
	else:
		return False

def convert2fastq(fasta_file,qual_file):
	print "\n\nConverting fasta and qual files to fastq file. The fastq file will be saved as filename: Intermediate_Fastq.fastq\n"
	fastq_intermediate=open("Intermediate_Fastq.fastq", "w")
	fastq_records=PairedFastaQualIterator(open(fasta_file), open(qual_file))
	SeqIO.write(fastq_records, fastq_intermediate, "fastq")
	fastq_intermediate.close()
	print "Convertion complete.....\n"

valid_params=['-m ','m ', 'mean_quality', 'meanquality', 'meanqual', '-l', 'l', 'min_length', 'minlength', '-minlen', '-fastq', 'fastq', '-fasta', 'fasta', '-qual', 'qual', '-out_format', 'out_format', '-mode', 'mode', '-w', 'w', '-window', 'window', '-hs', 'hs', '-low_qual_trim', 'low_qual_trim', '-v', 'v', 'verbose', '-V', 'V', '-seq_id_stat', '-plot', '-rk']


def GetOption(syscmd, option_construct):	# syscmd is sys.argv[1:] and option_construct eg: ['fastq=', 'fasta=', 'qual=', 'version',.......]

	#checking if two non options are together
	for non_option in range(2, len(sys.argv)):
		if not sys.argv[non_option-1][0]=='-' and not sys.argv[non_option][0]=='-':
			print "Unknown option", sys.argv[non_option]
			sys.exit(1)
			
	known_options=[]
	default_options=[]
	for each_opt in option_construct:
		if each_opt.endswith("="):
			known_options.append(each_opt[:-1])
		elif each_opt.endswith(":"):
			default_options.append(each_opt[:-1])
		else:
			known_options.append(each_opt)
	#print known_options		
	options_list=[]
	def option_required(option):
		for opt in option_construct:
			if opt[-1]=="=":
				if str(option)+"=" in option_construct:
					return True
			else:
				if str(option) in option_construct:
					return False
	
	def default_option_required(option):
		for opt in option_construct:
			if opt[-1]==":":
				if str(option)+":" in option_construct:
					return None
	
	

	for opt in range(len(syscmd)):
		if syscmd[opt].startswith('-'):
			if syscmd[opt][1:] in known_options:
			
				if option_required(syscmd[opt][1:]):
					try:
						if syscmd[opt+1].startswith("-"):
							print "Option", syscmd[opt], "require argument"
							sys.exit(1)
						else:
							#print syscmd[opt], "added"
							options_list.append((syscmd[opt],syscmd[opt+1]))
					except IndexError:
						if syscmd[opt].startswith("-"):
							print "Option", syscmd[opt], "require argument"
							sys.exit(1)
						else:
							#print syscmd[opt], "added"
							options_list.append((syscmd[opt],syscmd[opt+1]))
					#options_list.append((syscmd[opt],syscmd[opt+1]))
				else:
					try:
						if syscmd[opt+1].startswith("-"):
							#print syscmd[opt],"added"
							options_list.append((syscmd[opt],''))
						else:
							print "Option", syscmd[opt], "do not require argument"
							sys.exit(1)
					except IndexError:
						if syscmd[opt].startswith("-"):
							#print syscmd[opt], "added"
							options_list.append((syscmd[opt],''))
						else:
							print "Option", syscmd[opt], "do not require argument"
							sys.exit(1)

			elif syscmd[opt][1:] in default_options:
				try:
					if syscmd[opt+1].startswith("-"):
						options_list.append((syscmd[opt],''))
					else:
						
						options_list.append((syscmd[opt],syscmd[opt+1]))
				except IndexError:
						options_list.append((syscmd[opt],''))
			else:
				print "Unknown option", syscmd[opt], "in the command"
				sys.exit(1)

	return options_list			

for opt, arg in GetOption(sys.argv[1:],['version','fastq=', 'fasta=', 'qual=', 'output=', 'o=', 'm=', 'meanquality=', 'l=', 'minlen=', 'mode=', 'w=', 'window=', 'out_format=', 'seq_id_stat', 'plot:', 'v', 'verbose', 'rk=']):
	if opt in ('-version'):
		print "QTrim version" , version
		sys.exit(0)
	elif opt in ('-fastq'):
		fastq_file=arg
		fastq_input_check=True
		get_filetype(arg)
		#print "fastq file ", fastq_file
	elif opt in ('-fasta'):
		fasta_file=arg
		fasta_input_check=True
		get_filetype(arg)
		#print "fasta file ", fasta_file
	elif opt in ('-qual'):
		qual_file=arg
		qual_input_check=True
		#print "qual file ", qual_file
	elif opt in ('-o','-output'):
		output_fastq_file=arg
		
		#print output_fastq_file
		outputfile_dirname=os.path.dirname(os.path.abspath(os.path.expanduser(output_fastq_file)))
		if not os.path.exists(outputfile_dirname):
			print "The path to the supplied file directory does not exist. Please check the path of file"
			sys.exit(1)
			
		outputfile_basename=os.path.basename(output_fastq_file)
		output_fastq,ext=os.path.splitext(outputfile_basename)
		output_fastq=outputfile_dirname + "/" + output_fastq
		#print output_fastq
		output_file_supplied=True
		#print "output file ", output_fastq_file
	elif opt in ('-m', '-meanquality'):
		if check_value(arg):
			mean_quality=int(arg)
			meanqualityusersupply=True
			if mean_quality<1:
				print "Mean quality cannot be zero or less"
				sys.exit(1)
		else:
			print "Error in Mean quality input. Error type: Value not supplied or incorrect value.\nPlease check your command.\nNow sys.exiting the program...\n"
			sys.exit(1)
		#print "mean qual ", mean_quality
	elif opt in ('-l', '-minlen'):
		if check_value(arg):
			min_length=int(arg)
			minlengthusersupply=True
			if min_length<1:
				print "Minimum Length cannot be zero or less"
				sys.exit(1)
		else:
			print "Error in Minimum read length input. Error type: Value not supplied or incorrect value.\nPlease check your command.\nNow sys.exiting the program...\n"
			sys.exit(1)
		#print "Min len ", min_length
	elif opt in ('-mode'):
		if check_value(arg):
			user_mode_option=int(arg)
			user_mode_option_supplied=True
			if user_mode_option<1 or user_mode_option > 4:
				print "Mode should be 1 2 3 or 4"
				sys.exit(1)
		else:
			print "Error in mode input. Error type: Value not supplied or incorrect value.\nPlease check your command.\nNow sys.exiting the program...\n"
			sys.exit(1)
		#print "Mode ", user_mode_option
	elif opt in ('-w', '-window'):
		if check_value(arg):
			window_size=int(arg)
			window_size_supplied=True
			if window_size<1:
				print "Zero or negative window size not allowed"
				sys.exit(1)
		else:
			print "Error in Window Size input. Error type: Value not supplied or incorrect value.\nPlease check your command.\nNow sys.exiting the program...\n"
			sys.exit(1)
		#print "Window size ", window_size
	elif opt in ('-hs'):
		if check_value(arg):
			homopolymer_size=int(arg)
			homopolymer_size_supplied=True
			if homopolymer_size<2:
				print "Homopolymer Size less than 2 not allowed. Default value will be used"
				sys.exit(1)
				
		else:
			print "Homopolymer size should be a number"
			sys.exit(1)
		#print "homopolymer size ", homopolymer_size
	elif opt in ('-out_format'):
		user_output_format=int(arg)
		output_format_supplied=True
		#print "Output fomrat ", user_output_format
	elif opt in ('-low_qual_trim'):
		low_qual_trim_limit_score,low_qual_trim_length=arg.split(",")
		low_qual_trim_limit_score=int(low_qual_trim_limit_score)
		low_qual_trim_length=int(low_qual_trim_length)
		low_qual_trim_supplied=True
		#print "low qual trim ", arg
	elif opt in ('-seq_id_stat'):
		seq_id_stat_supplied=True
		#print "seq id stat True"
	elif opt in ('-plot'):
		plot_supplied=True
		if arg!='':
			format_plot=arg
		
		#supported_formats=['emf', 'eps', 'pdf', 'png', 'ps', 'raw', 'rgba', 'svg', 'svgz']
		supported_formats=['eps', 'pdf', 'svg', 'svgz']
		if not format_plot.lower() in supported_formats:
			print "The plot format", arg, "is not supported. Supported formats are", supported_formats
			sys.exit(1)
		try:
			import numpy as np
		except ImportError:
			print "Numpy is not installed or is not in the path of python"
			sys.exit(1)

		try:
			import matplotlib
	
		except ImportError:
			print "Matplotlib is not installed or is not in the path of python. Install numpy, before installing matplotlib. Simple Version can be used if no Matplot is installed."
			print "If Matplotlib is not installed, try running without the option -plot" 
			sys.exit(1)
		#print "Plot True"
	elif opt in ('-rk'):
		remove_key_supplied=True
		key_length=int(arg)
	elif opt in ('-v', '-verbose'):
		verbose_supplied=True
		#print "verbose True"
	else:
		print "Unknown option given:", opt
		sys.exit(1)
		

if verbose_supplied==False:
	print "Cleaning the sequence reads...please be patient"

def Verbose(toprint):
	if verbose_supplied==True:
		print "\r"+toprint

if fasta_input_check==True and fastq_input_check==True:
	#print "Fastq and Fasta files detected. Fastq will be processed\n"
	Verbose("Fastq and Fasta files detected. Fastq will be processed\n")
elif fastq_input_check==True and fasta_input_check==False and qual_input_check==False:
	#print "Fastq file will be processed\n"
	Verbose("Fastq file will be processed\n")
elif fastq_input_check==False and fasta_input_check==True and qual_input_check==True:
	#print "fasta file and qual files will be processed"
	Verbose("fasta file and qual files will be processed")
	convert2fastq(fasta_file,qual_file)
	fastq_file="Intermediate_Fastq.fastq"
elif fastq_input_check==False and fasta_input_check==True and qual_input_check==False:
	#print "Qual file not supplied for fasta file input:"
	Verbose("Qual file not supplied for fasta file input:")
	#print "Program sys.exiting now... "
	Verbose("Program sys.exiting now...")
	sys.exit(1)
elif fastq_input_check==False and fasta_input_check==False and qual_input_check==False:
	#print "No input file detected. sys.exiting the program now..."
	Verbose("No input file detected. sys.exiting the program now...")
	sys.exit(1)
else:
	#print "File input wrong. Check command line supply of fastq or fasta or qual file. sys.exiting the program now..."
	Verbose("File input wrong. Check command line supply of fastq or fasta or qual file. sys.exiting the program now...")
	sys.exit(1)

#print "Your inputs are:\n"
Verbose("Your inputs are:\n")

if meanqualityusersupply==True:
	#print "Mean quality score: ", mean_quality
	Verbose("Mean quality score: "+str(mean_quality))
else:
	#print "Mean quality score: "+str(mean_quality) +" [By Default]"
	Verbose("Mean quality score: "+str(mean_quality) +" [By Default]")

if minlengthusersupply==True:
	#print "Minimum read length: ", min_length
	Verbose("Minimum read length: "+ str(min_length))
else:
	#print "Minimum read length: "+str(min_length)+" [By Default]"
	Verbose("Minimum read length: "+str(min_length)+" [By Default]")

if user_mode_option_supplied==True:
	#print "Mode: ", user_mode_option
	Verbose("Mode: "+str(user_mode_option))
else:
	#print "Mode: "+str(user_mode_option)+ "[By Default]"
	Verbose("Mode: "+str(user_mode_option)+ "[By Default]")

if low_qual_trim_supplied==True:
	#print "Low quality region trimming: Limit Score = "+str(low_qual_trim_limit_score)+ " Length = "+str(low_qual_trim_length)
	Verbose("Low quality region trimming: Limit Score = "+str(low_qual_trim_limit_score)+ " Length = "+str(low_qual_trim_length))

if window_size_supplied==True:
	#print "Window size: ", window_size
	Verbose("Window size: "+str(window_size))
else:
	#print "Window size: ", window_size, "[By Default]"
	Verbose("Window size: "+ str(window_size) + "[By Default]")
if user_mode_option==6 and homopolymer_size_supplied==True:
	#print "Homopolymer Size: ", homopolymer_size
	Verbose("Homopolymer Size: "+ str(homopolymer_size))
else:
	#print "Homopolymer Size: ", str(2), "[By Default]"
	Verbose("Homopolymer Size: " +str(2)+"[By Default]\t[Only for mode 5 and 6]")
	
if output_format_supplied==True:
	if user_output_format==1:
		#print "Output File Format: ", user_output_format, "[hivhiv will be FASTQ with quality scores in integer value]"
		Verbose("Output File Format: "+str(user_output_format)+ "[Your output file format will be FASTQ with quality scores in integer value]")
	elif user_output_format==2:
		#print "Output File Format: ", user_output_format, "[Your output file format will be FASTQ with quality scores in ASCII character]"
		Verbose("Output File Format: "+str(user_output_format)+ "[Your output file format will be FASTQ with quality scores in ASCII character]")
	elif user_output_format==3:
		#print "Output File Format: ", user_output_format, "[Your output file format will be FASTA. Sequence and quality scores (in integer value) will be in two different files]"
		Verbose("Output File Format: "+str(user_output_format)+"[Your output file format will be FASTA. Sequence and quality scores (in integer value) will be in two different files]")





#print "Your Inputs: \nMinimum Quality score: "+ sys.argv[2] + "\nMean Quality Score: " + sys.argv[4] + "\nMinimum Read Length; " + sys.argv[3] + "\n454 Reads Cleaning in progress....."

#print "You have following options for processing: \n \
#1) Trim 3' end and remove Ns in middle of reads\n \
#2) Trim 3' end low quality scores\n \
#3) Trim 5' and 3' end and remove any Ns in the middle of reads\n \
#4) Trim 5' and 3' low quality scores\n \
#0) sys.exit the program without any data cleaning\n \
#Users are adviced to use all options one at a time to get the best result\n \
#Your Option Please [1,2,3,4,0]: ",

#user_mode_option=input()

try:
	file_handle=open(fastq_file)
except IOError:
	print "The file \""+ fastq_file + "\" cannot be found. Please check the path of the file."
	sys.exit(1)

if user_output_format==1 or user_output_format==2:
	outfile_handle=open(output_fastq+".fastq", "w")
if user_output_format==3:
	outfile_seq_handle=open(output_fastq+".fasta", 'w')
	outfile_qual_handle=open(output_fastq+".qual", 'w')
	#discarded_seq_handle=open(output_fastq+"_discarded_reads.fasta","w")
	#discarded_qual_handle=open(output_fastq+"_discarded_reads.qual","w")

if user_mode_option!=5 and user_mode_option!=6:
	discarded_handle=open(output_fastq+"_discarded_reads","w")
	read_with_zeros_handle=open(output_fastq+"_Reads_with_zeros_in_middle", "w")
	position_of_zero_handle=open(output_fastq+"_ZeroPositions", "w")
	trimming_handle=open(output_fastq+"_Trimming_report", "w")

stat_handle=open(output_fastq+"_stat.txt", "w")


if user_mode_option==2 or user_mode_option==4:
	position_of_zero_handle.write("Positions of zeros in the middle of read will not be published for option 2 or 4. Mean quality is calculated without trimming them off.\n")

#Processing your input

while(True):
	if user_mode_option in [0,1,2,3,4,5,6]:
		break
	else:
		print "Please input valid option [0,1,2,3,4,5,6]: "
		user_mode_option=input()

while(True):
	if user_output_format in [1,2,3]:
		break
	else:
		print "Please input valid output format [1,2,3]: "
		user_output_format=input()


if user_mode_option==0:
	#print "\nsys.exiting the program....\n"
	Verbose("\nsys.exiting the program....\n")
	sys.exit(0)

if user_mode_option==1 or user_mode_option==2:
	trimming_handle.write("The number to the right to identifier is the number of low quality bases removed including zeros\n")
if user_mode_option==3 or user_mode_option==4:
	trimming_handle.write("The number to the left is bases removed at 5'end and the number to the right is bases removed at 3' end including zeros\n")

from datetime import datetime
#print "Program started at :", datetime.now()
Verbose("Program started at :"+str(datetime.now()))



def get_readlen_info(lengths):	#function to classify the read lengths with integral multiplicaiton of min length
	max_len=max(lengths)
	minlen=min_length
	bins=[]
	if (max_len%minlen)==0:
		totalbins=max_len/minlen
	else:
		totalbins=(max_len/minlen)+1


	for eachbin in range(totalbins+1):
		bins.append(0)

	for eachlen in lengths:
		#print eachlen/minlen
		bins[eachlen/minlen]+=1
		#if eachlen/minlen==0:
		#	bins[eachlen/minlen]+=1
		#else:
		#	bins[(eachlen/minlen)]+=1

	#print "Output Reads Length range and count info"
	Verbose("Output Reads Length range and count info")
	for eachbin in range(len(bins)):
		#print str(eachbin*minlen)+" - "+str((eachbin+1)*minlen-1)+"    \t: "+ str(bins[eachbin])+"\n"
		Verbose(str(eachbin*minlen)+" - "+str((eachbin+1)*minlen-1)+"    \t: "+ str(bins[eachbin]))



#print "Counting the number of reads in your file..."
Verbose("Counting the number of reads in your input file...")
total_count=0
inputfile_total_seqs=0
origin_read_lengths=[]

for rec in SeqIO.parse(file_handle, "fastq"):
	total_count+=1
	inputfile_total_seqs+=1
	if remove_key_supplied==True:
		origin_read_lengths.append(len(str(rec.seq)[key_length:]))
	else:
		origin_read_lengths.append(len(str(rec.seq)))

file_handle.close()	#file closing after counting the number of reads

file_handle=open(fastq_file)	#opening after for cleaning
#print "Total reads in your file: ", total_count
Verbose("Total reads in your file: "+str(total_count))

get_readlen_info(origin_read_lengths)
#print "Please Look at your read length distribution from your input file while quality trimming" 
Verbose("Please Look at your read length distribution from your input file while quality trimming")

#print "\n\nQuality Trimming of the reads started at :", datetime.now()
Verbose("\n\nQuality Trimming of the reads started at :"+str(datetime.now()))

#print "Now Quality Trimming and Quality filtering/cleaning fastq file "+fastq_file+" .....\n"
Verbose("Now Quality Trimming and Quality filtering/cleaning fastq file "+fastq_file+" .....\n")

identifier=''
nt_seq=''
qual_identifier=''

total_reads_output=0
zero_in_middle=0
output_lengths=[]
output_means=[]
trimmed_at_5=[]
trimmed_at_3=[]

def convert_to_ascii_char(int_val):
	ascii_val=""
	scores=int_val.split()
	for each_score in scores:
		ascii_val+=chr(int(each_score)+33)	

	return ascii_val

def convert_int_to_ascii_char(scores):
	ascii_val=""
	
	for each_score in scores:
		ascii_val+=chr(int(each_score)+33)	

	return ascii_val
def convert_fastq_int_2_fastq_ascii(output_fastq):
	try:
		file_handle=open(output_fastq)
	except:
		file_handle=open(output_fastq + ".fastq")
		
	out_handle=open("Intermediate_Fastq.fastq",'w')
	for readline in file_handle:
		#print readline
		if readline.strip()=="":
			continue	
		elif readline[0]=='@':
			out_handle.write(readline.strip()+"\n")
		
		elif readline.strip().isalpha():
			out_handle.write(readline.strip()+"\n")
		
		elif readline[0]=='+':
			out_handle.write(readline.strip()+"\n")
		
		elif readline.strip().replace(" ","").isdigit():
			out_handle.write(convert_to_ascii_char(readline.strip())+"\n")
	

def convert_fastq_int_2_fasta_qual(output_fastq):
	import os
	file_in_handle=open(output_fastq+".fastq")
	fasta_out_handle=open(output_fastq+".fasta", "w")
	qual_out_handle=open(output_fastq+".qual", "w")

	identity=""
	for record in file_in_handle:
		if record[0]=="@":
			identity=record[1:]
			fasta_out_handle.write(">"+identity)
			qual_out_handle.write(">"+identity)		
		if record.strip().isalpha():
			fasta_out_handle.write(record)		
		if record.strip().replace(" ","").isdigit():
			qual_out_handle.write(record)
	file_in_handle.close()
	fasta_out_handle.close()
	qual_out_handle.close()
	os.system("rm \""+output_fastq+".fastq\"")
	
	
def convert_2_fastq_ascii1(fasta_file):
	import os
	fastq_output=open(output_fastq+"_ascii.fastq",'w')
	#fasta_input=open(fasta_file+".fasta")
	#qual_input=open(fasta_file+".qual")

	#records=PairedFastaQualIterator(open(sys.argv[1]), open(sys.argv[2]))
	#count=SeqIO.write(records,fastq_output, "fastq") #The last option "fastq" can be replaced with "fastq-sanger" for sanger format output, "fastq-illumina" for illumina format
	#records=PairedFastaQualIterator(fasta_input, qual_input)
	#SeqIO.write(records, fastq_output, "fastq")

	file_handle=open(output_fastq+".fastq")

	for readline in file_handle:
		#print readline
		if readline.strip()=="":
			continue	
		elif readline[0]=='@':
			fastq_output.write(readline.strip()+"\n")
		
		elif readline.strip().isalpha():
			fastq_output.write(readline.strip()+"\n")
		
		elif readline[0]=='+':
			fastq_output.write(readline.strip()+"\n")
		
		elif readline.strip().replace(" ","").isdigit():
			fastq_output.write(convert_to_ascii_char(readline.strip())+"\n")
		else:
			fastq_output.write("\n")
	
	fastq_output.close()
	#fasta_input.close()
	#qual_input.close()
	file_handle.close()
	#os.system("rm \""+output_fastq+".fasta \""+output_fastq+".qual\"")
	os.system("rm \""+output_fastq+".fastq\"")
	os.system("mv \""+output_fastq+"_ascii.fastq\"" + " " + "\""+output_fastq+ ".fastq\"")


	
def convert_2_fastq_ascii(fasta_file):
	import os
	fastq_output=open(output_fastq+".fastq","w")
	fasta_input=open(fasta_file+".fasta")
	qual_input=open(fasta_file+".qual")

	#records=PairedFastaQualIterator(open(sys.argv[1]), open(sys.argv[2]))
	#count=SeqIO.write(records,fastq_output, "fastq") #The last option "fastq" can be replaced with "fastq-sanger" for sanger format output, "fastq-illumina" for illumina format
	records=PairedFastaQualIterator(fasta_input, qual_input)
	SeqIO.write(records, fastq_output, "fastq")
	
	fastq_output.close()
	fasta_input.close()
	qual_input.close()
	os.system("rm \""+output_fastq+".fasta\"" +" " + "\"" + output_fastq + ".qual\"")


def find_mean(array_passed):
	if len(array_passed)==0:
		return 0
	else:
		return sum(array_passed)/len(array_passed)
		

def good_bad_scores(scores):

	bad=0
	good=0
	for each_score in scores:
		if each_score>20:
			good+=1
		else:
			bad+=1
	return(bad>0)




def find_good_score_position(scores, start):
	count=start
	for score_position in range(count,len(scores)):
		if scores[score_position]>=mean_quality:
			count+=1
		else:
			break
	return(count)




def find_good_score_position_keep_homopolymer(nt_seq, scores, start):
	count=start
	for score_position in range(count,len(scores)):
		if scores[score_position]>=mean_quality:
			count+=1
		else:
			if nt_seq[count-homopolymer_size+1:count+1].upper()==homopolymer_size*"A" or nt_seq[count-homopolymer_size+1:count+1]==homopolymer_size*"C" or nt_seq[count-homopolymer_size+1:count+1]==homopolymer_size*"G" or nt_seq[count-homopolymer_size+1:count+1]==homopolymer_size*"T":
				count+=1
			else:
				break
	return(count)

def split_good_seq(seq_id, nt_seq, scores):
	good_seq_list=[]
	good_seq_fragment=[]
	start=0
	while(start!=len(scores)):
		if scores[start]>=mean_quality:
			if user_mode_option==5:
				end=find_good_score_position(scores,start)
			if user_mode_option==6:
				end=find_good_score_position_keep_homopolymer(nt_seq, scores, start)

			good_seq=nt_seq[start:end]
			good_scores=scores[start:end]

			if len(good_seq)>=min_length:
				good_seq_fragment.append(good_seq)
				good_seq_fragment.append(good_scores)
				if len(good_seq)!=len(good_scores):
					print "Sequence length and score length not equal"
					sys.exit(1)
				good_seq_list.append(good_seq_fragment)
				good_seq_fragment=[]

			start=end
		else:
			start+=1

	return good_seq_list








#function to print the read length info with classification of integral multiplication of min length

#function to find the positions of zero
def get_zero_positions(read_id,sequence,qual_id,check_array_for_zero):

	if len(sequence)!=len(check_array_for_zero):
		#print "Sequence length and quality score length are not same. Something is wrong with identifier: "+ read_id
		Verbose("Sequence length and quality score length are not same. Something is wrong with identifier: "+ read_id)
		sys.exit(1)
	zero_positions=[]
	for search_zero in range(len(check_array_for_zero)-1):
		if check_array_for_zero[search_zero]==0:
			zero_positions.append(search_zero)
	if len(zero_positions)!=0:
		#write the sequence with zero in the middle to a file
		read_with_zeros_handle.write(read_id+"\n")
		read_with_zeros_handle.write(sequence+"\n")
		read_with_zeros_handle.write(qual_id+"\n")
		read_with_zeros_handle.write(convert_intqual_to_stringqual(check_array_for_zero)+"\n")
		position_of_zero_handle.write(read_id+"\n"+convert_intqual_to_stringqual(zero_positions)+"\n")
		position_addition=0
		#position_of_zero_handle.write(read_id+": "+str(zero_positions))
		for position in zero_positions:
			del check_array_for_zero[position+position_addition]
			sequence=sequence[0:position+position_addition]+sequence[position+position_addition+1:]
			position_addition-=1
	if len(sequence)!=len(check_array_for_zero):
		#print "Sequence length and quality score length after removing zeros are not same. Something is wrong after removing zero in id: "+ read_id
		Verbose("Sequence length and quality score length after removing zeros are not same. Something is wrong after removing zero in id: "+ read_id)
		sys.exit(1)
	return (sequence, check_array_for_zero)

def check_low_quality_region(nt_seq, scores, first_removal, last_removal):
	count=1
	clean_point=0
	myscore_pos=0
	total_low_scores=0
	while(myscore_pos<len(scores)):
		#print myscore_pos
		if scores[myscore_pos] < low_qual_trim_limit_score:
			count+=1
			total_low_scores+=1
			if count>=low_qual_trim_length and myscore_pos-total_low_scores<min_length: #or total_low_scores==1/100*len(scores)
				#if myscore_pos-low_qual_trim_length>len(scores)-myscore_pos:
				if total_good_scores(scores[:myscore_pos])>total_good_scores(scores[myscore_pos:]):
					last_removal+=len(scores[myscore_pos:])
					nt_seq=nt_seq[:myscore_pos]
					scores=scores[:myscore_pos]
				else:
					first_removal+=len(scores[:myscore_pos])
					nt_seq=nt_seq[myscore_pos+1:]
					scores=scores[myscore_pos+1:]

				myscore_pos=0
				nt_seq, scores,first_removal=trim_until_mean_at_left(nt_seq, scores,first_removal)
				nt_seq, scores,last_removal=trim_until_mean(nt_seq, scores,last_removal)

		else:
			count=0
			#clean_point=myscore_pos

		myscore_pos+=1
	return (nt_seq, scores,first_removal, last_removal)


def convert_intqual_to_stringqual(toConvert):
	if len(toConvert)==0:
		return ""
	elif len(toConvert)==1:
		return str(toConvert[0])
	else:
		str_qual=str(toConvert[0])
		for int_qual in toConvert[1:]:
			str_qual=str_qual+" "+ str(int_qual)
		return str_qual

def trim_until_mean_at_left(nt_seq, scores, first_removal):	#trims at the left end until the high qual is found
	if len(scores)<=50:
		return(nt_seq, scores,first_removal)

	while(scores[0]<30):
		del scores[0]
		nt_seq=nt_seq[1:]
		first_removal+1

	return(nt_seq, scores,first_removal)

def trim_until_mean(nt_seq, scores, last_removal):	#trims at the end until the mean of array is ok
	last_removal_low_quality=0
	if len(scores)<=min_length:
		return(nt_seq, scores, last_removal)

	while(scores[len(scores)-1]<mean_quality):
		del scores[len(scores)-1]
		last_removal_low_quality+=1
		#total_length=len(scores)
		mean_in_array=find_mean(scores)
	nt_seq=nt_seq[0:len(scores)]
	last_removal+=last_removal_low_quality
	return(nt_seq, scores, last_removal)
	



def check_mean_and_trim_3(nt_seq_half, scores_half, last_removal):
	if len(scores)<=min_length:
		return(nt_seq_half, scores_half, last_removal)

	last_removal_low_quality=0
	while(True):
		if find_mean(scores_half)<mean_quality and len(scores_half)>2:
			
			del scores_half[len(scores_half)-1]
			last_removal_low_quality+=1
	
			mean_in_array=find_mean(scores_half)
			
		else:
			break

	nt_seq_half=nt_seq_half[0:len(scores_half)]
	last_removal+=last_removal_low_quality

	#(nt_seq_half,scores_half,last_removal)=trim_until_mean(nt_seq_half, scores_half, last_removal)	#trim scores less than mean qual again from 3' end

	return(nt_seq_half, scores_half, last_removal)


		
def window_method_trimming(nt_seq, scores,last_removal):
	if len(scores)<=min_length:
		return(nt_seq, scores, last_removal)

	while(True):
		if find_mean(scores[-window_size:])<mean_quality and len(scores)>window_size:
			del scores[-1:]
			nt_seq=nt_seq[:-1]
			last_removal+=1
			#(nt_seq,scores,last_removal)=trim_until_mean(nt_seq,scores,last_removal)
		else:
			break
	(nt_seq,scores,last_removal)=trim_until_mean(nt_seq,scores,last_removal)
	if len(nt_seq)!=len(scores):		
		print "len of seq and scores not same in remove_from end with mean"
		sys.exit(1)

	return(nt_seq, scores,last_removal)


def window_method_trimming_bothside(nt_seq, scores,last_removal, first_removal):	#unused function
	if len(scores)<=min_length:
		return(nt_seq, scores, last_removal)

	while(True):
		if find_mean(scores[-window_size:])<mean_quality and len(scores)> window_size:
			del scores[-1:]
			nt_seq=nt_seq[:-1]
			last_removal+=1
			(nt_seq,scores,last_removal)=trim_until_mean(nt_seq,scores,last_removal)
		else:
			break
	while(True):
		if find_mean(scores[:window_size+1])<mean_quality and len(scores)>window_size+1:
			del scores[-1:]
			nt_seq=nt_seq[:-1]
			first_removal+=1
			(nt_seq,scores,first_removal)=trim_until_mean(nt_seq,scores,first_removal)
		else:
			break
	return(nt_seq,scores,last_removal, first_removal)

def remove_inner_lowqual_3(nt_seq,scores,last_removal):
	if len(scores)<=min_length:
		return(nt_seq, scores, last_removal)

	#readlen=len(scores)
	'''if readlen>2*min_length:
		half_3=int(readlen/2)
	elif readlen>min_length and readlen<=2*min_length:
		half_3=min_length
	else:
		return (nt_seq, scores,last_removal)

	#three_end_seq=nt_seq[half_3:]
	#three_end_scores=scores[half_3:]

	#(three_end_seq, three_end_scores, last_removal)=window_method_trimming(three_end_seq, three_end_scores, last_removal)
	(nt_seq,scores,last_removal)=trim_until_mean(nt_seq,scores,last_removal)'''

	(nt_seq, scores, last_removal)=window_method_trimming(nt_seq, scores, last_removal)
	
	if low_qual_trim_supplied==True:
		nt_seq, scores, nowhere, last_removal=check_low_quality_region(nt_seq, scores, 0 , last_removal)
	
	if len(nt_seq)!=len(scores):
		#print len(nt_seq), len(scores)
		#print "len of seq and scores not same in remove_from end with mean"
		Verbose("len of seq and scores not same in remove_from end with mean")
		sys.exit(1)
	#if len(scores)>2*min_length
	return (nt_seq, scores,last_removal)


def remove_inner_lowqual_53(identifier, nt_seq, scores, first_removal, last_removal):
	if len(scores)<=min_length:
		return (nt_seq, scores, first_removal, last_removal)
	
	#(nt_seq, scores, last_removal, first_removal)=window_method_trimming_bothside(nt_seq, scores, last_removal, first_removal)
	#(nt_seq,scores,last_removal)=trim_until_mean(nt_seq,scores,last_removal)
	
	(nt_seq, scores, last_removal)=window_method_trimming(nt_seq, scores, last_removal)
	nt_seq=nt_seq[::-1]
	scores.reverse()

	#(nt_seq,scores,last_removal)=trim_until_mean(nt_seq,scores,last_removal)
	(nt_seq, scores, first_removal)=window_method_trimming(nt_seq, scores, first_removal)
	nt_seq=nt_seq[::-1]
	scores.reverse()

	#(nt_seq,scores,last_removal)=trim_until_mean(nt_seq,scores,last_removal)

	if low_qual_trim_supplied==True:
		nt_seq, scores, first_removal, last_removal=check_low_quality_region(nt_seq, scores, first_removal, last_removal)

	if len(nt_seq)!=len(scores):		
		print "len of seq and scores not same in remove_from end with mean"
		sys.exit(1)

	return (nt_seq, scores, first_removal, last_removal)


def three_prime_low_quality_trim(identifier, nt_seq, qual_identifier, scores,last_removal):
	if len(scores)<min_length:
		return(last_removal, nt_seq, scores)

	last_removal_low_quality=0
	mean_in_array=find_mean(scores)
	total_length=len(scores)
	while(True):
		
		if len(scores)==0:
			break
		elif mean_in_array<mean_quality:
			del scores[len(scores)-1]
			last_removal_low_quality+=1
			total_length=len(scores)
			mean_in_array=find_mean(scores)
			continue
		#elif scores[len(scores)-1]<mean_quality:
		#	del scores[len(scores)-1]
		#	last_removal_low_quality+=1
		#	total_length=len(scores)
		#	mean_in_array=find_mean(scores)
		else:
			break
			
	nt_seq=nt_seq[0:len(scores)]

	last_removal+=last_removal_low_quality
	(nt_seq,scores,last_removal)=trim_until_mean(nt_seq,scores,last_removal)
	(nt_seq,scores,last_removal)=remove_inner_lowqual_3(nt_seq,scores,last_removal)
	if len(nt_seq)!=len(scores):
		print "After trimming low qual from function remove_inner_lowqual_3, seq len and qual len are not same"
		sys.exit(1)
	if len(nt_seq)<min_length:
		return(last_removal, nt_seq, scores)	

	#(nt_seq,scores,last_removal)=trim_until_mean(nt_seq,scores,last_removal)
		
	if len(nt_seq)!=len(scores):
		print "After returning from function trim_until_mean, length of sequence and quality scores are not same"
		sys.exit(1)
	
	return(last_removal, nt_seq, scores)


def any_below_meanqual(quals):
	returnval=False
	for x in quals:
		if x in range(1,mean_quality):
			returnval=True
			break
	return returnval
	
def five_three_prime_low_quality_trim(identifier, nt_seq,scores,first_removal, last_removal,seq_length):

	if len(scores)<min_length:
		return(first_removal, last_removal, nt_seq, scores)

	last_removal_low_quality=0
	first_removal_low_quality=0
	mean_scores=find_mean(scores)
	total_length=len(scores)
	#print "before", len(nt_seq), len(scores)
	while(mean_scores<mean_quality):
		if len(scores)==0:
			break
	#	elif scores[len(scores)-1]<mean_quality or scores[len(scores)-1]<=scores[0]:#compare first and last scores and remove lower one
	#		del scores[len(scores)-1]
	#		last_removal_low_quality+=1
	#		mean_scores=find_mean(scores)
	#		total_length=len(scores)
	#		if len(scores)<min_length:
	#			break
	#	elif scores[len(scores)-1]>scores[0] or scores[0]<mean_quality:#remove first one if lower than last one
	#		del scores[0]
	#		first_removal_low_quality+=1
	#		mean_scores=find_mean(scores)
	#		total_length=len(scores)
	#		if len(scores)<min_length:
	#			break
		else:
			del scores[len(scores)-1]
			last_removal_low_quality+=1
			mean_scores=find_mean(scores)
			total_length=len(scores)
			#if len(scores)<min_length:
			#	break			
	
	
	nt_seq=nt_seq[first_removal_low_quality:seq_length-last_removal_low_quality]

	if len(nt_seq)!=len(scores):
		#print "seq len and qual length not same after trimming low quals from 5 and 3.Point A" 
		Verbose("Sequence length and quality length are not equal")
		sys.exit(1)
	
	last_removal+=last_removal_low_quality
	first_removal+=first_removal_low_quality

	if len(scores)<min_length:
		return(first_removal, last_removal, nt_seq, scores)
	

		
	
	(nt_seq, scores, first_removal, last_removal)=remove_inner_lowqual_53(identifier, nt_seq, scores, first_removal, last_removal)
	#(nt_seq,scores,last_removal)=trim_until_mean(nt_seq,scores,last_removal)

	base_remove_at_ends_handle.write(str(first_removal_low_quality)+"\t"+str(last_removal_low_quality))
	first_removal+=first_removal_low_quality
	last_removal+=last_removal_low_quality
	if len(nt_seq)!=len(scores):
		#print "seq len and qual length not same after trimming low quals from 5 and 3. Point B" 
		Verbose("Sequence length and quality length are not equal")
		sys.exit(1)

	return(first_removal, last_removal, nt_seq, scores)

def user_output_format_1(nt_seq_id, nt_seq, qual_id, scores):
	outfile_handle.write(nt_seq_id)
	outfile_handle.write(nt_seq)
	outfile_handle.write(qual_id)
	outfile_handle.write(scores)

def user_output_format_2 (nt_seq_id, nt_seq, qual_id, scores):
	outfile_handle.write(nt_seq_id)
	outfile_handle.write(nt_seq)
	outfile_handle.write(qual_id)
	outfile_handle.write(scores)

def user_output_format_3(nt_seq_id, nt_seq, qual_id, scores):
	outfile_seq_handle.write(nt_seq_id)
	outfile_seq_handle.write(nt_seq)
	outfile_qual_handle.write(qual_id)
	outfile_qual_handle.write(scores)

def decision_making(identifier, nt_seq, qual_identifier, scores, first_removal, last_removal):

	
	total_length=len(scores)
	discarded=True
	this_score_mean=find_mean(scores)
	
	
	
	if (find_mean(scores)>=mean_quality and total_length>=min_length):	
		
		if user_output_format==1:
			if seq_id_stat_supplied==True:
				user_output_format_1(identifier+"\tMean: "+str(find_mean(scores))+"\tLen: "+str(len(scores))+"\tTrimmed at 5': " +str(first_removal)+ "\tTrimmed at 3': "+str(last_removal)+ "\n", nt_seq+"\n", "+\n", convert_intqual_to_stringqual(scores)+"\n")
			else:
				user_output_format_1(identifier+"\n", nt_seq+"\n", "+\n", convert_intqual_to_stringqual(scores)+"\n")

			#outfile_handle.write(identifier+"\tMean: "+str(this_score_mean)+"\tLen: "+str(len(scores))+"\tTrimmed at 5': " +str(first_removal)+ "\tTrimmed at 3': "+str(last_removal)+ "\n")
			#outfile_handle.write(nt_seq+"\n")
			#outfile_handle.write("+\n")
			#outfile_handle.write(convert_intqual_to_stringqual(scores)+"\n")

		if user_output_format==2:
			if seq_id_stat_supplied==True:
				user_output_format_2(identifier+"\tMean: "+str(find_mean(scores))+"\tLen: "+str(len(scores))+"\tTrimmed at 5': " +str(first_removal)+ "\tTrimmed at 3': "+str(last_removal)+ "\n", nt_seq+"\n", "+\n", convert_int_to_ascii_char(scores)+"\n")
			else:
				user_output_format_2(identifier+"\n", nt_seq+"\n", "+\n", convert_int_to_ascii_char(scores)+"\n")

			#outfile_handle.write(identifier+"\tMean: "+str(this_score_mean)+"\tLen: "+str(len(scores))+"\tTrimmed at 5': " +str(first_removal)+ "\tTrimmed at 3': "+str(last_removal)+ "\n")
			#outfile_handle.write(nt_seq+"\n")
			#outfile_handle.write("+\n")
			#outfile_handle.write(convert_int_to_ascii_char(scores)+"\n")

		if user_output_format==3:
			if seq_id_stat_supplied==True:
				user_output_format_3(">"+identifier[1:]+"\tMean: "+str(find_mean(scores))+"\tLen: "+str(len(scores))+"\tTrimmed at 5': " +str(first_removal)+ "\tTrimmed at 3': "+str(last_removal)+ "\n", nt_seq+"\n", ">"+identifier[1:]+"\n", convert_intqual_to_stringqual(scores)+"\n")
			else:
				user_output_format_3(">"+identifier[1:] + "\n", nt_seq+"\n", ">"+identifier[1:]+"\n", convert_intqual_to_stringqual(scores)+"\n")				

			#outfile_seq_handle.write(">"+identifier[1:]+"\tMean: "+str(this_score_mean)+"\tLen: "+str(len(scores))+"\tTrimmed at 5': " +str(first_removal)+ "\tTrimmed at 3': "+str(last_removal)+ "\n")
			#outfile_seq_handle.write(nt_seq+"\n")
			#outfile_qual_handle.write(">"+identifier[1:]+"\n")
			#outfile_qual_handle.write(convert_intqual_to_stringqual(scores)+"\n")

		
		discarded=False
		output_lengths.append(total_length)
		output_means.append(this_score_mean)
	
		
	else:
		if user_output_format==1:
			discarded_handle.write(identifier+"\tMean: "+"\t\tLength: "+str(len(scores))+"\n")
			discarded_handle.write(nt_seq +"\n")
			discarded_handle.write("+\n")
			discarded_handle.write(convert_intqual_to_stringqual(scores)+"\n")

			#output_seq_len_discarded.append(total_length)
			#output_seq_mean_discarded.append(this_score_mean)

		if user_output_format==2:
	
			discarded_handle.write(identifier+"\tMean: "+"\t\tLength: "+str(len(scores))+"\n")
			discarded_handle.write(nt_seq+"\n")
			discarded_handle.write(qual_identifier+"\n")
			discarded_handle.write(convert_int_to_ascii_char(scores)+"\n")

			#output_seq_len_discarded.append(total_length)
			#output_seq_mean_discarded.append(this_score_mean)

		if user_output_format==3:
			discarded_seq_handle.write(">"+identifier[1:]+"\tMean: "+"\tLength: "+str(len(nt_seq))+"\n")
			discarded_seq_handle.write(nt_seq  +"\n")
			discarded_qual_handle.write(">"+identifier[1:]+"\n")
			discarded_qual_handle.write(convert_intqual_to_stringqual(scores)+"\n")

			#output_seq_len_discarded.append(total_length)
			#output_seq_mean_discarded.append(this_score_mean)

	return (discarded)

#######################################################################################################################
#Main Program starts from here


length_dict={}

percent_count=0
inputfile_max_length=0
for record in SeqIO.parse(file_handle,'fastq'):
	percent_count+=1
	if percent_count==int(total_count*10/100):
		print "10% of reads cleaning completed..."
		#Verbose("10% of reads cleaning completed...")
	elif percent_count==int(total_count*20/100):
		print "20% of reads cleaning completed..."
		#Verbose("20% of reads cleaning completed...")
	elif percent_count==int(total_count*30/100):
		print "30% of reads cleaning completed..."
		#Verbose("30% of reads cleaning completed...")
	elif percent_count==int(total_count*40/100):
		print "40% of reads cleaning completed..."
		#Verbose("40% of reads cleaning completed...")
	elif percent_count==int(total_count*50/100):
		print "50% of reads cleaning completed..."
		#Verbose("50% of reads cleaning completed...")
	elif percent_count==int(total_count*60/100):
		print "60% of reads cleaning completed..."
		#Verbose("60% of reads cleaning completed...")
	elif percent_count==int(total_count*70/100):
		print "70% of reads cleaning completed..."
		#Verbose("70% of reads cleaning completed...")
	elif percent_count==int(total_count*80/100):
		print "80% of reads cleaning completed..."
		#Verbose("80% of reads cleaning completed...")
	elif percent_count==int(total_count*90/100):
		print "90% of reads cleaning completed..."
		#Verbose("90% of reads cleaning completed...")
		
	elif percent_count==int(total_count):
		#Verbose("100% of reads cleaning completed...")
		print "100% of reads cleaning completed..."

	#total_reads_input+=1
	#identifier="@"+record.id.upper()
	identifier="@"+record.description.upper()
	qual_identifier="+"+record.id
	
	if remove_key_supplied==True:
		nt_seq=str(record.seq).upper()[key_length:]
		scores=record.letter_annotations["phred_quality"][key_length:]
	else:	
		nt_seq=str(record.seq).upper()
		scores=record.letter_annotations["phred_quality"]
	
	seq_len_original=len(nt_seq)
	if seq_len_original > inputfile_max_length:
		inputfile_max_length=seq_len_original
		
	seq_mean_original=find_mean(scores)
	
	every_seq_len_original.append(seq_len_original)
	every_seq_mean_original.append(seq_mean_original)



	last_removal=0
	first_removal=0
	last_removal_low_quality=0
	first_removal_low_quality=0
	
	if user_mode_option==1 or user_mode_option==2:
		
		while(True):	#This while loop just to trim 3' 0 values from the end. 0 is for non identified bases N
			if scores[len(scores)-1]==0:
				del scores[len(scores)-1]
				last_removal+=1
			else:
				break
		nt_seq=nt_seq[0:len(scores)]
		if 0 in scores:
			zero_in_middle+=1
		if user_mode_option==1:
			(nt_seq, scores)=get_zero_positions(identifier,nt_seq, qual_identifier, scores)	#trim zeros in the middle
			if len(nt_seq)!=len(scores):
				#print "After trimming zeros in middle the seq len and qual len is not same"
				Verbose("Sequence length and quality length are not same")
				sys.exit(1)
		
		#print "calling function three prime low quality trim"
		original_len=len(nt_seq)
		while(True):
			(last_removal,nt_seq,scores)=three_prime_low_quality_trim(identifier, nt_seq, qual_identifier, scores, last_removal)
			#(nt_seq,scores,last_removal)=trim_until_mean(nt_seq,scores,last_removal)
			if original_len!=len(nt_seq):
				original_len=len(nt_seq)
			else:
				break
		#print "returning from function three prime low quality trim"
		
		
		if str(len(nt_seq)) in length_dict.keys():	#storing the len and no. of reads with that len in dictionary
			length_dict[str(len(nt_seq))]+=1
		else:
			length_dict[str(len(nt_seq))]=1
			
		trimmed_at_3.append(last_removal)
		trimming_handle.write(identifier+ ":\t"+ str(last_removal)+"\n")

		decision_discard=decision_making(identifier, nt_seq, qual_identifier, scores, first_removal, last_removal)
		if decision_discard==False:
			total_reads_output+=1
			
		
	
	if user_mode_option==3 or user_mode_option==4:	#removes 5' and 3' zeros
		base_remove_at_ends_handle=open(output_fastq+"_end_trimmings", "w")
		base_remove_at_ends_handle.write("5'\t3'\n")
		length_before_trimming=len(nt_seq)
		both_side_zeros=True
		first_removal=0
		last_removal=0
		#print len(nt_seq),len(scores)
		while(both_side_zeros):
			if scores[0]==0 or scores[len(scores)-1]==0:
				if scores[0]==0:
					del scores[0]
					first_removal+=1
				if scores[len(scores)-1]==0:
					del scores[len(scores)-1]
					last_removal+=1
			else:
				both_side_zeros=False
				
		nt_seq=nt_seq[first_removal:length_before_trimming-last_removal]
		if 0 in scores:
			zero_in_middle+=1
		#print first_removal, last_removal
		#print len(nt_seq), len(scores)
		if user_mode_option==3:
			(nt_seq, scores)=get_zero_positions(identifier,nt_seq, qual_identifier, scores)	#trim zeros in the middle
			if len(nt_seq)!=len(scores):
				#print "after trimming zero positions, len of seq and scores not same"
				Verbose("After trimming zero positions, len of seq and scores not same")
				sys.exit(1)
		
		#print "quality scores: ",scores
		original_len=len(nt_seq)
		while(True):
			(first_removal, last_removal,nt_seq, scores)=five_three_prime_low_quality_trim(identifier, nt_seq,scores,first_removal, last_removal,len(nt_seq))
			#(nt_seq,scores,last_removal)=trim_until_mean(nt_seq,scores,last_removal)
			if original_len!=len(nt_seq):
				original_len=len(nt_seq)
			else:
				break

		trimmed_at_5.append(first_removal)
		trimmed_at_3.append(last_removal)		
		
		if str(len(nt_seq)) in length_dict.keys():	#storing the len and no. of reads with that len in dictionary
			length_dict[str(len(nt_seq))]+=1
		else:
			length_dict[str(len(nt_seq))]=1
			
		trimming_handle.write(identifier+"\n"+str(first_removal)+ "\t" + str(last_removal)+ "\n")
			#nt_seq=nt_seq[first_removal_low_quality:length_before_trimming-last_removal_low_quality]
			


		decision_discard=decision_making(identifier, nt_seq, qual_identifier, scores,first_removal, last_removal)

		if decision_discard==False:
			total_reads_output+=1



	if user_mode_option==5 or user_mode_option==6:	#split the reads selecting only the sub sequence with all scores greater than or equal to the mean quality
		seq_id=record.id
		good_seq_list=split_good_seq(seq_id, nt_seq, scores)
		if len(good_seq_list)>1:
			for count in range(len(good_seq_list)):
				
				if user_output_format==1:
					user_output_format_1("@"+seq_id+"_fragment_"+str(count)+"\n", good_seq_list[count][0]+ "\n", "+\n", convert_intqual_to_stringqual(good_seq_list[count][1])+"\n")
				if user_output_format==2:
					user_output_format_2("@"+seq_id+"_fragment_"+str(count)+"\n", good_seq_list[count][0]+ "\n", "+\n", convert_int_to_ascii_char(good_seq_list[count][1])+"\n")
				if user_output_format==3:
					user_output_format_3(">"+seq_id+"_fragment_"+str(count)+"\n", good_seq_list[count][0]+ "\n", ">"+seq_id+"_fragment_"+str(count)+"\n", convert_intqual_to_stringqual(good_seq_list[count][1])+"\n")
				
				#outfile_handle.write("@"+seq_id+"_fragment_"+str(count)+"\n")
				#outfile_handle.write(good_seq_list[count][0]+ "\n")
				#outfile_handle.write("+\n")
				#outfile_handle.write(convert_intqual_to_stringqual(good_seq_list[count][1])+"\n")
				#output_lengths.append(len(good_seq_list[count][1]))
				total_reads_output+=1
				output_means.append(find_mean(good_seq_list[count][1]))
		elif len(good_seq_list)==1:
			#decision_making("@"+seq_id, good_seq_list[0][0],"+", good_seq_list[0][1], 0,0, original_nt_seq, original_scores)
			if user_output_format==1:
				user_output_format_1("@"+seq_id+"\n", good_seq_list[0][0]+ "\n", "+\n", convert_intqual_to_stringqual(good_seq_list[0][1])+"\n")
			if user_output_format==2:
				user_output_format_2("@"+seq_id+"\n", good_seq_list[0][0]+ "\n", "+\n", convert_int_to_ascii_char(good_seq_list[0][1])+"\n")
			if user_output_format==3:
				user_output_format_3(">"+seq_id+"_fragment_"+str(count)+"\n", good_seq_list[0][0]+ "\n", ">"+seq_id+"_fragment_"+str(count)+"\n", convert_intqual_to_stringqual(good_seq_list[0][1])+"\n")

			#outfile_handle.write("@"+seq_id+"\n")
			#outfile_handle.write(good_seq_list[0][0]+ "\n")
			#outfile_handle.write("+\n")
			#outfile_handle.write(convert_intqual_to_stringqual(good_seq_list[0][1])+"\n")
			#output_lengths.append(len(good_seq_list[0][1]))
			total_reads_output+=1
			output_means.append(find_mean(good_seq_list[0][1]))

		#else:
		#	return total_reads_output




#*****************************************************************************************************
if user_output_format==1 or user_output_format==2:
	outfile_handle.close()
if user_output_format==3:
	outfile_seq_handle.close()
	outfile_qual_handle.close()
	
	
#print "Total reads input: ", total_count
Verbose("Total reads input: "+ str(total_count))
stat_handle.write("Total reads input: "+str(total_count)+"\n")
#print "Total reads output: ", total_reads_output
Verbose("Total reads output: "+ str(total_reads_output))
stat_handle.write("Total reads output: "+str(total_reads_output)+"\n")
if user_mode_option==5 or user_mode_option==6:
	sys.exit(0)
if len(output_lengths)==0:
	print "The total output reads is Zero...No further statistics and analysis"
	sys.exit(0)
	
#print "Maximum read length in output: ", max(output_lengths)
Verbose("Maximum read length in output: " +str(max(output_lengths)))
stat_handle.write("Maximum read length in output: "+ str(max(output_lengths))+"\n")
#print "Minimum read length in output: ", min(output_lengths)
Verbose("Minimum read length in output: " + str(min(output_lengths)))
stat_handle.write("Minimum read length in output: "+str(min(output_lengths))+"\n")
#print "Mean read length in output: ", sum(output_lengths)/len(output_lengths)
Verbose("Mean read length in output: "+ str(sum(output_lengths)/len(output_lengths)))
stat_handle.write("Mean read length in output: "+ str(sum(output_lengths)/len(output_lengths))+"\n")

	
#*****************************************************************************************************


#In this section, the readlength and the number of reads with that readlength is counted

len_read_pair=[]
for len_read in length_dict.keys():
	if int(len_read)>=min_length:
		len_read_pair.append((int(len_read), length_dict[len_read]))

sorted_len_read_pair=sorted(len_read_pair, key=itemgetter(0))
#writing the readlength and number of reads pair to file
readlength_read_pair_handle=open(output_fastq+"_readlength_read_pair.txt", "w")
readlength_read_pair_handle.write("Read_Length\tNumber of Reads\n")
for x, y in [(each_pair) for each_pair in sorted_len_read_pair]:
	readlength_read_pair_handle.write(str(x)+"\t"+str(y)+"\n")

readlength_read_pair_handle.close()
#*****************************************************************************************************



sorted_output_lengths=sorted(output_lengths)
cumulative_output_lengths=[0,0]
for score in sorted_output_lengths:
	cumulative_output_lengths.append(cumulative_output_lengths[len(cumulative_output_lengths)-1]+score)

del cumulative_output_lengths[0:2]
n50_number=cumulative_output_lengths[len(cumulative_output_lengths)-1]/2
n80_number=cumulative_output_lengths[len(cumulative_output_lengths)-1]*8/10
n90_number=cumulative_output_lengths[len(cumulative_output_lengths)-1]*9/10
N50_read_length=0
N80_read_length=0
N90_read_length=0
for val in range(len(cumulative_output_lengths)-1):
	if cumulative_output_lengths[val]>n50_number:
		N50_read_length=sorted_output_lengths[val]
		break

'''for val in range(len(cumulative_output_lengths)-1):
	if cumulative_output_lengths[val]>n80_number:
		N80_read_length=sorted_output_lengths[val]
		break

for val in range(len(cumulative_output_lengths)-1):
	if cumulative_output_lengths[val]>n90_number:
		N90_read_length=sorted_output_lengths[val]
		break
'''

#print "More statistics on Output Read lengths:\n"
#print "N50 read length = ", N50_read_length
Verbose("N50 read length = "+str(N50_read_length))
#print "N80 read length = ", N80_read_length
#print "N90 read length = ", N90_read_length, "\n\n"



if len(trimmed_at_5)==0:
	#print "Zero bases trimmed at 5' end "
	Verbose("Zero bases trimmed at 5' end ")
else:
	#print "Minimum trimming at 5'end: "+str(sorted(trimmed_at_5)[0])+ " Maximum trimming at 5' end: "+str(sorted(trimmed_at_5)[len(trimmed_at_5)-1])
	Verbose("Minimum trimming at 5'end: "+str(sorted(trimmed_at_5)[0])+ " Maximum trimming at 5' end: "+str(sorted(trimmed_at_5)[len(trimmed_at_5)-1]))

if len(trimmed_at_3)==0:
	#print "Zero bases trimmed at 3' end"
	Verbose("Zero bases trimmed at 3' end")
else:
	#print "Minimum trimming at 3' end: "+str(sorted(trimmed_at_3)[0])+ " Maximum trimming at 3' end: "+str(sorted(trimmed_at_3)[len(trimmed_at_3)-1])
	Verbose("Minimum trimming at 3' end: "+str(sorted(trimmed_at_3)[0])+ " Maximum trimming at 3' end: "+str(sorted(trimmed_at_3)[len(trimmed_at_3)-1]))

#print "\n\nGetting the read length distribution after the quality trimming" 
Verbose("\n\nGetting the read length distribution after the quality trimming")
get_readlen_info(output_lengths)

def get_dict_from_numarray(array):
	mydict={}
	for eachnum in array:
		if eachnum in mydict:
			mydict[eachnum]+=1
		else:
			mydict[eachnum]=1
	
	return mydict		

def boxplotting(fastqfilename, interval, plot_title, plot_name):

	from Bio import SeqIO
	input_handle=open(fastqfilename)
	input_reads_scores=[]
	minimum_length=10000	#at least a value is less than 10000
	maximum_length=0
	average_length=0
	sum=0
	total=0
	for record in SeqIO.parse(input_handle, 'fastq'):
		record_score=record.letter_annotations["phred_quality"]
		input_reads_scores.append(record_score)
		if len(record_score)>maximum_length:
			maximum_length=len(record_score)
		if len(record_score)<minimum_length:
			minimum_length=len(record_score)
		sum+=len(record_score)
		total+=1

	average_length=sum/total

	#import operator
	#input_reads_scores.sort(key=lambda x: x[1], reverse=True)


	#for each in input_reads_scores:
	#	print each


	#find 1st quartile and 3 quartile values
	def first_quartile(scores):
		scores=sorted(scores)
		first_quartile_position=len(scores)/4
		return scores[first_quartile_position]
	
	def third_quartile(scores):
		scores=sorted(scores)
		third_quartile_position=len(scores)*3/4
		return scores[third_quartile_position]
	
	
	#find max length array
	max_len=0
	for array in input_reads_scores:
		if len(array)>max_len:
			max_len=len(array)

	#print max_len
	#interval=int(sys.argv[2])
	first_quartile_list=[]
	third_quartile_list=[]
	average_list=[]
	scores_at_position=[]
	array_of_first_quartiles=[]
	array_of_third_quartiles=[]
	array_of_max_score=[]
	array_of_min_score=[]
	all_scores_at_position=[]
	xlabels=[]
	num_of_seqs=[]
	for position in range(0, max_len, interval):
		sum=0
		total_num=0
		avg=0
		scores_at_position=[]
		xlabels.append(str(position))
		for myarray in input_reads_scores:
			if position<len(myarray):
				sum+=myarray[position]
				total_num+=1
				scores_at_position.append(myarray[position])
	
		all_scores_at_position.append(scores_at_position)
		num_of_seqs.append(len(scores_at_position))
	
		avg=sum/total_num
		average_list.append(avg)
		#array_of_first_quartiles.append(first_quartile(scores_at_position))
		array_of_third_quartiles.append(third_quartile(scores_at_position))
		array_of_max_score.append(max(scores_at_position))
		array_of_min_score.append(min(scores_at_position))
	
	

	#print "average scores",average_list
	####plotting average base score
	#plot_average_base_score(xlabels,all_scores_at_position, array_of_max_score, array_of_min_score, array_of_third_quartiles, sys.argv[3], "eps")
	#print "First quartile", array_of_first_quartiles
	#print "Third quartile", array_of_third_quartiles
	#print "Minimum length" , minimum_length
	#print "Maximum length", maximum_length
	#print "Average length", average_length





	#x = [22883,16525,20577,570,11947,12943]		#reads output
	#y1 = [162,89,175,106,75,113]		#mean_length
	#y1=[90.2164,30.4302,91.7823, 74.2899,19.7737,34.467]	#std deviation of read lengths
	#y2 = [216332,539,252979,616,0,0]		#total bad score bases
	#cpu_time=[2383.8094,4052.476,976.5002,13508.2533,20262.38,11578.5029]	#Bases per 1/100 of sec
	#algorithms_by_time=['R', 'C', 'P', 'F', 'H', 'G']
	#algorithms=['Ram', 'Clean_reads', 'PRINSEQ', 'Fastx', '454HIV','Geneious']

	if 1:
		fig = plt.figure(figsize=(15,5),facecolor='w')
		#plt.xticks(range(0,len(all_scores_at_position),10), xlabels)
		#plt.xticks(range(0,120,10), range(0,1200,100))
	 	plt.title(plot_title)
		# Pick some colors
		#palegreen = matplotlib.colors.colorConverter.to_rgb('#8CFF6F')
		#paleblue = matplotlib.colors.colorConverter.to_rgb('#708DFF')
	 
		# Plot response time
		ax1 = fig.add_subplot(111)
	
		bp = ax1.boxplot(all_scores_at_position, notch=0, sym='', vert=1, whis=1.5,positions=None, widths=None, patch_artist=True)
	
	
		#ax2.set_xticklabels(range(0,1200,100))
	
	
	
		#ax1.set_yscale('log')
	
		majorLocator   = MultipleLocator(10)	
		ax1.set_ylabel("Quality Score")
		ax1.set_xlabel("Read Base postion (Bin size : 10bp) ")
		ax2 = ax1.twinx()
		ax2.plot(num_of_seqs, 'o-', color='g', linewidth=2, markersize=4)
		
		#ax2.axis([0,inputfile_max_length/interval,0,inputfile_total_seqs])
		ax2.set_ylim(0, inputfile_total_seqs)
		ax2.set_ylabel("Number of sequences")
		ax2.xaxis.set_major_locator(majorLocator)
		#ax2.set_xlim(0,inputfile_max_length)
		#ax1.xaxis.set_major_locator(MaxNLocator(15))
		
		#ax1.axis([0,inputfile_max_length/interval,0,40])
		ax1.set_ylim(0,40)
		ax1.set_xlim(0,inputfile_max_length/interval)
		#ax1.axis('tight') #this does not allow using set_xlim
		
		
		ax1.xaxis.set_major_locator(majorLocator)
		tick_labels=[str(tickposition) for tickposition in range(0,inputfile_max_length, 100)]
		
		ax1.set_xticklabels(tick_labels)
		
		
		
		#ax1.set_xticklabels(xlabels)
		#ax1.xticks(range	(0,len(all_scores_at_position),10), xlabels)
	
		# Tweak colors on the boxplot
		#plt.setp(bp['boxes'], color='g')
		#plt.setp(bp['whiskers'], color='g')
		#plt.setp(bp['medians'], color='black')
		#plt.setp(bp['fliers'], color=palegreen, marker='+')


	 
		# Plot throughput
	
	 
	
	    	plt.draw()
	    	plt.savefig(plot_name+"."+format_plot, format=format_plot)





def plot_analysis(every_seq_len_original,every_seq_mean_original, output_lengths, output_means, output_fastq):
	
	#print "Drawing the histogram of read length vs counts for Original reads before trimming" 
	Verbose("Drawing the histogram of read length vs counts for Original reads before trimming")
	diff=max(every_seq_len_original)-min(every_seq_len_original)
	if diff!=0:
		majorLocator=MultipleLocator(50)
		mydict=get_dict_from_numarray(every_seq_len_original)
		fig=plt.figure(figsize=(15,5),facecolor='w')
		ax = fig.add_subplot(111)
		plt.title("Read Length Distribution before trimming")
		plt.bar(mydict.keys(), mydict.values(), color='black')
		plt.xlabel("Read Length")
		plt.ylabel("Number of Reads")
		ax.xaxis.set_major_locator(majorLocator)
		ax.axis('tight')
		plt.savefig(output_fastq+"_Before_trimming_Readlength_vs_counts."+format_plot, format=format_plot)
	else:
		#print "Readlength Vs Counts before trimming could not be drawn, probably all reads have same length" 
		Verbose( "Readlength Vs Counts before trimming could not be drawn, probably all reads have same length")
	plt.clf()
	

	#print "Drawing the histogram of output read length vs counts" 
	Verbose("Drawing the histogram of output read length vs counts")
	diff=max(output_lengths)-min(output_lengths)
	if diff!=0:
		majorLocator   = MultipleLocator(50)
		mydict=get_dict_from_numarray(output_lengths)
		fig=plt.figure(figsize=(15,5),facecolor='w')
		ax = fig.add_subplot(111)
		plt.title("Read Length Distribution after trimming")
		plt.bar(mydict.keys(), mydict.values(), color='black')
		plt.xlabel("Read Length")
		plt.ylabel("Number of Reads")
		ax.xaxis.set_major_locator(majorLocator)
		ax.axis('tight')
		plt.savefig(output_fastq+"_After_trimming_Readlength_vs_counts."+format_plot, format=format_plot)
		
	else:
		#print "Readlength Vs Counts after trimming could not be drawn, probably all reads have same length"
		Verbose("Readlength Vs Counts after trimming could not be drawn, probably all reads have same length")
	plt.clf()
	
	#print "Drawing the histogram of read mean vs counts for Original reads before trimming" 
	Verbose("Drawing the histogram of read mean vs counts for Original reads before trimming")
	diff=max(every_seq_mean_original)-min(every_seq_mean_original)
	if diff!=0:
		majorLocator   = MultipleLocator(1)
		mydict=get_dict_from_numarray(every_seq_mean_original)
		fig=plt.figure(figsize=(15,5),facecolor='w')
		ax = fig.add_subplot(111)
		plt.title("Read Mean Distribution before trimming")
		plt.bar(mydict.keys(), mydict.values(), align="center", color='black')
		plt.xlabel("Read Mean")
		plt.ylabel("Number of Reads")
		ax.xaxis.set_major_locator(majorLocator)
		ax.axis('tight')
		plt.savefig(output_fastq+"_Before_trimming_ReadMean_vs_ReadCount."+format_plot, format=format_plot)
	else:
		#print "Read Mean Vs Read Counts before trimming could not be drawn, probably all reads have same length"
		Verbose("Read Mean Vs Read Counts before trimming could not be drawn, probably all reads have same length")
	plt.clf()
	
	#print "Drawing the histogram of output read mean vs counts after trimming" 
	Verbose("Drawing the histogram of output read mean vs counts after trimming")
	diff=max(output_means)-min(output_means)
	if diff!=0:
		majorLocator   = MultipleLocator(1)
		mydict=get_dict_from_numarray(output_means)
		fig=plt.figure(figsize=(15,5),facecolor='w')
		ax = fig.add_subplot(111)
		plt.title("Read Mean Distribution after trimming")
		plt.bar(mydict.keys(), mydict.values(), align='center', color='black')
		plt.xlabel("Read Mean")
		plt.ylabel("Number of Reads")
		ax.xaxis.set_major_locator(majorLocator)
		ax.axis('tight')
		plt.savefig(output_fastq+"_After_trimming_ReadMean_vs_ReadCount."+format_plot, format=format_plot)
	else:
		#print "Read Mean Vs Read Counts after trimming could not be drawn, probably all reads have same length"
		Verbose("Read Mean Vs Read Counts after trimming could not be drawn, probably all reads have same length")
	
	plt.clf()


#*****************************************************************************************************
#closing all files

file_handle.close()
read_with_zeros_handle.close()
position_of_zero_handle.close()

#*****************************************************************************************************

if plot_supplied==True:
	print "Analytical Plotting now..."
		
	
	try:
		import numpy as np
	except ImportError:
		print "numpy module not found or not installed. Sorry cannot plot"
		sys.exit(1)
		

	try:
		import matplotlib
		matplotlib.use('Agg')
		import matplotlib.pyplot as plt
		from mpl_toolkits.axes_grid.parasite_axes import SubplotHost
		from matplotlib.ticker import MultipleLocator
		
	except ImportError:
		print "pyplot module not found in matplotlib or matplotlib not installed. Sorry cannot plot"
		sys.exit(1)
		
	plot_analysis(every_seq_len_original,every_seq_mean_original, output_lengths, output_means, output_fastq)
	#print inputfile_max_length
	#print inputfile_total_seqs
	
	if user_output_format==2:
		boxplotting(fastq_file, 10, "Trend of Quality Score before Trimming", outputfile_dirname + "/Quality_trend_before_trimming")

		boxplotting(output_fastq+".fastq", 10, "Trend of Quality Score after Trimming", outputfile_dirname+ "/Quality_trend_after_trimming")
	elif user_output_format==3:
		#print "Fasta file name ", output_fastq+".fasta"
		#print "Qual file name", output_fastq+ ".qual"
		
		fastq_intermediate=open("Intermediate_Fastq.fastq", "w")
		print "Converting output fasta and qual file to intermediate fastq file"
		fastq_records=PairedFastaQualIterator(open(output_fastq+".fasta"), open(output_fastq+ ".qual"))
		print "Done"
		SeqIO.write(fastq_records, fastq_intermediate, "fastq")
		fastq_intermediate.close()
		
		
		boxplotting(fastq_file, 10, "Trend of Quality Score before Trimming", outputfile_dirname + "/Quality_trend_before_trimming")
		boxplotting("Intermediate_Fastq.fastq", 10, "Trend of Quality Score after Trimming",outputfile_dirname + "/Quality_trend_after_trimming")
		import os
		os.system("rm Intermediate_Fastq.fastq")
	elif user_output_format==1:
		boxplotting(fastq_file, 10, "Trend of Quality Score before Trimming", outputfile_dirname+ "/Quality_trend_before_trimming")
		print "Converting fastq int file to fastq ascii file"
		convert_fastq_int_2_fastq_ascii(output_fastq)
		print "Done"
		boxplotting("Intermediate_Fastq.fastq", 10, "Trend of Quality Score after Trimming", outputfile_dirname + "/Quality_trend_after_trimming")
		import os
		os.system("rm Intermediate_Fastq.fastq")

#print "Program completed at ", datetime.now()
if user_mode_option==2 or user_mode_option==4:
	#print output_fastq + "_ZeroPositions"
	#print output_fastq + "_Reads_with_zeros_in_middle"
	os.system("rm " + "\"" + output_fastq + "_ZeroPositions\"")
	os.system("rm " + "\"" + output_fastq + "_Reads_with_zeros_in_middle\"")
	#os.system("rm " + output_fastq + "_end_trimmings")
	
Verbose("Program completed at "+ str(datetime.now()))

print "Completed....Terminating now"
#Finally sys.exit the program
sys.exit(0)
