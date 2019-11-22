import sys

## Get all files of the path passed as argument

from os import listdir
from os.path import isfile, join


# Check if all sequences have the same lenght        
def Check_Alignment(file): 
    import sys
    seq=''

    Alignment=True
    with open(file) as f:
            seq=''
            first=True
            header=f.readline()
            for line in f:
                if line.startswith(">"):
                    if first:
                        ref_len=len(seq)
                        first=False
                    else:
                        if (ref_len==len(seq)) == False:
                            Alignment=False

                    header=line

                    seq=''
                else:
                    seq=seq+line.rstrip()

    if (ref_len==len(seq)) == False:
        Alignment=False
    return(Alignment)

## Function to print a dictionnary in fasta format
def print_fasta_dict(fasta_dict):
    def insert_newlines(string, every=64):
        return '\n'.join(string[i:i+every] for i in range(0, len(string), every))
    
    for c,i in enumerate(list(fasta_dict.keys())):
        print(i)
        print(insert_newlines(list(fasta_dict.values())[c]))

## Function to print a dictionnary in Phylip format
def print_PHYLIP_dict(fasta_dict):
    print(len(fasta_dict),len(list(fasta_dict.values())[0]))
    for c,i in enumerate(list(fasta_dict.keys())):
        print(i[1:6] + '\t' +list(fasta_dict.values())[c]) 
        

def get_help_and_exit(errorMsg=""):
	help = """Usage: python Concat_alignment.py alignment_folder format > output
	This script will concatenate all alignment in a folder and output a file in the given format.
	Parameters:
		alignment_folder: a folder containing all alignments file to concatenate
		format: Format in which the concatenated alignment will be outputed. Must be fasta or phylip"""
	if errorMsg:	
		error = "\n----------------------------------------\nError:    "+errorMsg
	else:
		error=""
	message = help+error

	sys.exit(message)

if __name__=="__main__":        
    
    if len(sys.argv) < 3:
	if len(sys.argv)==2 and (sys.argv[1]=="--help" or sys.argv[1]=="-h"):
		get_help_and_exit()
	get_help_and_exit('Please specify a path to the folder containing the alignments and an output format: fasta or phylip')

    if sys.argv[2]!= 'fasta' and sys.argv[2]!='phylip':
        get_help_and_exit('Please specify an output format: fasta or phylip')

    onlyfiles = [f for f in listdir(sys.argv[1]) if isfile(join(sys.argv[1], f))]

    ## First it will go through all the files in order to get all headers
    IDs_SEQs_dict={}
    wrongFormat= list()
    for file in onlyfiles:
	alignFile= False
        with open(sys.argv[1] + file) as f:
            for line in f:
                if line.startswith('>'):
		    alignFile = True
                    if line.rstrip() not in IDs_SEQs_dict :
                        IDs_SEQs_dict[line.rstrip()] = ''
	if not alignFile:
		wrongFormat.append(file)
    for wFile in wrongFormat:
	onlyfiles.remove(wFile)
    if len(onlyfiles)==0:
	get_help_and_exit('No alignment file found. Please specify a folder containing alignment files in FASTA format.')

    for file in onlyfiles:  
        #print(Check_Alignment)
        if Check_Alignment(sys.argv[1] + file):
            seq=''
            with open(sys.argv[1] + file) as f:
    
                header=f.readline()
                for line in f:
                    if line.startswith(">"):
                        IDs_SEQs_dict[header.rstrip()] += seq

                        header=line
                        seq=''

                    else:
                        seq=seq+line.rstrip()

                IDs_SEQs_dict[header.rstrip()] += seq
            # end of fasta file
            #now check if all seqs have the same length and add '-' to the ones which are too short
        
            # get the longest seq
            ref=0
            for i in IDs_SEQs_dict.values():
                if len(i) > ref:
                    ref = len(i)
                
            # add '-' to be equal to all other
            for i in IDs_SEQs_dict.keys():
                if len(IDs_SEQs_dict[i])!=ref:
                    IDs_SEQs_dict[i]+=(ref-len(IDs_SEQs_dict[i]))*'-'
                

        
        else:
            get_help_and_exit('Non-standard alignments file. The aligned sequences do not have the same size. Please check the format of your files.')
        
        
    if sys.argv[2]=='fasta':
        print_fasta_dict(IDs_SEQs_dict) 
    if sys.argv[2]=='phylip':   
        print_PHYLIP_dict(IDs_SEQs_dict)
