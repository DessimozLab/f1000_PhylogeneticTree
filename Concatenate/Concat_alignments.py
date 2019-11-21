import sys

## Get all files of the path passed as argument

from os import listdir
from os.path import isfile, join
onlyfiles = [f for f in listdir(sys.argv[1]) if isfile(join(sys.argv[1], f))]


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
        
if __name__="__main__":        
    if len(sys.argv) < 3:
	sys.exit('please specify a path to the folder containing the alignments and an output format: fasta or phylip')

    if sys.argv[2]!= 'fasta' and sys.argv[2]!='phylip':
    sys.exit('please specify an output format: fasta or phylip')


    ## First it will go through all the files in order to get all headers
    IDs_SEQs_dict={}
    for file in onlyfiles:
        with open(sys.argv[1] + file) as f:
            for line in f:
                if line.startswith('>'):
                    if line.rstrip() not in IDs_SEQs_dict :
                        IDs_SEQs_dict[line.rstrip()] = ''


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
            sys.exit('files not aligned')
        
        
    if sys.argv[2]=='fasta':
        print_fasta_dict(IDs_SEQs_dict) 
    if sys.argv[2]=='phylip':   
        print_PHYLIP_dict(IDs_SEQs_dict)
