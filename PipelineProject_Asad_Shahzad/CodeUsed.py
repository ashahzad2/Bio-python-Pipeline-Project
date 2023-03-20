#Retrieving fasta records from NCBI database
from Bio import Entrez 
Entrez.email = "youremail.com"
handle = Entrez.efetch(db="nucleotide", id=["accession_id#"], rettype="fasta") #searching the specific accession id in the nucleotide database in fasta format
records = handle.read() #read in the file
print (records) 
with open ("file_name.fasta", "w") as f:
    f.write(records) #write the records in a fasta file

    
    

#Calculating the number of contigs with a length greater than 1000
from Bio import SeqIO #Import the SeqIO class to be able to parse the given fasta file
contigFile = "contigs.fasta" #use the contigs.fasta file given from SPAdes
contigNum = 0 #variable for tracking the contigs

for record in SeqIO.parse(contigFile, "fasta"): #iterate through each record in the file
    if len(record.seq) > 1000: #If the length of the record is greater than 1000...
             contigNum += 1 #then add to the contigs variable and move to the next record

print ("This is the number of contigs with a length > 1000: " , contigNum)
#Prints the number of contigs




#Calculating the length of the assembly
from Bio import SeqIO #Import the SeqIO class to be able to parse the given fasta file
contigFile = "contigs.fasta" #use the contigs.fasta file given from SPAdes
totalLength = 0 #variable for tracking the total length

for record in SeqIO.parse(contigFile, "fasta"): #iterate through each record in the file
     if len(record.seq) > 1000: #If the length of the record is greater than 1000...
             totalLength += len(record.seq) #then add it to the length of the variable
             
print("This is the number of bp in the assembly: ", totalLength)
# prints the number of base pairs in the assembly




#Finding the longest contig
from Bio import SeqIO #Import the SeqIO class to be able to parse the given fasta file
 
contigFile = "contigs.fasta" #use the contigs.fasta file given from SPAdes
contigLongest = 0 #variable for tracking the longest contig
contigRecord = None #initalize the records as None
 
for record in SeqIO.parse(contigFile, "fasta"): #iterate through each record in the file
     if len(record.seq) > contigLongest: #If the length of the record is greater than the longest contig..
             contigRecord = record #then that record becomes the current record
             contigLongest = len(record.seq) #and that record's length is now in the longest contig variable
 
print("The longest contig from the assembly is: ", contigRecord.seq)
#prints the entire sequence of the longest contig
print("The longest contig from the assembly is: ", contigRecord.id)
#prints id of the contig




#BLASTing from Python
input_file = "contig.fasta" #contructing a string for the contigs file
output_file = "results.csv" #contructing a string for the output file
blast_cmd = 'blastn -query ' +input_file+ -db Betaherpesvirinae -out '+output_file+' -outfmt "6 sacc pident length qstart qend sstart send bitscore evalue stitle" 
#this blast command specifies the type of blast we are running (blastn), what the query file is, the local database, output file and the outfmt output
#the outfmt output allows a tab delimited output and allows us to specify the information we want from BLAST for our top ten hits

#the blast command is run on an os system call, which takes the string and run through terminal
import os
os.system(blast_cmd)
#creates the BLAST results in the working directory as results.csv



