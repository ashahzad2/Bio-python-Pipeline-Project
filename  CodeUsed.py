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
