# Bio-python-Pipeline-Project
COMP 383 001
Project
-	Human herpesvirus 5 is also known as Human cytomegalovirus and is typically abbreviated as HCMV.

Track 2: Genome Assembly
-	We want to compare HCMV transcriptomes 2- and 6-days post infection (dpi)

	o A transcriptome is the sum total of all the messenger RNA molecules expressed from the genes of an organism.
	
-	We will achieve this through genome assembly (Track 2)

       o Mapping: Reads are mapped to the available known or reference genome rather than assembled to each other
       
       o De novo Assembly: Reads are constructed into longer sequences (contigs then scaffolds)

-	Which strains are most similar to these patient samples?

-	To compare to other strains, assemble these transcriptome reads.

-	All code and its comments found in main folder.


Commands:

1.	Retrieved the 4 HCMV transcriptomes.

a.	Look up each provided SRA file on NCBI.

b.	Click on the Run ID

c.	Click on Data Access

d.	Copy the AWS location link address.

e.	Wget the link address on desired directory

f.	Repeated these steps for all of the donor files.

g.	Uncompressed the data using fastq dump.

h.	Run the following command on terminal: fastq-dump -I –split-files SRR(File)

i.	Record the number of read pairs of each file.



2.	Genome Assembly
	
	This assembly will produce a genome enough to be useful in BLAST.

a.	Search the Accession ID: NC_006273.2 on NCBI which contains the genome of HCMV.

b.	Retrieve the HCMV genome using the following python code:


>>> from Bio import Entrez
>>> 
>>> Entrez.email = "youremail.com"
>>> 
>>> handle = Entrez.efetch(db="nucleotide", id=["accession_id#"], rettype="fasta")
>>> 
>>> records = handle.read()
>>> 
>>> print (records)
>>> 
>>> with open ("file_name.fasta", "w") as f:
>>> 
... 	f.write(records)

	

c.	This code retrieves the genome from NCBI in fasta format 

d.	Copy the file into desired terminal directory 

e.	Then to build the “index” to map to, we have to use Bowtie2, a software tool for creating reference assemblies. Can be found here (http://bowtie-bio.sourceforge.net/bowtie2/index.shtml). Use the bowtie2 command: bowtie2-build file_name.fasta index_name.

f.	This creates 6 files labeled:
index_name.1.bt2, index_name.2.bt2, index_name.3.bt2, index_name.4.bt2, index_name.rev.1.bt2, index_name.rev.2.bt2

g.	Now you can index with the following command: bowtie –quiet -x index_name -1 SRR5660030_1.fastq -2 SRR5660030_2.fastq -S index_name_map1.sam

h.	Then run a bowtie2 command that writes out just reads that map with the following command: bowtie –quiet -x index_name -1 SRR5660030_1.fastq -2 SRR5660030_2.fastq -S index_name_map2.sam –al -conc-gz SRR5660030_1fastq_%.fq.gz 

i.	Repeat for all data until all reads have been mapped and recorded the following from the assembly:

Donor 1 (2dpi): Had 2259287 read pairs before Bowtie2 filtering and had 1827537 read pairs after. (80.99% overall alignment rate)

Donor 1 (6dpi): Had 2004530 read pairs before Bowtie2 filtering and had 1423416 read pairs after. (71.01% overall alignment rate)

Donor 3 (2dpi): Had 2730258 read pairs before Bowtie2 filtering and had 1994180 read pairs after. (73.04% overall alignment rate)

Donor 3 (6dpi): Had 2476889 read pairs before Bowtie2 filtering and has 1319686 read pairs after. (53.28% overall alignment rate)


3.	Using Bowtie2 output reads, assemble all 4 transcriptomes together to produce 1 assembly via SPAdes. SPAdes is a fast de novo assembler that can handle sequencing reads. Can be download here(http://cab.spbu.ru/software/spades/)

Command:

spades.py -k 77,99,127 -t 2 --only-assembler 
--pe-1 1 SRR5660030_1.fastq --pe-2 1 SRR5660030_2.fastq 
--pe-1 2 SRR5660033_1.fastq --pe-2 2 SRR5660033_2.fastq 
--pe-1 3 SRR5660044_1.fastq --pe-2 3 SRR5660044_2.fastq 
--pe-1 4 SRR5660045_1.fastq --pe-2 4 SRR5660045_2.fastq 
-o HCMV_assembly/

4.	You will find a new directory made for SPAdes outputs in the same working directory you ran the command. SPAdes will create a "contigs.fasta" file. Use this python code to calculate the number of contigs with a length > 1000:

>>> from Bio import SeqIO
>>> 
>>> contigFile = "contigs.fasta"
>>> 
>>> contigNum = 0
>>> 
>>> for record in SeqIO.parse(contigFile, "fasta"):
>>> 
...     if len(record.seq) > 1000:
...
...             contigNum += 1
... 
>>> print ("This is the number of contigs with a length > 1000: " , contigNum)
>>> 
#This is the number of contigs with a length > 1000:  114

There are 114 contigs > 1000 bp in the assembly


Then use this python code to calculate the length of the assembly (total number of bp in all of the contigs > 1000 bp in length):

>>> from Bio import SeqIO
>>> 
>>> contigFile = "contigs.fasta"
>>> 
>>> totalLength = 0
>>> 
>>>
>>> for record in SeqIO.parse(contigFile, "fasta"):
>>> 
...     if len(record.seq) > 1000:
...
...             totalLength += len(record.seq)
...
... 
>>> print("This is the number of bp in the assembly: ", totalLength)
>>> 
#This is the number of bp in the assembly:  195531

There are 195531 bp in the assembly

5.	Does your assembly align with other virus strains? We will find the longest contig and use that as blast+ input to query the nr nucleotide database demlimited to members of the Betaherpesvirinae subfamily.

Use this python code to retrieve the longest contig from your SPAdes assembly.

>>> from Bio import SeqIO
>>> 
>>> contigFile = "contigs.fasta"
>>> 
>>> contigLongest = 0
>>> 
>>> contigRecord = None
>>> 
>>> 
>>> for record in SeqIO.parse(contigFile, "fasta"):
>>> 
...     if len(record.seq) > contigLongest:
...
...             contigRecord = record
...
...             contigLongest = len(record.seq)
...
... 
>>> print("The longest contig from the assembly is: ", contigRecord.seq)
>>> 
#The longest contig from the assembly is:  TGTGCGCTGTGTTTATTTTTTCTTCTGTGTCT…



Make a local database of just sequences from the Betaherpesvirinase sub family:
Retrieving sequences from NCBI:
	Searched txid10357[Organism:exp]
	
	Filter to show RefSeq Collection
	
	Download the sequences to computer in fasta format sorted by length.
	
Select the top 10 hits and download it as a fasta file.

Copy this file into desired directory and run the following commands to create a local database:

makeblastdb -in fileName.fasta -out dbnam -title dbname -dbtype nucl

Now run the following python commands to BLAST the longest contig(query sequence) delimited to the local database of just sequences from the Betaherpesvirinae subfamily:

>>> input_file = "contig.fasta"
>>> 
>>> output_file = "results.csv"
>>> 
>>> blast_cmd = 'blastn -query ' +input_file+ -db Betaherpesvirinae -out '+output_file+' -outfmt "6 sacc pident length qstart qend sstart send bitscore evalue stitle"
>>> 
>>> import os
>>> 
>>> os.system(blast_cmd)
>>> 

BLAST results will be created in the working directory as results.csv. It will display the 10 ten hits, which are the virus strains your assembly closest to.
The following is displayed for your top 10 hits:
Subject accession, Percent identity, Alignment length, Start of alignment in query, End of alignment in query, Start of alignment in subject, End of alignment in subject, Bit score, E-value, and Subject title.

