# Bio-python-Pipeline-Project
COMP 383 001
Project Info
-	Human herpesvirus 5 is also known as Human cytomegalovirus and is typically abbreviated as HCMV.

Track 2: Genome Assembly
-	We want to compare HCMV transcriptomes 2- and 6-days post infection (dpi)
o	A transcriptome is the sum total of all the messenger RNA molecules expressed from the genes of an organism.
-	We will achieve this through genome assembly (Track 2)
o	Mapping: Reads are mapped to the available known or reference genome rather than assembled to each other
o	De novo Assembly: Reads are constructed into longer sequences (contigs then scaffolds)
-	Which strains are most similar to these patient samples?
-	To compare to other strains, assemble these transcriptome reads.

Commands:
1.	Retrieved the 4 HCMV transcriptomes.
a.	Look up each provided SRA file on NCBI.
b.	Click on the Run ID
c.	Click on Data Access
d.	Copy the AWS location link address.
e.	Wget the link address on pipeline project directory
f.	Repeated these steps for all of the donor files.
g.	Uncompressed the data using fastq dump.
h.	Fastq-dump -I –split-files SRR(File)
i.	Record the number of read pairs of each file before Bowtie2.


2.	Genome Assembly
a.	Search the Accession ID: NC_006273.2 on NCBI which contains the genome of HCMV.
b.	Retrieve the HCMV genome using the following python code:

from Bio import Entrez
Entrez.email = "ashahzad2@luc.edu"
handle = Entrez.efetch(db="nucleotide", id=["NC_006273.2"], rettype="fasta")
records = handle.read()
print (records)
with open ("HCMV.fasta", "w") as f:
 f.write(records)

c.	This code retrieves the genome from NCBI fasta format 
d.	Copy the file into desired directory 
e.	Then to build the “index” to map to, use the bowtie2 command: bowtie2-build HCMV.fasta HCMV.
f.	This creates 6 files labeled:
HCMV.1.bt2, HCMV.2.bt2, HCMV.3.bt2, HCMV.4.bt2, HCMV.rev.1.bt2, HCMV.rev.2.bt2
g.	Now you can index with the following command: bowtie –quiet -x HCMV -1 SRR5660030_1.fastq -2 SRR5660030_2.fastq -S HCMVmap1.sam
h.	Then run a bowtie2 command that writes out just reads that map with the following command: bowtie –quiet -x HCMV -1 SRR5660030_1.fastq -2 SRR5660030_2.fastq -S HCMVmap2.sam –al -conc-gz SRR5660030_1fastq_%.fq.gz 
i.	Repeated for all data until all reads have been mapped and recorded the following

Donor 1 (2dpi): Had 2259287 read pairs before Bowtie2 filtering and had 1827537 read pairs after. (80.99% overall alignment rate)
Donor 1 (6dpi): Had 2004530 read pairs before Bowtie2 filtering and had 1423416 read pairs after. (71.01% overall alignment rate)
Donor 3 (2dpi): Had 2730258 read pairs before Bowtie2 filtering and had 1994180 read pairs after. (73.04% overall alignment rate)
Donor 3 (6dpi): Had 2476889 read pairs before Bowtie2 filtering and has 1319686 read pairs after. (53.28% overall alignment rate)

3.	Using Bowtie2 output reads, assemble all 4 transcriptomes together to produce 1 assembly via SPAdes.

Command:
spades.py -k 77,99,127 -t 2 --only-assembler 
--pe-1 1 SRR5660030_1.fastq --pe-2 1 SRR5660030_2.fastq 
--pe-1 2 SRR5660033_1.fastq --pe-2 2 SRR5660033_2.fastq 
--pe-1 3 SRR5660044_1.fastq --pe-2 3 SRR5660044_2.fastq 
--pe-1 4 SRR5660045_1.fastq --pe-2 4 SRR5660045_2.fastq 
-o HCMV_assembly/

4.	Python code to calculate the number of contigs with a length > 1000

>>> from Bio import SeqIO
>>> contigFile = "contigs.fasta"
>>> contigNum = 0
>>> 
>>> for record in SeqIO.parse(contigFile, "fasta"):
...     if len(record.seq) > 1000:
...             contigNum += 1
... 
>>> print ("This is the number of contigs with a length > 1000: " , contigNum)
This is the number of contigs with a length > 1000:  114

There are 114 contigs > 1000 bp in the assembly

Python code to calculate the length of the assembly (total number of bp in all of the contigs > 1000 bp in length)

>>> from Bio import SeqIO
>>> contigFile = "contigs.fasta"
>>> totalLength = 0
>>> 
>>> for record in SeqIO.parse(contigFile, "fasta"):
...     if len(record.seq) > 1000:
...             totalLength += len(record.seq)
... 
>>> print("This is the number of bp in the assembly: ", totalLength)
This is the number of bp in the assembly:  195531

There are 195531 bp in the assembly

5.	Does your assembly align with other virus strains?
Python code to retrieve the longest contig from your SPAdes assembly.

>>> from Bio import SeqIO
>>> 
>>> contigFile = "contigs.fasta"
>>> contigLongest = 0
>>> contigRecord = None
>>> 
>>> for record in SeqIO.parse(contigFile, "fasta"):
...     if len(record.seq) > contigLongest:
...             contigRecord = record
...             contigLongest = len(record.seq)
... 
>>> print("The longest contig from the assembly is: ", contigRecord.seq)
The longest contig from the assembly is:  TGTGCGCTGTGTTTATTTTTTCTTCTGTGTCT…

>>> print("The longest contig from the assembly is: ", contigRecord.id)
The longest contig from the assembly is:  NODE_1_length_9848_cov_138.362617

Take the longest contig and blast it
Using blastn

Make a local database of just sequences from the Betaherpesvirinase sub family:
Retrieving sequences from NCBI:
	Searched txid10357[Organism:exp]
	Filter to show RefSeq Collection
	Download the sequences to computer in fasta format sorted by length.
Select the top 10 hits and download it as a txt file.
