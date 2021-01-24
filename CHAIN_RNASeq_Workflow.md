# Check for required software on Bluewave
### Update software to latest version
- fastqc  
- MultiQC  
- fastp  
- HiSat2  
- Samtools  
- StringTie  
- gffcompare  
- Python  

# 1) Obtain Reference Genome
[Cunning et al 2018](https://www.nature.com/articles/s41598-018-34459-8)

/data/putnamlab/REFS/Pdam

[Vidal Dupiol et al ](https://www.biorxiv.org/content/10.1101/698688v3)
/data/putnamlab/REFS/Pact




### Functional annotation
files

# make folder structure
```
mkdir data
cd data

```

# Downloaded files

## path where we stored the RAW fastq.gz files
```/data/putnamlab/KITT/hputnam/CHAIN```

# 2) Check file integrity 

a) Count all files to make sure all downloaded

```
ls -1 | wc -l

35 files
```
/data/putnamlab/hputnam/CHAIN_RNASeq

b) Verify data transfer integrity with md5sum

```
nano /data/putnamlab/hputnam/CHAIN_RNASeq/scripts/check_transfer.sh
```

```
#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/hputnam/CHAIN_RNASeq

md5sum /data/putnamlab/KITT/hputnam/CHAIN/*.gz > CHAIN_URIcheckmd5.md5

```

```
sbatch /data/putnamlab/hputnam/CHAIN_RNASeq/scripts/check_transfer.sh
```

### Checksum from Sequencer
```

```


### Checksum from files on Bluewaves
```
8b8f59c8b92d3ea0f24bdd7c2caa573a  /data/putnamlab/KITT/hputnam/CHAIN/G021_S1_R1_001.fastq.gz
f911856f0094504644483edb60f244f2  /data/putnamlab/KITT/hputnam/CHAIN/G023_S8_R1_001.fastq.gz
d52ee304dba6f7e86744614addce8b57  /data/putnamlab/KITT/hputnam/CHAIN/G024_S22_R1_001.fastq.gz
a53a4279bcbd75ad59db9b5a288bb290  /data/putnamlab/KITT/hputnam/CHAIN/G028_S29_R1_001.fastq.gz
569363bdf0460defa1d4d78bf893cdf4  /data/putnamlab/KITT/hputnam/CHAIN/G029_S15_R1_001.fastq.gz
e34db2ed3f6e598dfed557232a86fdb4  /data/putnamlab/KITT/hputnam/CHAIN/G041_S9_R1_001.fastq.gz
3996b46dddf0742f5d3ac546d19d6577  /data/putnamlab/KITT/hputnam/CHAIN/G048_S2_R1_001.fastq.gz
3cde3f56722265129563ed06381c9fed  /data/putnamlab/KITT/hputnam/CHAIN/G053_S16_R1_001.fastq.gz
a15b0fc96d52ef78e36f6a8e97622af5  /data/putnamlab/KITT/hputnam/CHAIN/G056_S23_R1_001.fastq.gz
ffff1899d5b43613ce39a96945901f77  /data/putnamlab/KITT/hputnam/CHAIN/G060_S30_R1_001.fastq.gz
5f4b952085b11b04c79f9cb1069aca7e  /data/putnamlab/KITT/hputnam/CHAIN/G066_S10_R1_001.fastq.gz
41ddb02a786f8957c9e892b4e01de06b  /data/putnamlab/KITT/hputnam/CHAIN/G068_S3_R1_001.fastq.gz
fb537df24888a3fbf1d8e323efe3db12  /data/putnamlab/KITT/hputnam/CHAIN/G075_S35_R1_001.fastq.gz
13784bb7eea4bb858fd964c94c7cf4dc  /data/putnamlab/KITT/hputnam/CHAIN/G077_S24_R1_001.fastq.gz
dc465b43ad4c482256480b6fe73dc5fe  /data/putnamlab/KITT/hputnam/CHAIN/G079_S17_R1_001.fastq.gz
f3db83dfaab4725feef706d370a64801  /data/putnamlab/KITT/hputnam/CHAIN/G085_S33_R1_001.fastq.gz
fe1cf98ac44cde76dccacc52d03d5f21  /data/putnamlab/KITT/hputnam/CHAIN/G088_S4_R1_001.fastq.gz
80d0cae7fc8fe110e649aa051f699705  /data/putnamlab/KITT/hputnam/CHAIN/G091_S11_R1_001.fastq.gz
71db99ad7fc39bf9ca4f7e4769e31b3b  /data/putnamlab/KITT/hputnam/CHAIN/G094_S25_R1_001.fastq.gz
b7a7388487315fb7707e0c05dc0692b4  /data/putnamlab/KITT/hputnam/CHAIN/G096_S21_R1_001.fastq.gz
c1d615e9c5a22e813cec2a1fc3d9bbe5  /data/putnamlab/KITT/hputnam/CHAIN/R003_S26_R1_001.fastq.gz
b53d1afe6215aea62158abd33fa8bb94  /data/putnamlab/KITT/hputnam/CHAIN/R004_S12_R1_001.fastq.gz
c139c5fa609fdb24dd1a1f607f59c809  /data/putnamlab/KITT/hputnam/CHAIN/R009_S19_R1_001.fastq.gz
18dbe34ea12e104fcab2d2699e91de23  /data/putnamlab/KITT/hputnam/CHAIN/R015_S32_R1_001.fastq.gz
65970a9f945d125526a00ae7110c5246  /data/putnamlab/KITT/hputnam/CHAIN/R020_S5_R1_001.fastq.gz
e785879743a7f79a9860f5455129ee88  /data/putnamlab/KITT/hputnam/CHAIN/R021_S27_R1_001.fastq.gz
57182639c5b17230bcc68c8dfe705edf  /data/putnamlab/KITT/hputnam/CHAIN/R023_S6_R1_001.fastq.gz
438a6d9ea685c0589ac01f6cd32e8142  /data/putnamlab/KITT/hputnam/CHAIN/R025_S13_R1_001.fastq.gz
9f8beb1617ad70eb0a4f0ed5e00695ec  /data/putnamlab/KITT/hputnam/CHAIN/R026_S34_R1_001.fastq.gz
1b61c83bcf94445b110ec3f119886c9b  /data/putnamlab/KITT/hputnam/CHAIN/R038_S20_R1_001.fastq.gz
2467688126d7fd90e912ae29d882d132  /data/putnamlab/KITT/hputnam/CHAIN/R046_S31_R1_001.fastq.gz
01f46e6e0a52759ab09ae0b56a395aca  /data/putnamlab/KITT/hputnam/CHAIN/R048_S28_R1_001.fastq.gz
6da8abbfaee13cf209639eea78373698  /data/putnamlab/KITT/hputnam/CHAIN/R049_S14_R1_001.fastq.gz
cf983cc95ab66131d4398efa31e9cece  /data/putnamlab/KITT/hputnam/CHAIN/R059_S18_R1_001.fastq.gz
4908c57c39c93d1c46f2e9894029b839  /data/putnamlab/KITT/hputnam/CHAIN/R060_S7_R1_001.fastq.gz
```


c) Count number of reads per file 

check for code after @ in fastq.gz files(e.g.,@GWNJ).

```
zcat /data/putnamlab/KITT/hputnam/CHAIN/*.gz | echo $((`wc -l`/4)) > rawread.counts.txt

```



# 3) Run FastQC

a) Make folders for raw FastQC results and scripts

b) Write script for checking quality with FastQC and submit as job on bluewaves

```
nano /data/putnamlab/hputnam/CHAIN_RNASeq/scripts/fastqc_raw.sh
```

```  
#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/hputnam/CHAIN_RNASeq

module load all/FastQC/0.11.9-Java-11

for file in /data/putnamlab/KITT/hputnam/CHAIN/*.gz
do
fastqc $file --outdir /data/putnamlab/hputnam/CHAIN_RNASeq/qc
done
```

```
sbatch /data/putnamlab/hputnam/CHAIN_RNASeq/scripts/fastqc_raw.sh
```


c) Make sure all files were processed

```
ls -1 | wc -l 

70 files = 2 qc files per original 35 files
```

## Combined QC output into 1 file with MultiQC

```
module load MultiQC/1.9-intel-2020a-Python-3.8.2
multiqc .

```

c) Copy MultiQC files to local computer

```
scp -r hputnam@bluewaves.uri.edu:/data/putnamlab/hputnam/CHAIN_RNASeq/qc/multiqc_report.html /Users/hputnam/MyProjects/CHAIN_Transcriptomics
```

Samples RO20, RO23, RO60 have higher duplication, lower quality (still above X), and high sequence A content.


# 4) Trim and clean reads 

a) Make trimmed reads folder in all other results folders 

```

mkdir trimmed
cd trimmed

```

c) Write script for Trimming and run on bluewaves

```
nano /data/putnamlab/hputnam/CHAIN_RNASeq/scripts/trim.sh
```

```
#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/hputnam/CHAIN_RNASeq

module load fastp/0.19.7-foss-2018b

for file in "G021_S1" "G023_S8" "G029_S15" "G024_S22" "G028_S29" "G048_S2" "G041_S9" "G053_S16" "G056_S23" "G060_S30" "G068_S3" "G066_S10" "G079_S17" "G077_S24" "G075_S35" "G088_S4" "G091_S11" "G096_S21" "G094_S25" "G085_S33" "R020_S5" "R004_S12" "R009_S19" "R003_S26" "R015_S32" "R023_S6" "R025_S13" "R038_S20" "R021_S27" "R026_S34" "R060_S7" "R049_S14" "R059_S18" "R048_S28" "R046_S31"
do
fastp --in1 /data/putnamlab/KITT/hputnam/CHAIN/${file}_R1_001.fastq.gz  --out1 /data/putnamlab/hputnam/CHAIN_RNASeq/trimmed/${file}_R1_001.clean.fastq.gz --qualified_quality_phred 20 --length_required 50  
done
```
```
sbatch /data/putnamlab/hputnam/CHAIN_RNASeq/scripts/trim.sh
```


# 5) Check quality of trimmed files 

a) Check number of files 

```
ls -1 | wc -l
```


b) Run FastQC on trimmed data

mkdir trimmed_qc

```
nano /data/putnamlab/hputnam/CHAIN_RNASeq/scripts/fastqc_trimmed.sh
```

```  
#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=200GB
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/hputnam/CHAIN_RNASeq

module load all/FastQC/0.11.9-Java-11

for file in /data/putnamlab/hputnam/CHAIN_RNASeq/trimmed/*.gz
do
fastqc $file --outdir /data/putnamlab/hputnam/CHAIN_RNASeq/qc_trimmed
done
```

```
sbatch /data/putnamlab/hputnam/CHAIN_RNASeq/scripts/fastqc_trimmed.sh
```

## Combined QC output into 1 file with MultiQC

```
module load MultiQC/1.9-intel-2020a-Python-3.8.2
multiqc .

```

scp -r hputnam@bluewaves.uri.edu:/data/putnamlab/hputnam/CHAIN_RNASeq/qc_trimmed/multiqc_report.html /Users/hputnam/MyProjects/CHAIN_Transcriptomics

# 6) Align reads 

a) Generate genome build for the coral host

### Need to unzip genome files before running

```
gunzip /data/putnamlab/REFS/Pdam/GCF_003704095.1_ASM370409v1_genomic.fna.gz 
mv  /data/putnamlab/REFS/Pdam/GCF_003704095.1_ASM370409v1_genomic.fna /data/putnamlab/REFS/Pdam/GCF_003704095.1_ASM370409v1_genomic.fasta  
```
### HiSat2 Align reads to refernece genome
[HiSat2](https://daehwankimlab.github.io/hisat2/main/)
[HiSat2 Github](https://github.com/DaehwanKimLab/hisat2)

```
nano /data/putnamlab/hputnam/CHAIN_RNASeq/scripts/Hisat2_genome_build.sh
```

```
#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=200GB
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/hputnam/CHAIN_RNASeq

module load HISAT2/2.1.0-foss-2018b

hisat2-build -f /data/putnamlab/REFS/Pdam/GCF_003704095.1_ASM370409v1_genomic.fasta ./Pdam_ref
```

```
sbatch /data/putnamlab/hputnam/CHAIN_RNASeq/scripts/Hisat2_genome_build.sh
```

```
nano /data/putnamlab/hputnam/CHAIN_RNASeq/scripts/Hisat2_Sym_genome_build.sh
```

```
#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=200GB
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/hputnam/CHAIN_RNASeq

module load HISAT2/2.1.0-foss-2018b

hisat2-build -f /data/putnamlab/REFS/Sym_C1/SymbC1.Genome.Scaffolds.fasta ./C1_Sym_ref
```

```
sbatch /data/putnamlab/hputnam/CHAIN_RNASeq/scripts/Hisat2_Sym_genome_build.sh
```



b) Align reads to Coral genome

mkdir mapped

```
nano /data/putnamlab/hputnam/CHAIN_RNASeq/scripts/Hisat2_align2.sh
```

```
#!/bin/bash
#SBATCH -t 72:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --export=NONE
#SBATCH --mem=500GB
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/hputnam/CHAIN_RNASeq
#SBATCH -p putnamlab
#SBATCH --cpus-per-task=3

module load HISAT2/2.1.0-foss-2018b

for file in "G021_S1" "G023_S8" "G029_S15" "G024_S22" "G028_S29" "G048_S2" "G041_S9" "G053_S16" "G056_S23" "G060_S30" "G068_S3" "G066_S10" "G079_S17" "G077_S24" "G075_S35" "G088_S4" "G091_S11" "G096_S21" "G094_S25" "G085_S33" "R020_S5" "R004_S12" "R009_S19" "R003_S26" "R015_S32" "R023_S6" "R025_S13" "R038_S20" "R021_S27" "R026_S34" "R060_S7" "R049_S14" "R059_S18" "R048_S28" "R046_S31"
do
hisat2 -p 48 --dta -q -x /data/putnamlab/hputnam/CHAIN_RNASeq/Pdam_ref \
-U /data/putnamlab/hputnam/CHAIN_RNASeq/trimmed/${file}_R1_001.clean.fastq.gz \
-S /data/putnamlab/hputnam/CHAIN_RNASeq/mapped/${file}.sam
done
```

```
sbatch /data/putnamlab/hputnam/CHAIN_RNASeq/scripts/Hisat2_align2.sh
```

b) Align reads to Sym C1 genome

mkdir mapped

```
nano /data/putnamlab/hputnam/CHAIN_RNASeq/scripts/Hisat2_Sym_align2.sh
```

```
#!/bin/bash
#SBATCH -t 72:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --export=NONE
#SBATCH --mem=500GB
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/hputnam/CHAIN_RNASeq
#SBATCH -p putnamlab
#SBATCH --cpus-per-task=3

module load HISAT2/2.1.0-foss-2018b

for file in "G021_S1" "G023_S8" "G029_S15" "G024_S22" "G028_S29" "G048_S2" "G041_S9" "G053_S16" "G056_S23" "G060_S30" "G068_S3" "G066_S10" "G079_S17" "G077_S24" "G075_S35" "G088_S4" "G091_S11" "G096_S21" "G094_S25" "G085_S33" "R020_S5" "R004_S12" "R009_S19" "R003_S26" "R015_S32" "R023_S6" "R025_S13" "R038_S20" "R021_S27" "R026_S34" "R060_S7" "R049_S14" "R059_S18" "R048_S28" "R046_S31"
do
hisat2 -p 48 --dta -q -x /data/putnamlab/hputnam/CHAIN_RNASeq/C1_Sym_ref \
-U /data/putnamlab/hputnam/CHAIN_RNASeq/trimmed/${file}_R1_001.clean.fastq.gz \
-S /data/putnamlab/hputnam/CHAIN_RNASeq/Sym_mapped/${file}_Sym.sam
done
```

```
sbatch /data/putnamlab/hputnam/CHAIN_RNASeq/scripts/Hisat2_Sym_align2.sh
```

## Sort and convert sam to bam

### Coral

```
nano /data/putnamlab/hputnam/CHAIN_RNASeq/scripts/SAMtoBAM.sh
```


```
#!/bin/bash
#SBATCH -t 72:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --export=NONE
#SBATCH --mem=500GB
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/hputnam/CHAIN_RNASeq
#SBATCH -p putnamlab
#SBATCH --cpus-per-task=3

module load SAMtools/1.9-foss-2018b

for i in "G021_S1" "G023_S8" "G029_S15" "G024_S22" "G028_S29" "G048_S2" "G041_S9" "G053_S16" "G056_S23" "G060_S30" "G068_S3" "G066_S10" "G079_S17" "G077_S24" "G075_S35" "G088_S4" "G091_S11" "G096_S21" "G094_S25" "G085_S33" "R020_S5" "R004_S12" "R009_S19" "R003_S26" "R015_S32" "R023_S6" "R025_S13" "R038_S20" "R021_S27" "R026_S34" "R060_S7" "R049_S14" "R059_S18" "R048_S28" "R046_S31"
do
samtools view -b -S /data/putnamlab/hputnam/CHAIN_RNASeq/mapped/${i}.sam | samtools sort -o /data/putnamlab/hputnam/CHAIN_RNASeq/mapped/${i}.sorted.bam -O bam
done
```

```
sbatch /data/putnamlab/hputnam/CHAIN_RNASeq/scripts/SAMtoBAM.sh
```

### Sym

```
nano /data/putnamlab/hputnam/CHAIN_RNASeq/scripts/Sym_SAMtoBAM.sh
```


```
#!/bin/bash
#SBATCH -t 72:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --export=NONE
#SBATCH --mem=500GB
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/hputnam/CHAIN_RNASeq
#SBATCH --cpus-per-task=3

module load SAMtools/1.9-foss-2018b

for i in "G021_S1" "G023_S8" "G029_S15" "G024_S22" "G028_S29" "G048_S2" "G041_S9" "G053_S16" "G056_S23" "G060_S30" "G068_S3" "G066_S10" "G079_S17" "G077_S24" "G075_S35" "G088_S4" "G091_S11" "G096_S21" "G094_S25" "G085_S33" "R020_S5" "R004_S12" "R009_S19" "R003_S26" "R015_S32" "R023_S6" "R025_S13" "R038_S20" "R021_S27" "R026_S34" "R060_S7" "R049_S14" "R059_S18" "R048_S28" "R046_S31"
do
samtools view -b -S /data/putnamlab/hputnam/CHAIN_RNASeq/Sym_mapped/${i}_Sym.sam | samtools sort -o /data/putnamlab/hputnam/CHAIN_RNASeq/Sym_mapped/${i}.sorted.bam -O bam
done
```

```
sbatch /data/putnamlab/hputnam/CHAIN_RNASeq/scripts/Sym_SAMtoBAM.sh
```

### Remove Sam files to save space
```
rm /data/putnamlab/hputnam/CHAIN_RNASeq/mapped/*.sam
rm /data/putnamlab/hputnam/CHAIN_RNASeq/Sym_mapped/*.sam

```


# 7) Perform gene counts with stringTie

```
mkdir counts
cd counts

```


b) Assemble and estimate reads 

```
nano /data/putnamlab/hputnam/CHAIN_RNASeq/scripts/Coral_StringTie_Assemble.sh
```

```
#!/bin/bash
#SBATCH -t 72:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --export=NONE
#SBATCH --mem=500GB
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/hputnam/CHAIN_RNASeq
#SBATCH -p putnamlab
#SBATCH --cpus-per-task=3

module load StringTie/1.3.5-foss-2018b


for i in "G021_S1" "G023_S8" "G029_S15" "G024_S22" "G028_S29" "G048_S2" "G041_S9" "G053_S16" "G056_S23" "G060_S30" "G068_S3" "G066_S10" "G079_S17" "G077_S24" "G075_S35" "G088_S4" "G091_S11" "G096_S21" "G094_S25" "G085_S33" "R020_S5" "R004_S12" "R009_S19" "R003_S26" "R015_S32" "R023_S6" "R025_S13" "R038_S20" "R021_S27" "R026_S34" "R060_S7" "R049_S14" "R059_S18" "R048_S28" "R046_S31"
do
stringtie /data/putnamlab/hputnam/CHAIN_RNASeq/mapped/${i}.sorted.bam -p 48 -e -G /data/putnamlab/REFS/Pdam/GCF_003704095.1_ASM370409v1_genomic.gff -o /data/putnamlab/hputnam/CHAIN_RNASeq/mapped/${i}.gtf 
done
```

```
sbatch /data/putnamlab/hputnam/CHAIN_RNASeq/scripts/Coral_StringTie_Assemble.sh
```



c) Merge stringTie gtf results 

```
ls *gtf > mergelist.txt

```
 nano samplelist.txt
 
 C17	/data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/data/mapped/C17.gtf
 C18	/data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/data/mapped/C18.gtf
 C19	/data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/data/mapped/C19.gtf
 C20	/data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/data/mapped/C20.gtf
 C21	/data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/data/mapped/C21.gtf
 C22	/data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/data/mapped/C22.gtf
 C23	/data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/data/mapped/C23.gtf
 C24	/data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/data/mapped/C24.gtf
 C25	/data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/data/mapped/C25.gtf
 C26	/data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/data/mapped/C26.gtf
 C27	/data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/data/mapped/C27.gtf
 C28	/data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/data/mapped/C28.gtf
 C29	/data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/data/mapped/C29.gtf
 C30	/data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/data/mapped/C30.gtf
 C31	/data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/data/mapped/C31.gtf
 C32	/data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/data/mapped/C32.gtf
 E1		/data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/data/mapped/E1.gtf
 E2		/data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/data/mapped/E2.gtf
 E3		/data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/data/mapped/E3.gtf
 E4		/data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/data/mapped/E4.gtf
 E5		/data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/data/mapped/E5.gtf
 E6		/data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/data/mapped/E6.gtf
 E7		/data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/data/mapped/E7.gtf
 E8		/data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/data/mapped/E8.gtf
 E9		/data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/data/mapped/E9.gtf
 E10	/data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/data/mapped/E10.gtf  
 E11	/data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/data/mapped/E11.gtf     
 E12	/data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/data/mapped/E12.gtf    
 E13	/data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/data/mapped/E13.gtf    
 E14	/data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/data/mapped/E14.gtf    
 E15	/data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/data/mapped/E15.gtf  
 E16	/data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/data/mapped/E16.gtf    
   
f) Create gene matrix

```
nano /data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/scripts/GTFtoCounts.sh
```

```
#!/bin/bash
#SBATCH -t 72:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --export=NONE
#SBATCH --mem=500GB
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/data/raw
#SBATCH -p putnamlab
#SBATCH --cpus-per-task=3

module load StringTie/1.3.5-foss-2018b
module load Python/2.7.15-foss-2018b


python prepDE.py  -i samplelist.txt -g Poc_gene_count_matrix.csv
```

```
sbatch /data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/scripts/GTFtoCounts.sh
```


g) Secure-copy gene counts onto local computer

```
scp jillashey@bluewaves.uri.edu:/data/putnamlab/jillashey/Francois_data/Florida/stringTie/Ofav/GTF_merge/gene_count_ofav_matrix.csv /Users/jillashey/Desktop/Putnamlab/Repositories/SedimentStress/SedimentStress/Output/DESeq2/star/
```

