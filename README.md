# FSJ_Blood_metagenomics
Tram N. Nguyen
<br />
Updated: Dec 15, 2022
<br />

This repository contains scripts to isolate and identify microbes within host whole-blood samples.
<br />
Briefly, we start with Illumina paired-end sequences. We first trim adapter sequences and apply read quality filters (ILLUMINACLIP:NEB_Illumina_adapters.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:15). Please see the [Trimmomatic manual](http://www.usadellab.org/cms/?page=trimmomatic) for details. We then align these reads to our reference genome and discard host reads, retaining unmapped reads with Kneaddata. We then build the standard Kraken database and run the program to identify microbes in the remaining reads. A customized database can be created for your specific study following these instructions: https://ccb.jhu.edu/software/kraken/MANUAL.html#custom-databases. 

<br /> Please feel free to contact tn337@cornell.edu for any questions.


# Dependencies
## Software
- [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) v0.36
- [bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml#:~:text=Bowtie%202%20is%20an%20ultrafast,long%20(e.g.%20mammalian)%20genomes.) v2.4.5
- [kneaddata](https://huttenhower.sph.harvard.edu/kneaddata/) v0.12.0
- [Kraken2](https://ccb.jhu.edu/software/kraken2/) v2.1.0

# Step-by-step pipeline

### 1. Run Trimmomatic and apply read quality filters. 
This is executed through a SLURM scheduler, with each job as one sample. 

Required inputs:
- A text file containing names of your samples
- Your reference genome
- An adapter file containing sequences you want to trim

Run ```SLURM_trimmomatic.sh```

### 2. Run Kneaddata. 
We will align our reads to our host reference genome and discard these, keeping only the remaining reads. 

We first need to build our reference database.
```bowtie2-build /workdir/tn337/final.assembly.fasta FSJ_genome_db```

Then we run KneadData.
KneadData is a tool designed to perform quality control on metagenomic and metatranscriptomic sequencing data, especially data from microbiome experiments. See http://huttenhower.sph.harvard.edu/kneaddata

Because we already trimmed reads, we use --bypass-trim option in KneadData. The input for the script are your R1_paired_output and R2_paired_output that are trimmed of adapters. This step is implemented with SLURM and done per sample.

Run ```SLURM_kneaddata.sh```

### 3. Build Kraken database
Because Kraken database downloads files from several sources such as NCBI, if any changes are made to these files on their end, the build will have trouble completing. Thus, this portion of the workflow will take some time to troubleshoot and the script provided serves as more of a documentation of all the different troubleshooting steps I've had to try recently. I have included several github issue pages for common problems, but this step will constant need to be updated and tweaked. The latest version that successfully completed was August 2022 for this project.

Run ```kraken-build-commands.sh```


### 4. Run Kraken!
Now that we've hopefully built our microbe/pathogen database, let's run Kraken to assign taxonomic labels to those remaining reads to quantify microbes within the host whole-blood samples.

Kraken_pipeline_2022.sh) outputs: /fs/cbsuclarkfs1/storage/tn337/Kraken/kraken_results/2022/
4. Run Braken (run_Bracken2022.sh) outputs:/fs/cbsuclarkfs1/storage/tn337/Kraken/Bracken/2022
5. Plot a bunch of contamination/bracken QC and analyses on lm11 R.
6. Begin Qiime2 protocol -- FIRST, kraken-biom (kraken2 reports to .biom table)

