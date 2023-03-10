# FSJ_Blood_metagenomics
Tram N. Nguyen
<br />
Updated: Dec 15, 2022
<br />

This repository contains scripts to isolate and identify microbes within host whole-blood samples.
<br />
Briefly, we start with Illumina paired-end sequences. We first trim adapter sequences and apply read quality filters (ILLUMINACLIP:NEB_Illumina_adapters.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:15). Please see the [Trimmomatic manual](http://www.usadellab.org/cms/?page=trimmomatic) for details. We then align these reads to our reference genome and discard host reads, retaining unmapped reads with Kneaddata. We then build the standard Kraken database and run the program to identify microbes in the remaining reads. A customized database can be created for your specific study following these instructions: https://ccb.jhu.edu/software/kraken/MANUAL.html#custom-databases. Next, compute the abundance of species in our metagenomic samples using the program Bracken, then do inital exploratory plots and visualize the data. Then, use SourceTracker to pinpoint potential lab or sample preparation contamination sources. Finally, remove likely contaminants with the program Decontam.

<br /> Please feel free to contact tn337@cornell.edu for any questions.


# Dependencies
## Software
- [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) v0.36
- [bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml#:~:text=Bowtie%202%20is%20an%20ultrafast,long%20(e.g.%20mammalian)%20genomes.) v2.4.5
- [kneaddata](https://huttenhower.sph.harvard.edu/kneaddata/) v0.12.0
- [Kraken2](https://ccb.jhu.edu/software/kraken2/) v2.1.0
- [Bracken](https://ccb.jhu.edu/software/bracken/) v2.0
- [SourceTracker](https://github.com/biota/sourcetracker2) v1.0.1
- [phyloseq](https://joey711.github.io/phyloseq/) v
- [Decontam](https://benjjneb.github.io/decontam/vignettes/decontam_intro.html) 

# Step-by-step pipeline

### 1. Run Trimmomatic and apply read quality filters. 
This is executed through a SLURM scheduler, with each job as one sample. 

Required inputs:
- A text file containing names of your samples
- Your reference genome
- An adapter file containing sequences you want to trim

Run ```SLURM_trimmomatic.sh```


<br />


### 2. Run Kneaddata. 
We will align our reads to our host reference genome and discard these, keeping only the remaining reads. 

We first need to build our reference database.
```bowtie2-build /workdir/tn337/final.assembly.fasta FSJ_genome_db```

Then we run KneadData.
KneadData is a tool designed to perform quality control on metagenomic and metatranscriptomic sequencing data, especially data from microbiome experiments. See http://huttenhower.sph.harvard.edu/kneaddata

Because we already trimmed reads, we use --bypass-trim option in KneadData. The input for the script are your R1_paired_output and R2_paired_output that are trimmed of adapters. This step is implemented with SLURM and done per sample.

Run ```SLURM_kneaddata.sh```

<br />

### 3. Build Kraken database
Because Kraken database downloads files from several sources such as NCBI, if any changes are made to these files on their end, the build will have trouble completing. Thus, this portion of the workflow will take some time to troubleshoot and the script provided serves as more of a documentation of all the different troubleshooting steps I've had to try recently. I have included several github issue pages for common problems, but this step will constant need to be updated and tweaked. The latest version that successfully completed was August 2022 for this project.

Run ```kraken-build-commands.sh```

<br />

### 4. Run Kraken!
Now that we've hopefully built our microbe/pathogen database, let's run Kraken to assign taxonomic labels to those remaining reads to quantify microbes within the host whole-blood samples. Below is code to loop through our samples but this can also be parallelized as separate cluster runs. You will need your reads from the KneadData output and a text file containing your sample IDs.


```
source $HOME/miniconda3/bin/activate
source activate kraken2 #navigate into kraken2 conda environment

export PATH=/programs/bowtie2-2.4.5-linux-x86_64:$PATH
DBNAME=/fs/cbsuclarkfs1/storage/tn337/Kraken/kraken_db_2021

cd /workdir/tn337/Kraken/

for number in $(seq 1 294);
   do
         mySample=`sed -n "${number}p" ./full_unique_sample_list_2022.txt | cut -f 1` #look through a list of samples
         echo ${mySample}
         kraken2 --db ${DBNAME} --output /kraken_results/2022/${mySample}_stdb.krkn \
         --report /kraken_results/2022/${mySample}_stdb.rpt \
         --paired /kneaddata2022/${mySample}_merged_1P_kneaddata_paired_1.fastq \
         /kneaddata2022/${mySample}_merged_1P_kneaddata_paired_2.fastq
   done
```

<br />

### 5. Run Bracken
Because Kraken classifies reads to the best matching location in the taxonomic tree, but does not estimate abundances of species, we will now use Bracken to compute the abundance of species in DNA sequences from a metagenomics sample. You will once again need to build a Bracken database and use your Kraken outputs from the previous step.

Run ```run_Bracken2022.sh```

<br />

### 6. Exploratory Visualization and QC on Bracken outputs
This script provides code for creating several exploratory plots from Bracken results. For example, this code can visualize the 50 most prevalent species within your data. It can also be used to explore potential lab/handling contamination by including metadata from collectors and sequencing batches. These scripts are explorative and can be a start to doing some initial QC and analyses in R. Analyses are in ```plot_bracken_2022.R```

<br />

### 7. Run SourceTracker to identify more potential sources of lab/handling contamination
To ensure that the microbes we are detecting are truly host-specific and not a part of laboratory preparation or sampling handling, we will download metagenomes from the Human Microbiome Project to serve as a negative control for our dataset. Because these samples were opportunistically sampled, we could not produce a reliable negative blank control. Therefore, this strategy may also be helpful for trickier samples where a control was not possible.\

We will be using the program SourceTracker, which uses Bayesian approach to predict the source of microbial communities in a set of input samples (i.e., the sink samples). We want to know the proportion of microbes that are likely from the HMP dataset and thus, likely human-associated and not bird-associated.

In this study, we download mWGS from relevant body sites (in our case: nasal, mouth, skin) from the Human Microbiome Project.
![HMP logo](https://images.squarespace-cdn.com/content/538e5c3ce4b02add9dca1fd5/1410527133814-0BSTCQ4YHQ9F9O1TNOMN/?content-type=image%2Fpng)


