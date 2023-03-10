## Run Braken
## Nov 2022
source $HOME/miniconda3/bin/activate
conda activate kraken2

export PATH=/programs/kraken2.1.0:$PATH
export PATH=/programs/Bracken-2.0:/programs/Bracken-2.0/src:/programs/Bracken-2.0/analysis_scripts:$PATH

#This pipeline is for analyzing Kraken output results. You can access the instructions here https://github.com/jenniferlu717/Bracken

KRAKEN_DB=/workdir/tn337/Kraken/kraken_db_2021
cd /workdir/tn337/Kraken/bracken

'''Step 1: Generate the Bracken database file'''

bracken-build -d ${KRAKEN_DB} -t 24 -k 35 -l 150

#        `${KRAKEN_DB}`  = location of the built Kraken 1.0 or Kraken 2.0 database
#        `${THREADS}`    = number of threads to use with Kraken and the Bracken scripts
#        `${KMER_LEN}`   = length of kmer used to build the Kraken database
#                                Kraken 1.0 default kmer length = 31
#                                Kraken 2.0 default kmer length = 35
#                                Default set in the script is 35.
#        `${READ_LEN}`   = the read length of your data
#                                e.g., if you are using 100 bp reads, set it to `100`.

### LOG
[tn337@cbsulm11 bracken]$ bracken-build -d ${KRAKEN_DB} -t 24 -k 35 -l 150
 >> Selected Options:
       kmer length = 35
       read length = 150
       database    = /workdir/tn337/Kraken/kraken_db_2021
       threads     = 24
 >> Checking for Valid Options...
 >> Creating database.kraken [if not found]
          database.kraken exists, skipping creation....
          Finished creating database.kraken [in DB folder]
 >> Creating database150mers.kmer_distrib
        >>STEP 0: PARSING COMMAND LINE ARGUMENTS
                Taxonomy nodes file: /workdir/tn337/Kraken/kraken_db_2021/taxonomy/nodes.dmp
                Seqid file:          /workdir/tn337/Kraken/kraken_db_2021/seqid2taxid.map
                Num Threads:         24
                Kmer Length:         35
                Read Length:         150
        >>STEP 1: READING SEQID2TAXID MAP
                243142 total sequences read
        >>STEP 2: READING NODES.DMP FILE
                2312166 total nodes read
        >>STEP 3: READING DATABASE.KRAKEN FILE
                105660 total sequences read
        >>STEP 4: CONVERTING KMER MAPPINGS INTO READ CLASSIFICATIONS:
                150mers, with a database built using 35mers
                105660 sequences converted...
        Time Elaped: 62 minutes, 54 seconds, 0.00000 microseconds
        =============================
PROGRAM START TIME: 11-23-2022 19:12:41
...21283 total genomes read from kraken output file
...creating kmer counts file -- lists the number of kmers of each classification per genome
...creating kmer distribution file -- lists genomes and kmer counts contributing to each genome
PROGRAM END TIME: 11-23-2022 19:12:43
          Finished creating database150mers.kraken and database150mers.kmer_distrib [in DB folder]
          *NOTE: to create read distribution files for multiple read lengths,
                 rerun this script specifying the same database but a different read lengths



''' Run Bracken for Abundance Estimation '''
#Given the expected kmer distribution for genomes in a kraken database along with a kraken report file, the number of reads belonging to each species (or genus) is estimated using the estimate_abundance.py file, run with the following command line:


KRAKEN_DB=/workdir/tn337/Kraken/kraken_db_2021
SAMPLE=/fs/cbsuclarkfs1/storage/tn337/Kraken/kraken_results/2022/JDSP126_stdb.rpt
LEVEL=P

bracken -d ${KRAKEN_DB} -i ${SAMPLE}.rpt -o ${SAMPLE}.bracken -r ${READ_LEN} -l ${LEVEL} -t 10
## GIVES ERROR:
#[tn337@cbsulm11 bracken]$ bracken -d ${KRAKEN_DB} -i ${SAMPLE}.rpt -o ${SAMPLE}.bracken -r ${READ_LEN} -l ${LEVEL} -t 10
# >> Checking for Valid Options...
# ERROR: /workdir/tn337/Kraken/kraken_db_2021/database-lmers.kmer_distrib does not exist
#        Run bracken-build to generate the kmer distribution file.

https://github.com/jenniferlu717/Bracken/issues/109

''' TRY RUNNING EST ABUNDNCE DIRECTLY AS SUGGESTED -- NOV 28, 2022'''

cd /workdir/tn337/Kraken/bracken
git clone https://github.com/jenniferlu717/Bracken.git
cd Bracken
bash install_bracken.sh

cd /workdir/tn337/Kraken/bracken

for CTAX in `cat taxlevel.txt`; do  ## run bracken at 6 taxonomic levels: P,C,O,F,G,S # levels listed in taxlevel.txt
  for cfile in `cat full_unique_sample_list_2022.txt`; do  # list of sample names

      SAMPLE=$cfile
      KRAKEN_DB=/workdir/tn337/Kraken/kraken_db_2021
      OUTDIR=/fs/cbsuclarkfs1/storage/tn337/Kraken/Bracken/2022
      INPATH=/fs/cbsuclarkfs1/storage/tn337/Kraken/kraken_results/2022

      python /workdir/tn337/Kraken/bracken/Bracken/src/est_abundance.py -i ${INPATH}/${SAMPLE}_stdb.rpt -k ${KRAKEN_DB}/database150mers.kmer_distrib -l ${CTAX} -o ${OUTDIR}/${CTAX}/${SAMPLE}.bracken -t 10

    done
    echo $CTAX
    date
done
date
