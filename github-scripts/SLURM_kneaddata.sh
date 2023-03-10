#!/bin/bash -l
#SBATCH --job-name=kneaddata
#SBATCH --output=kneaddata-%A_%a.out
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=100000
#SBATCH --time=1-10:00:00
#SBATCH --partition=long7,long30
#SBATCH --exclude=cbsubscb18
#SBATCH --account=bscb02
#SBATCH --array=35,166,175,178,199-201,206-209

# run Kneaddata on all 291 samples
# get your array number from the fsj_wgs_sample_list_2022.txt

# date
d1=$(date +%s)

newdir=${SLURM_JOB_ID}

allSamples=( `cat full_unique_sample_list_2022.txt` )
mySample=${allSamples[${SLURM_ARRAY_TASK_ID}-1]}

myOutDir=/fs/cbsuclarkfs1/storage/tn337/Kraken/kneaddata2022

echo $HOSTNAME
echo sample name
echo $mySample
echo output directory
echo $myOutDir
echo task index
echo $SLURM_ARRAY_TASK_ID
echo $newdir

mkdir -p /workdir/$USER/$newdir
cd /workdir/$USER/$newdir

echo running analyses in directory:
pwd

## navigate into the miniconda3 environment "kraken2" with kneaddata program installed
echo navigating into miniconda3 kraken2 environment
source $HOME/miniconda3/bin/activate
source activate kraken2

## copy required files
echo copying files
/programs/bin/labutils/mount_server cbsuclarkfs1 /storage
cp /fs/cbsuclarkfs1/storage/tn337/Kraken/trimmed_seq/merged/${mySample}_merged_1P.fq.gz . #copy all the 1P and 2P paired trimmed seq
cp /fs/cbsuclarkfs1/storage/tn337/Kraken/trimmed_seq/merged/${mySample}_merged_2P.fq.gz .
cp -r /fs/cbsuclarkfs1/storage/tn337/Kraken/FSJgenome_db_2022/ . #copy database


export PATH=/programs/bowtie2-2.4.5-linux-x86_64:$PATH

## kneaddata
echo running kneaddata
kneaddata -i ${mySample}_merged_1P.fq.gz -i ${mySample}_merged_2P.fq.gz -o ./ -db FSJgenome_db_2022/ -t 8 --bowtie2-options "--very-sensitive --dovetail" --bypass-trim --remove-intermediate-output
echo finished kneaddata

echo deactivating from conda environment kraken2
conda deactivate

## copy output to storage
echo copying files to outDir

cp ${mySample}*.fastq $myOutDir
cp *_kneaddata.log $myOutDir

## clean up
cd ..
rm -r ./$newdir

# date
d2=$(date +%s)
sec=$(( ( $d2 - $d1 ) ))
hour=$(echo - | awk '{ print '$sec'/3600}')
echo Runtime: $hour hours \($sec\s\)
