#!/bin/bash -l
#SBATCH --job-name=fsjtrim
#SBATCH --output=fsjtrim-%A_%a.out
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=20000
#SBATCH --time=4:00:00
#SBATCH --partition=short
#SBATCH --account=bscb02
#SBATCH --array=296-357

# Run trimmomatic on all samples for kneaddata later on

# get your array number from the fsj_wgs_sample_list_2022.txt
# may need to exclude=cbsubscb02
# date
d1=$(date +%s)

newdir=${SLURM_JOB_ID}

mySample=`sed -n "${SLURM_ARRAY_TASK_ID}p" fsj_wgs_sample_list_2022.txt | cut -f 1`
myRun=`sed -n "${SLURM_ARRAY_TASK_ID}p" fsj_wgs_sample_list_2022.txt | cut -f 2`
myOutStr=${mySample}_${myRun}
myOutDir=/fs/cbsuclarkfs1/storage/tn337/Kraken/trimmed_seq

echo $HOSTNAME
echo sample name
echo $mySample
echo run ID
echo $myRun
echo output directory
echo $myOutDir
echo output file prefix
echo $myOutStr
echo task index
echo $SLURM_ARRAY_TASK_ID
echo $newdir

mkdir -p /workdir/$USER/$newdir
cd /workdir/$USER/$newdir

/programs/bin/labutils/mount_server cbsuclarkfs1 /storage
#cp /fs/cbsuclarkfs1/storage/ejc87/fsj/tram_wgs/rawdata/${mySample}_C*${myRun}*.fq.gz . ## the first 1-295 samples
cp /fs/cbsuclarkfs1/storage/tn337/fsj_wgs20x/rawdata/batch3/usftp21.novogene.com/raw_data/${mySample}/*.fq.gz . ## the last batch of samples 296-357
cp /fs/cbsuclarkfs1/storage/tn337/fsj_wgs20x/pre-processing/2022-Jan/NEB_Illumina_adapters.fa .
cp /fs/cbsuclarkfs1/storage/tn337/fsj_wgs20x/pre-processing/2021-04-14/FSJv2_reference/*.* .

#mkdir -p $myOutDir

pwd

## trimmomatic
echo running trimmomatic
java -jar /programs/trimmomatic/trimmomatic-0.36.jar PE -threads ${SLURM_NTASKS} ./*_1.fq.gz ./*_2.fq.gz -baseout ${myOutStr}.fq.gz ILLUMINACLIP:NEB_Illumina_adapters.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:15
echo finish trimmomatic


## copy output to storage
mv *1U.fq.gz $myOutDir
mv *2U.fq.gz $myOutDir
mv *1P.fq.gz $myOutDir
mv *2P.fq.gz $myOutDir

cd ..
rm -r ./$newdir

# date
d2=$(date +%s)
sec=$(( ( $d2 - $d1 ) ))
hour=$(echo - | awk '{ print '$sec'/3600}')
echo Runtime: $hour hours \($sec\s\)
