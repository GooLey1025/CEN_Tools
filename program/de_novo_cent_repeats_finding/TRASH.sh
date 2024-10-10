#!/bin/bash

source ~/.bashrc

data_dir="$gulei/genome_data"
data_files=("Nip_genome.fa" "Zpal_genome.fa")
for data_file in ${data_files[@]}; do
  data_name="${data_file%%_*}"
  prefix=$data_name

  sbatch --export=prefix=$data_name,data_name=$data_name,data_file=$data_file,data_dir=$data_dir <<-EOF
#!/bin/bash
#SBATCH --job-name=TRAS_${data_name}
#SBATCH --output=/public/home/cszx_huangxh/qiujie/collabrators/gulei/cen_analysis/CenTools/cen/${data_name}.log
#SBATCH --error=/public/home/cszx_huangxh/qiujie/collabrators/gulei/cen_analysis/CenTools/cen/${data_name}.log
#SBATCH --time=30-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=128G
#SBATCH --partition=CPU


#####  1\\  de novo TRASH repeat finding
source ~/.bashrc
expana
source /public/home/cszx_huangxh/biosoftwares/Anaconda3_8/etc/profile.d/conda.sh
conda activate trash
        
cd /public/home/cszx_huangxh/qiujie/collabrators/gulei/cen_analysis/CenTools/cen
mkdir -p ${data_name}/${prefix}_repeat
${gulei}/biosoftwares/TRASH/TRASH_run.sh $data_dir/${data_file} --simpleplot --horclass CEN145 --frep 100 --par 24 --o /public/home/cszx_huangxh/qiujie/collabrators/gulei/cen_analysis/CenTools/cen/${data_name}/${prefix}_repeat


#####  2\\  get each_centromere_consensus for accessions
python consensus.py -i ./${data_name}/${prefix}_repeat/Summary.of.repetitive.regions.${data_file}.csv -o ./${data_name}/TRASH_${prefix}_consensus.txt


#####  3\\  最大可能切割潜在的HOR,以获得尽可能多的monomer
mkdir -p ./${data_name}/${prefix}_HOR 
${gulei}/biosoftwares/TRASH/TRASH_run.sh $data_dir/${data_file} --seqt ./${data_name}/TRASH_${prefix}_consensus.txt --simpleplot --horclass `less ./${data_name}/TRASH_${prefix}_consensus.txt | cut -d "," -f 1 | sed -n "2p"` --par 24 --o /public/home/cszx_huangxh/qiujie/collabrators/gulei/cen_analysis/CenTools/cen/${data_name}/${prefix}_HOR  --N.max.div 5

EOF
done


