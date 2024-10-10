#!/bin/bash

date

#source /public/home/cszx_huangxh/biosoftwares/Anaconda3_8/etc/profile.d/conda.sh
#conda activate /public/home/cszx_huangxh/biosoftwares/Anaconda3_8/envs/gene_anno
#conda env list
source ~/.bashrc
data_dir="${gulei}/genome_data"
fa_files=("Nip_genome.fa" "Zpal_genome.fa")
#fa_files=("Zlat_genome.fa" "Zpal_genome.fa")
for fa_file in "${fa_files[@]}"; do
  species_name="${fa_file%%_*}"
  sbatch <<-EOF
#!/bin/bash
#SBATCH --job-name=EDTA_${species_name}
#SBATCH --output=$gulei/cen_analysis/EDTA/${species_name}.log
#SBATCH --error=$gulei/cen_analysis/EDTA/${species_name}.log
#SBATCH --time=30-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=168G

source ~/.bashrc
expmini
conda activate EDTA2.2

cd $gulei/cen_analysis/EDTA
mkdir -p $species_name && cd $species_name
cp -rf $gulei/biosoftwares/EDTA/* .
ln -sf $data_dir/${fa_file} .

if [ "$species_name" == "Nip" ]; then
	perl EDTA.pl --genome $fa_file --species rice --step all --overwrite 1 --sensitive 1 --anno 1 --threads 24 --force 1 --evaluate 1
else 
	perl EDTA.pl --genome $fa_file --species others --step all --overwrite 1 --sensitive 1 --anno 1 --threads 24 --force 1 --evaluate 1
fi

date
                                                 
EOF
done
