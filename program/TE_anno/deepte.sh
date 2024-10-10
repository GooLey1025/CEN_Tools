#!/bin/bash
source ~/.bashrc

data_dir="${gulei}/genome_data"
fa_files=("Nip_genome.fa")
for fa_file in "${fa_files[@]}"; do
  species_name="${fa_file%%_*}"
  sbatch <<-EOT
#!/bin/bash
#SBATCH --job-name=DpTe_${species_name}
#SBATCH --output=$gulei/cen_analysis/DeepTE/DeepTE_${species_name}.log
#SBATCH --error=$gulei/cen_analysis/DeepTE/DeepTE_${species_name}.log
#SBATCH --time=30-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=128G
#SBATCH --partition=GPU
#SBATCH --gres=gpu:1
date

source ~/.bashrc
#source /public/home/cszx_huangxh/biosoftwares/Miniconda/install/etc/profile.d/conda.sh
expmini
conda activate DeepTE

cd /public/home/cszx_huangxh/qiujie/collabrators/gulei/cen_analysis/DeepTE
mkdir -p ${species_name} && cd ${species_name}

cp -rf /public/home/cszx_huangxh/qiujie/collabrators/gulei/cen_analysis/EDTA/${species_name}/* .

grep "LTR/unknown" ${species_name}_genome.fa.mod.EDTA.TElib.fa | seqtk subseq ${species_name}_genome.fa.mod.EDTA.TElib.fa - >LTR_unknown.fa || exit 1

grep -v "LTR/unknown" ${species_name}_genome.fa.mod.EDTA.TElib.fa | seqtk subseq ${species_name}_genome.fa.mod.EDTA.TElib.fa - >LTR_known.fa || exit 1


python /public/home/cszx_huangxh/qiujie/collabrators/gulei/biosoftwares/DeepTE/DeepTE.py -i LTR_unknown.fa -sp P -m_dir /public/home/cszx_huangxh/qiujie/collabrators/gulei/biosoftwares/DeepTE/Plants_model -fam LTR 

sed 's/LTR\/unknown__ClassI_LTR_Copia/LTR\/Copia/' opt_DeepTE.fasta | sed 's/LTR\/unknown__ClassI_LTR_Gypsy/LTR\/Gypsy/' | sed 's/LTR\/unknown__ClassI_LTR/LTR\/unknown/' > LTR_unknown_DeepTE.fa

cat LTR_known.fa LTR_unknown_DeepTE.fa > ${species_name}_genome.fa.mod.EDTA.TElib.fa

conda activate EDTA2.2

if [ "$species_name" == "Nip" ]; then
    perl EDTA.pl --genome $fa_file --species Rice --step all --overwrite 1  --anno 1 --threads 24 --force 1 --evaluate 1 
else
    perl EDTA.pl --genome $fa_file --species others --step all --overwrite 1  --anno 1 --threads 24 --force 1 --evaluate 1            
fi

date
EOT
done
