#!/bin/bash

source ~/.bashrc
gulei=/public/home/cszx_huangxh/qiujie/collabrators/gulei
data_dir="${gulei}/genome_data"
output_dir="${gulei}/cen_analysis/cen_define_3.0"
data_files=("Zpal_genome.fa" "Nip_genome.fa")

for data_file in "${data_files[@]}"; do
  data_name="${data_file%%_*}"
  sbatch <<-EOF
#!/bin/bash
#SBATCH --job-name=Cen_${data_name}
#SBATCH --output=$gulei/cen_analysis/cen_define_3.0/${data_name}.log
#SBATCH --error=$gulei/cen_analysis/cen_define_3.0/${data_name}.log
#SBATCH --time=30-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=256G
#SBATCH --partition=CPU

source ~/.bashrc
expana 
conda activate gene_anno

cd "${gulei}/cen_analysis/cen_define_3.0"
echo "进入cen_define_3.0目录"
mkdir -p "${data_name}"
cd "${data_name}"


cp -rf $gulei/biosoftwares/CentIER/* .
ln -sf $data_dir/${data_file} .
./centIER.py ${data_file} 

conda deactivate

EOF
done
