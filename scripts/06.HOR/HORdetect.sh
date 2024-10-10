#!/bin/bash

# 初始化变量
input_csv=""
output_bed="./dimer_sf_HOR.bed"
output_fasta="./dimer_sf_sequences.fa"
mkdir -p results
mkdir -p data
output="./results"
skip_timeout=false  # 跳过超时的开关，默认关闭
timeout_duration=60  # 设置超时时间为60秒

# 使用 getopts 解析命令行参数
while getopts ":i:b:f:s" opt; do
  case ${opt} in
    i )
      input_csv=$OPTARG  # -i 后的参数
      ;;
    b )
      output_bed=$OPTARG  # -b 后的参数
      ;;
    f )
      output_fasta=$OPTARG  # -f 后的参数
      ;;
    s )
      skip_timeout=true  # -s 开关用于跳过超时
      ;;
    \? )
      echo "Invalid option: -$OPTARG" 1>&2
      exit 1
      ;;
    : )
      echo "Option -$OPTARG requires an argument." 1>&2
      exit 1
      ;;
  esac
done
echo "Input CSV: $input_csv"
echo "Output BED: $output_bed"
echo "Output FASTA: $output_fasta"

if [ -z "$input_csv" ]; then
  echo "Usage: $0 
  -i <input_csv> 
  [-o1 <output_bed>]
  [-o2 <output_fasta>]
  [-s] if only want to get De Bruijn plots , please use this.Because too much data may cause very low speed for running."
  exit 1
fi


# for example:  Rscript 06.HOR.r -i data/Zpal_CentZp_subfamily.csv -o1 ./data/Zpal_dimer_sf_HOR.bed -o2 ./data/Zpal.fa

Rscript 06.HOR.r -i "$input_csv" -o1 "$output_bed" -o2 "$output_fasta"

grep ">" "$output_fasta" | sed 's/>//' > ./data/chr.list

conda activate HOR_env
while read chr; do
  echo "Processing $chr..."

  if $skip_timeout; then
    # 如果启用超时功能，使用 timeout 限制每次运行时间为 60 秒
    timeout $timeout_duration python HORdetect.py -fi "$output_fasta" -pos "$output_bed" -p "$chr" -O $output
    exit_status=$?

    if [ $exit_status -eq 124 ]; then
      echo "Timeout occurred for $chr, skipping to the next..."
      continue  # 跳到下一个循环
    fi
  else
    # 如果没有启用超时功能，直接执行命令
    python HORdetect.py -fi "$output_fasta" -pos "$output_bed" -p "$chr" -O $output
  fi

done < ./data/chr.list


