if (!require("argparse")) {
    install.packages("argparse")  # 安装 argparse 包
    library(argparse)  # 加载 argparse 包
} else {
    library(argparse)  # 如果已安装，直接加载
}
if (!require("dplyr")) {
    install.packages("dplyr") 
    library(dplyr)  
} else {
    library(dplyr)  
}

parser <- ArgumentParser(description = 'Process input and output files')

parser$add_argument('-i', '--input_csv', required=TRUE, help='Input subfamily CSV file path')
parser$add_argument('-o1', '--output_bed', default="./dimer_sf_HOR.bed", help='Output BED file path')
parser$add_argument('-o2', '--output_fasta', default="./dimer_sf_sequences.fa", help='Output FASTA file path')
args <- parser$parse_args()

csv_file <- args$input_csv
output_bed_file <- args$output_bed
output_fasta_file <- args$output_fasta


temp_cen_sf <- read.csv(csv_file,sep = "\t" )
temp_cen_sf <- temp_cen_sf %>%
  arrange(Chr, start)  # 按染色体和起始位置排序

# 初始化变量，用于存储最终结果
result <- data.frame()

# 遍历每条染色体
for (chr in unique(temp_cen_sf$Chr)) {
  
  # 获取当前染色体的数据
  chr_data <- temp_cen_sf %>% filter(Chr == chr)
  
  # 初始化变量用于合并同样的 SF 值
  current_SF <- chr_data$SF[1]  # 初始 SF 值
  start_pos <- chr_data$start[1]  # 起始位置
  end_pos <- chr_data$end[1]  # 结束位置
  count <- 1  # 连续 SF 的计数
  
  # 遍历当前染色体的数据，从第二行开始
  for (i in 2:nrow(chr_data)) {
    # 判断当前行的 SF 是否与上一行相同，以及 start 和 end 的差值是否小于等于 50
    if (chr_data$SF[i] == current_SF && abs(chr_data$start[i] - end_pos) <= 50) {
      # 如果 SF 值相同且差值在 50 以内，继续合并，更新结束位置和计数
      end_pos <- chr_data$end[i]
      count <- count + 1
    } else {
      # 如果不满足条件，将之前的单元记录到结果中
      result <- rbind(result, data.frame(
        Chr = chr,
        start = start_pos,
        end = end_pos,
        strand = chr_data$strand[i - 1],
        SF = current_SF,
        unit_count = count  # 新增一列记录单元数
      ))
      
      # 开始处理新的 SF 值
      current_SF <- chr_data$SF[i]
      start_pos <- chr_data$start[i]
      end_pos <- chr_data$end[i]
      count <- 1
    }
  }
  
  # 将最后一个单元记录到结果中
  result <- rbind(result, data.frame(
    Chr = chr,
    start = start_pos,
    end = end_pos,
    strand = chr_data$strand[nrow(chr_data)],
    SF = current_SF,
    unit_count = count  # 新增一列记录单元数
  ))
}

# 查看结果
head(result)

print("write bed file starting...")

# 初始化一个新的结果数据框
df <- result
# 初始化一个新的结果数据框
result <- data.frame()

# 遍历数据框中的行
i <- 1
while (i <= nrow(df)) {
  if (i < nrow(df)) {
    # 检查相邻行是否满足条件：1. end 和 start 差值在 50 以内，2. SF 值不同
    if (abs(df$start[i + 1] - df$end[i]) <= 50 && df$SF[i] != df$SF[i + 1]) {
      # 如果满足条件，合并这两行
      new_row <- data.frame(
        Chr = df$Chr[i],                      # 保留染色体
        start = df$start[i],                  # 保留第一行的起始位置
        end = df$end[i + 1],                  # 保留第二行的结束位置
        strand = df$strand[i],                # 保留链方向
        dimer = paste0(df$SF[i], df$SF[i + 1]),  # 形成双聚体，例如 "JK"
        unit_count = df$unit_count[i] + df$unit_count[i + 1],  # 合并单元计数
        SF_count = paste0("(", df$unit_count[i], ",", df$unit_count[i + 1], ")")  # 记录两个 SF 的计数
      )
      # 将新行添加到结果中
      result <- rbind(result, new_row)
      
      # 跳过下一行，因为它已经被合并
      i <- i + 2
    } else {
      # 如果不满足条件，直接保留当前行
      new_row <- data.frame(
        Chr = df$Chr[i],
        start = df$start[i],
        end = df$end[i],
        strand = df$strand[i],
        dimer = df$SF[i],                    # 保留原始的 SF 值
        unit_count = df$unit_count[i],
        SF_count = paste0("(", df$unit_count[i], ")")  # 单个 SF 计数
      )
      result <- rbind(result, new_row)
      
      # 继续处理下一行
      i <- i + 1
    }
  } else {
    # 处理最后一行（如果是单独一行，无法与下一行合并）
    new_row <- data.frame(
      Chr = df$Chr[i],
      start = df$start[i],
      end = df$end[i],
      strand = df$strand[i],
      dimer = df$SF[i],                      # 保留原始的 SF 值
      unit_count = df$unit_count[i],
      SF_count = paste0("(", df$unit_count[i], ")")  # 单个 SF 计数
    )
    result <- rbind(result, new_row)
    
    i <- i + 1
  }
}

# 新增一列，记录 dimer 的倒序
result$reversed_dimer <- sapply(result$dimer, function(x) {
  # 分割字符，反转后再组合
  paste(rev(strsplit(x, NULL)[[1]]), collapse = "")
})

# 初始化一个新的结果数据框
HOR_sf <- result
result <- data.frame()

# 初始化索引
i <- 1
while (i <= nrow(HOR_sf)) {
  # 初始化合并的行数
  repeats_unit <- 1
  
  # 检查相邻行是否满足条件
  while (i < nrow(HOR_sf) &&
         abs(HOR_sf$start[i + 1] - HOR_sf$end[i]) < 5 &&         # 判断start和end差值
         HOR_sf$dimer[i] == HOR_sf$dimer[i + 1] &&               # 判断dimer相同
         HOR_sf$SF_count[i] == HOR_sf$SF_count[i + 1]) {         # 判断SF_count相同
    # 如果满足条件，合并行，更新结束位置并增加合并的行数
    HOR_sf$end[i] <- HOR_sf$end[i + 1]
    repeats_unit <- repeats_unit + 1
    i <- i + 1  # 跳到下一行继续合并
  }
  
  # 将当前行保存到结果数据框，并记录repeats_unit
  new_row <- data.frame(
    Chr = HOR_sf$Chr[i],
    start = HOR_sf$start[i - repeats_unit + 1],  # 合并后的起始位置
    end = HOR_sf$end[i],                        # 合并后的结束位置
    strand = HOR_sf$strand[i],
    dimer = HOR_sf$dimer[i],
    unit_count = HOR_sf$unit_count[i],
    SF_count = HOR_sf$SF_count[i],
    reversed_dimer = HOR_sf$reversed_dimer[i],
    repeats_unit = repeats_unit                 # 新增列记录合并的行数
  )
  
  # 添加新行到结果中
  result <- rbind(result, new_row)
  
  # 继续下一行
  i <- i + 1
}


# 查看结果
#print(result)



dimer_sf <- HOR_sf %>% select(Chr,start,end,strand,reversed_dimer,repeats_unit,SF_count,dimer) %>%
 filter( !nchar(dimer)==1)

write.table(dimer_sf, file = output_bed_file, row.names = FALSE,sep = "\t",col.names = FALSE,quote = FALSE)

print("write bed file done.")
print("write fa file starting...")

SF_sequences <- dimer_sf %>%
  group_by(Chr) %>%
  summarise(SF_sequence = paste(dimer, collapse = ""))

# 2. 按照 fasta 格式输出 
# 创建 fasta 文件
output_file <- output_fasta_file
fasta_lines <- c()

for (i in 1:nrow(SF_sequences)) {
  # 创建 fasta 的 header
  header <- paste0(">", SF_sequences$Chr[i])
  
  sequence <- SF_sequences$SF_sequence[i]
  
  # 将 header 和 sequence 加入到输出列表
  fasta_lines <- c(fasta_lines, header, sequence)
}

# 写入 fasta 文件
writeLines(fasta_lines, output_file)

