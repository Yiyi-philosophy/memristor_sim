#!/bin/bash  
  
# 设置日志文件路径  
output_file="memory_usage.log"  
  
# 写入日志文件的标题  
echo "Date Time, Memory Used (GB)" > "$output_file"  
  
# 循环，每隔1秒记录一次内存使用情况  
while true; do  
    # 使用top命令获取内存使用情况，并通过grep和awk提取used内存部分  
    # 注意：我们需要将MiB转换为GiB，因此要除以1024  
    used_mem_mib=$(top -b -n 1 | grep "MiB Mem" | awk '{print $7}')  
    used_mem_gib=$(echo "scale=2; $used_mem_mib / 1024" | bc)  
      
    # 获取当前日期和时间  
    datetime=$(date '+%Y-%m-%d %H:%M:%S')  
      
    # 将日期时间和内存使用情况写入日志文件  
    echo "$datetime, $used_mem_gib" >> "$output_file"  
      
    # 等待1秒  
    sleep 3
done  
