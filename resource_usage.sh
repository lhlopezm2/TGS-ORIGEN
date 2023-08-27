#!/usr/bin/bash

# START - Function to measure Time, RAM, Disk, and CPU
measure_resource_usage() {
    local path_out="$1"  # El primer argumento es el directorio de salida
    local sleep_duration="$2"  # El segundo argumento es la duración del intervalo de medición
    local logfile="$3"  # Archivo de registro

    max_memory_used_gb=0
    max_memory_used_mb=0
    min_memory_used_mb=$(free -m | awk 'NR==2 {print $3}')
    min_memory_used_gb=$(awk "BEGIN { printf \"%.4f\", $min_memory_used_mb * 1024 }") 

    max_disk_usage_GB=0
    max_disk_usage_MB=0
    min_disk_usage_MB=$(du -sm "${path_out}" | cut -f1)
    min_disk_usage_GB=$(awk "BEGIN {printf \"%.5f\", $min_disk_usage_MB / 1024}")

    cpu_max=0

    start_time=$(date +%s)

    while : # This will run forever until we explicitly kill it
    do
        # Time
        duration_seconds=$(($(date +%s) - start_time))
        duration_minutes=$(awk "BEGIN {printf \"%.5f\", $duration_seconds / 60}")

        # RAM Memory
        memory_used_mb=$(free -m | awk 'NR==2 {print $3}')
        memory_used_gb=$(awk "BEGIN { printf \"%.4f\", $memory_used_mb / 1024 }") 
        if (( $(awk "BEGIN { print ($memory_used_gb >= $max_memory_used_gb) ? 1 : 0 }") )); then
            max_memory_used_gb=$memory_used_gb
            max_memory_used_mb=$memory_used_mb # Updated to use the value in MB
        fi
        if (( $(awk "BEGIN { print ($memory_used_gb <= $min_memory_used_gb) ? 1 : 0 }") )); then
            min_memory_used_gb=$memory_used_gb
            min_memory_used_mb=$memory_used_mb # Updated to use the value in MB
        fi
        MinMaxMB=$((max_memory_used_mb - min_memory_used_mb))
        MinMaxGB=$(awk "BEGIN { printf \"%.4f\", $MinMaxMB / 1024 }")

        # Disk Usage
        disk_usage_MB=$(du -sm "${path_out}" | cut -f1)
        disk_usage_GB=$(awk "BEGIN {printf \"%.5f\", $disk_usage_MB / 1024}")
        if (( $(awk "BEGIN { print ($disk_usage_GB >= $max_disk_usage_GB) ? 1 : 0 }") )); then
            max_disk_usage_GB=$disk_usage_GB
            max_disk_usage_MB=$disk_usage_MB            
        fi
        if (( $(awk "BEGIN { print ($disk_usage_GB <= $min_disk_usage_GB) ? 1 : 0 }") )); then
            min_disk_usage_GB=$disk_usage_GB
            min_disk_usage_MB=$disk_usage_MB                        
        fi
        
        # CPU Usage
        total_cpus=$(nproc)
        cpu_usage=()
        cpu_mean=0
        for ((cpu=1; cpu<=total_cpus; cpu++)); do
            cpu_info=$(top -bn1 -p $cpu | grep "Cpu(s)" | awk '{print $2 + $4}')
            cpu_usage[$cpu]=${cpu_info%.*}  # Remove decimal part if any
            cpu_mean=$((cpu_mean + cpu_usage[$cpu]))
        done
        cpu_mean=$((cpu_mean / total_cpus))
        if (( $(awk "BEGIN { print ($cpu_mean >= $cpu_max) ? 1 : 0 }") )); then
            cpu_max=$cpu_mean
        fi
        
        # Redirigir la salida al archivo de registro especificado por 'logfile'
        {
            echo -e "\n\n Resources used: "
            echo "Time: ${duration_seconds} seconds (${duration_minutes} minutes)"
            echo "RAM Memory Actual: ${memory_used_mb} MB (${memory_used_gb} GB)"
            echo "RAM Memory Min: ${min_memory_used_mb} MB (${min_memory_used_gb} GB)"
            echo "RAM Memory Max: ${max_memory_used_mb} MB (${max_memory_used_gb} GB)"
            echo "RAM Memory MinMax: ${MinMaxMB} MB (${MinMaxGB} GB)"
            echo "Disk Actual: ${disk_usage_MB} MB (${disk_usage_GB} GB)"
            echo "Disk Min: ${min_disk_usage_MB} MB (${min_disk_usage_GB} GB)"
            echo "Disk Max: ${max_disk_usage_MB} MB (${max_disk_usage_GB} GB)"
            echo "CPU Mean Max: ${cpu_max}%"
            echo "CPU Mean: ${cpu_mean}%"
            echo -e "CPUs: "
            for ((cpu=1; cpu<=total_cpus; cpu++)); do
                # Mostrar CPU usage
                echo "CPU $cpu Usage: ${cpu_usage[$cpu]}%"
            done
        } >> "$logfile"  # Usar >> para agregar al archivo existente

        sleep "$sleep_duration" # Ajusta este valor para cambiar el intervalo de medición
    done
}
# END - Function to measure Time, RAM, Disk, and CPU
