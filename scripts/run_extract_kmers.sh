#!/bin/bash
set -e  # detener si ocurre un error
cd "$(dirname "$0")/.." # ir a la raíz del proyecto

# --- Compilar extract_kmers.cpp ---
g++ -std=c++17 -O2 src/extract_kmers.cpp -o src/extract_kmers

# Carpeta de entrada y salida
DATA_DIR="data"
OUT_DIR="results/kmers"
LOG_FILE="results/memoria_extract_kmers.csv"

# Encabezado del log
echo "archivo,k,memoria_KB,tiempo_segundos" > "$LOG_FILE"


# Valores de k a usar
K_LIST=(21 31)

for K in "${K_LIST[@]}"; do
  echo ">>> Procesando k=$K"

  # Carpeta específica para cada k
  K_OUT_DIR="$OUT_DIR/k${K}"
  mkdir -p "$K_OUT_DIR"

  # Procesar todos los .fna en la carpeta /data
  for f in "$DATA_DIR"/*.fna; do
    BASENAME=$(basename "$f" .fna)

    COUNTS_FILE="$K_OUT_DIR/${BASENAME}_counts.csv"
    KMERS_FILE="$K_OUT_DIR/${BASENAME}_kmers.txt"

    # Ejecutar y capturar memoria y tiempo
    /usr/bin/time -v ./src/extract_kmers "$f" $K "$COUNTS_FILE" "$KMERS_FILE" \
      2> tmp.log

    # Extraer valores relevantes
    MEM=$(grep "Maximum resident set size" tmp.log | awk '{print $6}')
    TIME=$(grep "Elapsed (wall clock) time" tmp.log | awk '{print $8}')

    # Guardar en CSV
    echo "${BASENAME},${K},${MEM},${TIME}" >> "$LOG_FILE"
  done
done
