#!/bin/bash
set -e  # detener si ocurre un error
cd "$(dirname "$0")/.." # ir a la raíz del proyecto

# --- Compilar extract_kmers.cpp ---
g++ -std=c++17 -O2 src/extract_kmers.cpp -o src/extract_kmers

# Carpeta de entrada y salida
DATA_DIR="data"
OUT_DIR="results/kmers"

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

    echo "Archivo: $f"
    # Limpiar salidas anteriores
    rm -f "$COUNTS_FILE" "$KMERS_FILE"

    ./src/extract_kmers "$f" $K "$COUNTS_FILE" "$KMERS_FILE"
  done
done
