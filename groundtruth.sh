#!/bin/bash
set -e  # detener si ocurre un error

# Carpeta de entrada y salida
DATA_DIR="data"
OUT_DIR="results/ground_truth"

# Valores de k a usar
K_LIST=(21 31)

# Valores de phi (puedes ajustar según memoria)
PHI_LIST=(1e-3 1e-4 1e-5)

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
    HH_FILE="$K_OUT_DIR/${BASENAME}_HH.csv"

    echo "    Archivo: $f"
    # Limpiar salidas anteriores
    rm -f "$COUNTS_FILE" "$KMERS_FILE" "$HH_FILE"

    ./src/extract_kmers "$f" $K "$COUNTS_FILE" "$KMERS_FILE"
  done
done
