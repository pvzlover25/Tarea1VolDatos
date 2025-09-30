#!/bin/bash
set -e #detener si ocurre un error
cd "$(dirname "$0")/.." # ir a la raíz del proyecto

#compilar
g++ -std=c++17 -O2 /home/pvzlover25/Desktop/Tarea1VolDatos/src/sketches.cpp -o sketches

#Carpeta de entrada y salida
DATA_DIR_21="results/kmers/k21"
DATA_DIR_31="results/kmers/k31"
OUT_DIR_21="results/sketches/k21"
OUT_DIR_31="results/sketches/k31"

#crear carpetas de salida si no existen
mkdir -p "$OUT_DIR_21"
mkdir -p "$OUT_DIR_31"

#procesar todos los .fna en DATA_DIR_21
for f in "$DATA_DIR_21"/*.csv; do
  BASENAME=$(basename "$f" .csv)
  CSV_FILE="$OUT_DIR_21/${BASENAME}_sketch_result.csv"
  echo "Archivo: $f"
  ./sketches "$f" "$CSV_FILE"
done

#procesar todos los .fna en DATA_DIR_31
for f in "$DATA_DIR_31"/*.csv; do
  BASENAME=$(basename "$f" .csv)
  CSV_FILE="$OUT_DIR_31/${BASENAME}_sketch_result.csv"
  echo "Archivo: $f"
  ./sketches "$f" "$CSV_FILE"
done
