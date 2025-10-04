#!/bin/bash
set -e  # detener si ocurre un error
cd "$(dirname "$0")/.." # ir a la raíz del proyecto

# --- Compilar extract_hh.cpp ---
echo ">>> Compilando extract_hh.cpp..."
g++ -std=c++17 -O2 src/extract_hh.cpp -o src/extract_hh

PHI="1e-6"

# Carpetas de entrada y salida
INPUT_K21="results/kmers/k21"
OUTPUT_K21="results/hh/k21"

INPUT_K31="results/kmers/k31"
OUTPUT_K31="results/hh/k31"

# Crear carpetas de salida si no existen
mkdir -p "$OUTPUT_K21"
mkdir -p "$OUTPUT_K31"

# Archivo de log único
LOG_FILE="results/memoria_extract_hh.csv"
echo "k,phi,memoria_KB,tiempo_segundos" > "$LOG_FILE"

# Ejecutar extract_hh para k=21 con medición
echo "Procesando k21..."
/usr/bin/time -v ./src/extract_hh "$INPUT_K21" "$OUTPUT_K21" "$PHI" 2> tmp.log
MEM=$(grep "Maximum resident set size" tmp.log | awk '{print $6}')
TIME=$(grep "Elapsed (wall clock) time" tmp.log | awk '{print $8}')
echo "21,${PHI},${MEM},${TIME}" >> "$LOG_FILE"

# Ejecutar extract_hh para k=31 con medición
echo "Procesando k31..."
/usr/bin/time -v ./src/extract_hh "$INPUT_K31" "$OUTPUT_K31" "$PHI" 2> tmp.log
MEM=$(grep "Maximum resident set size" tmp.log | awk '{print $6}')
TIME=$(grep "Elapsed (wall clock) time" tmp.log | awk '{print $8}')
echo "31,${PHI},${MEM},${TIME}" >> "$LOG_FILE"

rm -f tmp.log

echo "Heavy Hitters extraídos exitosamente. Log guardado en $LOG_FILE"
