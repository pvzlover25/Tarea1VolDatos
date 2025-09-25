#!/bin/bash
set -e  # detener si ocurre un error
cd "$(dirname "$0")/.." # ir a la raíz del proyecto

# --- Compilar extract_hh.cpp ---
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

# Ejecutar extract_hh para k21
echo "Procesando k21..."
./src/extract_hh "$INPUT_K21" "$OUTPUT_K21" "$PHI"

# Ejecutar extract_hh para k31
echo "Procesando k31..."
./src/extract_hh "$INPUT_K31" "$OUTPUT_K31" "$PHI"

echo "Heavy Hitters extraídos exitosamente."