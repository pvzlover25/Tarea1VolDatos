# Tarea1VolDatos

---
## Actividad 1
### Requisitos

- Compilador **C++17** o superior.  
- Estructura recomendada:
  ```
  data/                 # Archivos .fna
  results/kmers/k21/   # Counts.csv para k=21
  results/kmers/k31/   # Counts.csv para k=31
  results/hh/k21/       # HH exactos (phi elegido)
  results/hh/k31/       # HH exactos (phi elegido)
  results/phi_selection # Tablas para decidir phi
  scripts/              # Scripts auxiliares (.sh)
  ```
---

### Programas incluidos

#### 1. `extract_kmers`

- **Descripción:**  
  Lee un archivo `.fna` y extrae todos los k-mers **canónicos** de tamaño `k`.  
  Genera:  
  - `*_counts.csv` → cada k-mer único y su frecuencia.  
  - `*_kmers.txt` → lista de cada kmer canonico encontrado (con repeticiones).  

- **Compilación:**

- Para un archivo unico:
  ```bash
  g++ -std=c++17 -O2 extract_kmers.cpp -o extract_kmers
  ```

    - **Uso:**
  ```bash
  ./extract_kmers <input.fna> <k> <output_prefix>
  ```

  Ejemplo:
  ```bash
  ./extract_kmers ../data/genoma1.fna 21 ../results/kmers/k21/genoma1
  ```

---

#### 2. `calc_phi_table`

- **Descripción:**  
  Recorre una carpeta con archivos `_counts.csv` y calcula el número de HH encontrados para diferentes valores de ϕ.

- **Compilación:**
  ```bash
  g++ -std=c++17 -O2 calc_phi_table.cpp -o calc_phi_table
  ```

- **Uso:**
  ```bash
  ./calc_phi_table <carpeta_counts> <salida.csv> <phi1,phi2,...>
  ```

  Ejemplo:
  ```bash
  ./calc_phi_table ../results/kmers/k21 ../results/phi_selection/k21_phi_table.csv "1e-4,1e-5,1e-6,1e-7"
  ```

---

#### 3. `extract_hh`

- **Descripción:**  
  Filtra los Heavy Hitters de cada archivo `_counts.csv` dado un ϕ.  
  Genera un archivo `*_HH.csv` por genoma con los k-mers que cumplen la condición.  

- **Compilación:**
  ```bash
  g++ -std=c++17 -O2 extract_hh.cpp -o extract_hh
  ```

- **Uso:**
  ```bash
  ./extract_hh <carpeta_counts> <carpeta_salida> <phi>
  ```

  Ejemplo:
  ```bash
  ./extract_hh ../results/kmers/k21 ../results/hh/k21 1e-6
  ```

---

### Flujo de trabajo sugerido

1. **Extraer k-mers y frecuencias**  
   ```bash
   bash run_extract_kmers.sh 
   ```

2. **Calcular tabla de HH vs φ en genomas representativos**  
   ```bash
   ./calc_phi_table ../results/kmers/k21 ../results/phi_selection/k21_phi_table.csv "1e-4,1e-5,1e-6,1e-7"
   ./calc_phi_table ../results/kmers/k31 ../results/phi_selection/k31_phi_table.csv "1e-4,1e-5,1e-6,1e-7"
   ```

3. **Elegir φ global**   

4. **Ajustar φ en `run_extract_hh.sh`**

5. **Extraer HH exactos para todos los genomas**  
   ```bash
   bash run_extract_hh.sh 
   ```

---