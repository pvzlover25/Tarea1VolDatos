// Otro intento de implementar los sketches por si hay un error con el primer codigo
// los errores dan todos 0 por ahora igual que el anterior e itera menos veces
// No se cual será mejor

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <cstdint>
#include <limits>
#include <functional>
#include <random>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include "utils.cpp"

using namespace std;

class CountSketch {
private:
    int d;  // numero de filas
    int w;  // ancho de fila
    vector<vector<int64_t>> C;  // matriz de contadores
    vector<hash<string>> hash_funcs;
    mt19937 rng;
    
    // funcion hash h_j: mapea k-mer a [0, w-1]
    inline int h(int j, const string& kmer) {
        size_t hash_val = hash_funcs[0](kmer);
        hash_val ^= (hash_val << 13) ^ j;
        hash_val ^= (hash_val >> 7);
        hash_val ^= (hash_val << 17);
        return hash_val % w;
    }
    
    // funcion signo g_j: mapea k-mer a {-1, +1}
    inline int g(int j, const string& kmer) {
        size_t hash_val = hash_funcs[0](kmer);
        hash_val ^= (hash_val >> 11) ^ (j * 0x9e3779b9);
        return ((hash_val & 1) == 0) ? 1 : -1;
    }
    
public:
    CountSketch(int rows, int width) : d(rows), w(width), rng(42) {
        C.resize(d, vector<int64_t>(w, 0));
        hash_funcs.resize(1); 
    }
    
    void insert(const string& kmer, int64_t c = 1) {
        string cano = canonical(kmer);
        for (int j = 0; j < d; j++) {
            int pos = h(j, cano);
            int sign = g(j, cano);
            C[j][pos] += sign * c;
        }
    }
    
    int64_t estimate(const string& kmer) {
        string cano = canonical(kmer);
        vector<int64_t> estimates;
        
        for (int j = 0; j < d; j++) {
            int pos = h(j, cano);
            int sign = g(j, cano);
            estimates.push_back(sign * C[j][pos]);
        }
        
        sort(estimates.begin(), estimates.end());
        return estimates[d / 2];
    }
    
    size_t memory_usage() const {
        return d * w * sizeof(int64_t);
    }
    
    void print_stats() const {
        cout << "CountSketch: d=" << d << ", w=" << w 
                  << ", memoria=" << (memory_usage() / 1024.0) << " KB\n";
    }
};


class TowerSketch {
private:
    struct Layer {
        int bits;           // bits por contador
        int w;              // ancho de la capa
        int num_rows;       // numero de filas x capa
        vector<vector<uint64_t>> counters;
        uint64_t max_val;   // valor max
        
        Layer(int b, int width, int rows) : bits(b), w(width), num_rows(rows) {
            counters.resize(num_rows, vector<uint64_t>(w, 0));
            max_val = (1ULL << bits) - 1;
        }
    };
    
    vector<Layer> layers;
    hash<string> hash_func;
    
    inline int h(int j, const string& kmer, int w) {
        size_t hash_val = hash_func(kmer);
        hash_val ^= (hash_val << 13) ^ j;
        hash_val ^= (hash_val >> 7);
        hash_val ^= (hash_val << 17);
        return hash_val % w;
    }
    
public:
    TowerSketch(const vector<tuple<int, int, int>>& config) {
        for (const auto& [bits, width, rows] : config) {
            layers.emplace_back(bits, width, rows);
        }
    }
    
    static TowerSketch create_default(int base_width = 10000) {
        vector<tuple<int, int, int>> config = {
            {64, base_width / 10, 2},      // 64 bits
            {32, base_width / 5, 2},       // 32 bits
            {16, base_width / 2, 2},       // 16 bits
            {8, base_width, 1},            // 8 bits
            {4, base_width * 2, 1}         // 4 bits
        };
        return TowerSketch(config);
    }
    
    // Insertar k-mer con incremento c
    void insert(const string& kmer, uint64_t c = 1) {
        string cano = canonical(kmer);
        
        for (auto& layer : layers) {
            for (int j = 0; j < layer.num_rows; j++) {
                int pos = h(j, cano, layer.w);
                
                if (layer.counters[j][pos] + c <= layer.max_val) {
                    layer.counters[j][pos] += c;
                } else {
                    layer.counters[j][pos] = layer.max_val;  // Saturar
                }
            }
        }
    }
    
    uint64_t estimate(const string& kmer) {
        string cano = canonical(kmer);
        uint64_t min_estimate = numeric_limits<uint64_t>::max();
        
        for (const auto& layer : layers) {
            for (int j = 0; j < layer.num_rows; j++) {
                int pos = h(j, cano, layer.w);
                uint64_t val = layer.counters[j][pos];
                
                if (val < layer.max_val) {
                    min_estimate = min(min_estimate, val);
                }
            }
        }
        
        return min_estimate;
    }
    
    size_t memory_usage() const {
        size_t total = 0;
        for (const auto& layer : layers) {
            total += layer.num_rows * layer.w * (layer.bits / 8);
        }
        return total;
    }
    
    void print_stats() const {
        cout << "TowerSketch:\n";
        for (size_t i = 0; i < layers.size(); i++) {
            cout << "  Layer " << i << ": " << layers[i].bits << " bits, "
                      << layers[i].w << " width, " << layers[i].num_rows << " rows\n";
        }
        cout << "  Memoria total: " << (memory_usage() / 1024.0) << " KB\n";
    }
};

unordered_map<string, uint64_t> read_kmer_counts(const string& filename) {
    unordered_map<string, uint64_t> counts;
    ifstream file(filename);
    
    if (!file.is_open()) {
        cerr << "Error: No se pudo abrir el archivo " << filename << "\n";
        return counts;
    }
    
    string line;
    getline(file, line);
    
    while (getline(file, line)) {
        istringstream ss(line);
        string kmer;
        uint64_t count;
        
        if (getline(ss, kmer, ',')) {
            ss >> count;
            counts[kmer] = count;
        }
    }
    
    file.close();
    return counts;
}

struct EvaluationMetrics {
    double avg_absolute_error;
    double avg_relative_error;
    double max_error;
    int num_evaluated;
};

template<typename SketchType>
EvaluationMetrics evaluate_sketch(SketchType& sketch, const unordered_map<string, uint64_t>& ground_truth) {
    EvaluationMetrics metrics = {0.0, 0.0, 0.0, 0};
    
    for (const auto& [kmer, true_freq] : ground_truth) {
        auto estimated = sketch.estimate(kmer);
        double abs_error = abs((int64_t)estimated - (int64_t)true_freq);
        double rel_error = (true_freq > 0) ? abs_error / true_freq : 0.0;
        
        metrics.avg_absolute_error += abs_error;
        metrics.avg_relative_error += rel_error;
        metrics.max_error = max(metrics.max_error, abs_error);
        metrics.num_evaluated++;
    }
    
    if (metrics.num_evaluated > 0) {
        metrics.avg_absolute_error /= metrics.num_evaluated;
        metrics.avg_relative_error /= metrics.num_evaluated;
    }
    
    return metrics;
}


int main(int argc, char** argv) {
    if (argc < 3) {
        cerr << "Usage: " << argv[0] << " <kmer_counts.csv> <output_prefix>\n";
        cerr << "  Lee k-mers y sus frecuencias desde CSV\n";
        cerr << "  Implementa CountSketch y TowerSketch para evaluacion\n";
        return 1;
    }
    
    string counts_file = argv[1];
    string output_prefix = argv[2];
    
    cerr << "Leyendo k-mers desde CSV...\n";
    auto freq = read_kmer_counts(counts_file);
    
    if (freq.empty()) {
        cerr << "Error: No se pudieron leer k-mers del archivo\n";
        return 1;
    }
    
    uint64_t total_kmers = 0;
    for (const auto& [kmer, count] : freq) {
        total_kmers += count;
    }
    
    cerr << "Ground truth: " << freq.size() << " k-mers unicos, "
              << total_kmers << " k-mers totales\n\n";
    
    cerr << "CountSketch\n";
    vector<pair<int, int>> cs_configs = {
        {3, 5000}, {5, 3000}, {7, 2000}, {5, 5000}
    };
    
    for (auto [d, w] : cs_configs) {
        CountSketch cs(d, w);
        
        for (const auto& [kmer, count] : freq) {
            cs.insert(kmer, count);
        }
        
        auto metrics = evaluate_sketch(cs, freq);
        
        cs.print_stats();
        cout << "  Error absoluto promedio: " << metrics.avg_absolute_error << "\n";
        cout << "  Error relativo promedio: " << (metrics.avg_relative_error * 100) << "%\n";
        cout << "  Error maximo: " << metrics.max_error << "\n\n";
    }
    
    cerr << "TowerSketch\n";
    TowerSketch ts1 = TowerSketch::create_default(5000);
    for (const auto& [kmer, count] : freq) {
        ts1.insert(kmer, count);
    }
    auto metrics1 = evaluate_sketch(ts1, freq);
    ts1.print_stats();
    cout << "  Error absoluto promedio: " << metrics1.avg_absolute_error << "\n";
    cout << "  Error relativo promedio: " << (metrics1.avg_relative_error * 100) << "%\n";
    cout << "  Error maximo: " << metrics1.max_error << "\n\n";
    
    vector<tuple<int, int, int>> config2 = {
        {64, 500, 2}, {32, 1000, 2}, {16, 2000, 3}, {8, 5000, 2}, {4, 10000, 1}
    };
    TowerSketch ts2(config2);
    for (const auto& [kmer, count] : freq) {
        ts2.insert(kmer, count);
    }
    auto metrics2 = evaluate_sketch(ts2, freq);
    ts2.print_stats();
    cout << "  Error absoluto promedio: " << metrics2.avg_absolute_error << "\n";
    cout << "  Error relativo promedio: " << (metrics2.avg_relative_error * 100) << "%\n";
    cout << "  Error maximo: " << metrics2.max_error << "\n\n";
    
    cerr << "Calibracion completada\n";

    ofstream results(output_prefix + ".csv");
    if (!results.is_open()) {
        cerr << "no se pudo abrir el archivo de salida\n";
        return 1;
    }

    results << "sketch,params,avg_abs_error,avg_rel_error,max_error\n";

    { 
        int idx = 0;
        for (auto [d, w] : cs_configs) {
            CountSketch cs(d, w);
            for (const auto& [kmer, count] : freq) {
                cs.insert(kmer, count);
            }
            auto metrics = evaluate_sketch(cs, freq);

            results << "CountSketch" << ",d=" << d << ";w=" << w << "," << metrics.avg_absolute_error
            << "," << metrics.avg_relative_error << "," << metrics.max_error << "\n";
            idx++;
        }
    }

    {
        TowerSketch ts1 = TowerSketch::create_default(5000);
        for (const auto& [kmer, count] : freq) ts1.insert(kmer, count);
        auto metrics1 = evaluate_sketch(ts1, freq);
        results << "TowerSketch,default,"
                << metrics1.avg_absolute_error
                << "," << metrics1.avg_relative_error
                << "," << metrics1.max_error << "\n";

        vector<tuple<int, int, int>> config2 = {
            {64, 500, 2}, {32, 1000, 2}, {16, 2000, 3}, {8, 5000, 2}, {4, 10000, 1}
        };
        TowerSketch ts2(config2);
        for (const auto& [kmer, count] : freq) ts2.insert(kmer, count);
        auto metrics2 = evaluate_sketch(ts2, freq);
        results << "TowerSketch,custom,"
                << metrics2.avg_absolute_error
                << "," << metrics2.avg_relative_error
                << "," << metrics2.max_error << "\n";
    }

    results.close();
    cerr << "resultados en " << output_prefix << ".csv\n";

        
    return 0;
}

// compila con: g++ -std=c++17 -O3 sketchtest.cpp -o sketches
// sketches kmer_counts.csv output_prefix

// extract_kmers GCA_023315275.2_PDT001299634.2_genomic 3 counts.csv kmers.csv