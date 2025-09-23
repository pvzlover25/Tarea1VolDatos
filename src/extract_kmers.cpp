#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <unordered_map>
#include <algorithm>
#include <cstdint>

#include "utils.cpp"

// Helper function to check if a character is a valid DNA base
inline bool is_valid_base(char c) {
    return c == 'A' || c == 'C' || c == 'G' || c == 'T';
}

int main(int argc, char** argv) {
    if (argc < 5) {
        std::cerr << "Usage: " << argv[0] << " <input.fasta> <k> <output_counts.csv> <output_kmers.csv>\n";
        return 1;
    }

    std::string input = argv[1];
    int k = std::stoi(argv[2]);
    std::string output_counts = argv[3];
    std::string output_kmers = argv[4];
    std::ofstream out_kmers(output_kmers);

    std::cerr << "Extrayendo k-mers de largo " << k << " desde " << input << "\n";

    auto seqs = read_fasta(input);

    std::unordered_map<std::string, uint64_t> freq;
    uint64_t total_kmers = 0;

    for (auto &seq : seqs) {
        // limpiar Ns
        std::string clean;
        for (char c : seq) {
            if (is_valid_base(c)) clean.push_back(c);
        }

        for (size_t i = 0; i + k <= clean.size(); i++) {
            std::string sub = clean.substr(i, k);
            // if (sub.find('N') != std::string::npos) continue; // saltar si hay Ns
            std::string cano = canonical(sub);
            freq[cano]++;
            total_kmers++;

            out_kmers << cano << "\n";
        }
    }

    std::ofstream out_counts(output_counts);
    out_counts << "kmer,count\n";
    for (auto &kv : freq) {
        out_counts << kv.first << "," << kv.second << "\n";
    }
    out_counts.close();
    out_kmers.close();

    std::cerr << "------Extracción de K-mers completada.------\n";

    std::cerr << "Total de k-mers: " << total_kmers << "\n";
    std::cerr << "K-mers únicos: " << freq.size() << "\n";
    
    std::cerr << "------------------------------------------------------------\n";
    return 0;
}