#ifndef UTILS
#define UTILS

#include <string>
#include <vector>
#include <fstream>
#include <algorithm>
#include <filesystem>
using namespace std;

// Función de reverso complementario
string revcomp(const string &kmer) {
    string rc(kmer.size(), 'N');
    for (int i = 0; i < (int)kmer.size(); i++) {
        char c = kmer[kmer.size() - 1 - i];
        switch(c) {
            case 'A': rc[i] = 'T'; break;
            case 'C': rc[i] = 'G'; break;
            case 'G': rc[i] = 'C'; break;
            case 'T': rc[i] = 'A'; break;
            default: rc[i] = 'N'; break;
        }
    }
    return rc;
}

// Función para obtener k-mer canónico
string canonical(const string &kmer) {
    string rc = revcomp(kmer);
    return min(kmer, rc); // lexicográficamente menor
}

// Lectura FASTA
vector<string> read_fasta(const string &path) {
    ifstream in(path);

    if (!in.is_open()) {
        throw runtime_error("Error: Unable to open file " + path);
    }

    vector<string> seqs;
    string line, seq;

    while (getline(in, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') {
            if (!seq.empty()) {
                seqs.push_back(seq);
                seq.clear();
            }
        } else {
            for (char &c : line) c = toupper(c);
            seq += line;
        }
    }
    if (!seq.empty()) seqs.push_back(seq);
    
    in.close();
    return seqs;
}

// Función para dividir string por delimitador
vector<string> split(const string &s, char delim) {
    vector<string> elems;
    string item;
    stringstream ss(s);
    while (getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

// Lee conteos desde archivo CSV (ignora header), devuelve vector de pares <kmer, counts>
vector<pair<string,uint64_t>> read_counts(const string &path, uint64_t &N) {
    ifstream in(path);
    if (!in.is_open()) {
        cerr << "Error abriendo archivo: " << path << "\n";
        exit(1);
    }
    string line;
    getline(in, line); // saltar header
    vector<pair<string,uint64_t>> counts;
    N = 0;
    while (getline(in, line)) {
        if (line.empty()) continue;
        auto parts = split(line, ',');
        if (parts.size() != 2) continue;
        string kmer = parts[0];
        uint64_t c = stoull(parts[1]);
        counts.emplace_back(kmer, c);
        N += c;
    }
    return counts;
}

#endif
