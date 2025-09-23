#include <bits/stdc++.h>
#include <filesystem>
using namespace std;
namespace fs = std::filesystem;

#include "utils.cpp"

int main(int argc, char** argv) {
    if (argc < 4) {
        cerr << "Usage: " << argv[0] << " <carpeta_counts> <carpeta_salida> <phi>\n";
        return 1;
    }

    string dir_path = argv[1];
    string out_dir = argv[2];
    double phi = stod(argv[3]);

    if (!fs::exists(out_dir)) {
        fs::create_directories(out_dir);
    }

    // Recorrer carpeta
    for (auto &entry : fs::directory_iterator(dir_path)) {
        if (!entry.is_regular_file()) continue;
        if (entry.path().extension() != ".csv") continue;

        string path = entry.path().string();
        uint64_t N;
        auto counts = read_counts(path, N);

        uint64_t U = static_cast<uint64_t>(ceil(phi * N));
        cout << "Procesando " << path << " con phi=" << phi << " (U=" << U << ")\n";

        // Archivo de salida: mismo nombre pero con sufijo _HH.csv
        string fname = entry.path().stem().string(); // sin extensión
        string out_file = out_dir + "/" + fname + "_HH.csv";

        ofstream out(out_file);
        if (!out.is_open()) {
            cerr << "No se pudo abrir archivo de salida: " << out_file << "\n";
            continue;
        }

        out << "kmer,count\n";
        size_t hh_count = 0;
        for (auto &p : counts) {
            if (p.second >= U) {
                out << p.first << "," << p.second << "\n";
                hh_count++;
            }
        }
        out.close();

        cout << "  → Guardados " << hh_count << " HH en " << out_file << "\n";
    }

    return 0;
}
