#include <bits/stdc++.h>
#include <filesystem>
using namespace std;
namespace fs = std::filesystem;

#include "utils.cpp"

int main(int argc, char** argv) {
    if (argc < 4) {
        cerr << "Usage: " << argv[0] << " <carpeta_counts> <salida.csv> <phi1,phi2,...>\n";
        return 1;
    }

    string dir_path = argv[1];
    string out_file = argv[2];
    vector<string> phi_strs = split(argv[3], ',');
    vector<double> phis;
    for (auto &s : phi_strs) {
        phis.push_back(stod(s));
    }

    ofstream out(out_file);
    if (!out.is_open()) {
        cerr << "No se pudo abrir archivo de salida\n";
        return 1;
    }

    // Escribir header
    out << "archivo,N";
    for (auto &phi : phis) {
        out << ",phi=" << phi;
    }
    out << "\n";

    // Recorrer carpeta
    for (auto &entry : fs::directory_iterator(dir_path)) {
        if (!entry.is_regular_file()) continue;
        string path = entry.path().string();

        // solo procesar archivos .csv
        if (entry.path().extension() != ".csv") continue;

        uint64_t N;
        auto counts = read_counts(path, N);

        out << entry.path().filename().string(); // nombre de archivo, no ruta completa
        out << "," << N; // total de k-mers
        for (auto phi : phis) {
            uint64_t U = static_cast<uint64_t>(std::ceil(phi * N));
            cout << "Procesando " << path << " con phi=" << phi << " (U=" << U << ")\n";
            size_t hh_count = count_if(counts.begin(), counts.end(),
                           [&](const pair<string,uint64_t>& p){ return p.second >= U; });

            out << "," << hh_count;
        }
        out << "\n";
    }

    out.close();
    return 0;
}
