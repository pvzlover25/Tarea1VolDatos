// sketches.cpp
// Compilar: g++ -O3 -std=c++17 sketches.cpp -o sketches
// Uso: ./sketches counts.csv kmers.csv [phi]
// Requiere utils.cpp en el mismo directorio (con canonical(), read_counts(), read_fasta()).

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>
#include <stdexcept>
#include <cstdint>
#include <cmath>
#include <utility>
#include <limits>
#include "utils.cpp" // archivo con canonical(), read_counts(), read_fasta()

using namespace std;

static inline uint64_t hash_u64(const string &s, uint64_t salt = 1469598103934665603ULL) {

    uint64_t h = std::hash<std::string>{}(s);
    h ^= salt + 0x9e3779b97f4a7c15ULL + (h<<6) + (h>>2);
    return h;
}

static inline uint64_t hash_u64_from_u64(uint64_t x, uint64_t salt) {
    x ^= salt + 0x9e3779b97f4a7c15ULL + (x<<6) + (x>>2);
    uint64_t z = (x + 0x9e3779b97f4a7c15ULL);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    z = z ^ (z >> 31);
    return z;
}

// CountSketch
// d fila  x  w colum, counters int32_t
class CountSketch {
public:
    CountSketch(int d = 5, int w = 1<<14, uint64_t seed = 12345ULL)
        : d_(d), w_(w), seed_(seed), table_((size_t)d*(size_t)w, 0) {
        if (d_ <= 0 || w_ <= 0) throw runtime_error("Invalid CS dims");
    }

    void insert(const string &kmer) {
        string cano = canonical(kmer);
        uint64_t key = hash_u64(cano, seed_);
        for (int i = 0; i < d_; ++i) {
            uint64_t h = hash_u64_from_u64(key, (uint64_t)i ^ seed_);
            uint32_t idx = (uint32_t)(h % (uint64_t)w_);
            int s = ((h >> 1) & 1) ? 1 : -1; // pseudo sign from hash bit
            size_t pos = (size_t)i * (size_t)w_ + idx;
            int32_t val = table_[pos];
            val += s;
            table_[pos] = val;
        }
    }

    int64_t estimate(const string &kmer) const {
        string cano = canonical(kmer);
        uint64_t key = hash_u64(cano, seed_);
        vector<int64_t> estimates; estimates.reserve(d_);
        for (int i = 0; i < d_; ++i) {
            uint64_t h = hash_u64_from_u64(key, (uint64_t)i ^ seed_);
            uint32_t idx = (uint32_t)(h % (uint64_t)w_);
            int s = ((h >> 1) & 1) ? 1 : -1;
            size_t pos = (size_t)i * (size_t)w_ + idx;
            int32_t v = table_[pos];
            estimates.push_back((int64_t)s * (int64_t)v);
        }

        size_t m = estimates.size()/2;
        nth_element(estimates.begin(), estimates.begin()+m, estimates.end());
        return estimates[m];
    }

    size_t memory_bytes() const {
        return table_.size() * sizeof(int32_t) + sizeof(*this);
    }

private:
    int d_;
    int w_;
    uint64_t seed_;
    vector<int32_t> table_;
};

//TowerSketch
// r filas, w columns (r >= 1). 
// Al insertar: calcula r posiciones, obtiene el mínimo m, incrementa solo las celdas iguales a m (Actualización Conservadora).
class TowerSketch {
public:
    TowerSketch(const vector<pair<int,int>> &levels, uint64_t seed = 54321ULL)
        : seed_(seed) {
        if (levels.empty()) throw runtime_error("TowerSketch debe tener 1 nivel");
        for (size_t i = 0; i < levels.size(); ++i) {
            int r = levels[i].first;
            int w = levels[i].second;
            if (r <= 0 || w <= 0) throw runtime_error("Invalid level params");
            Level L;
            L.r = r; L.w = w;
            L.seed = seed_ + (uint64_t)i * 0x9e3779b97f4a7c15ULL;
            L.table.assign((size_t)r * (size_t)w, 0u);
            levels_.push_back(std::move(L));
        }
    }

    void insert(const string &kmer) {
        string cano = canonical(kmer);
        uint64_t key = hash_u64(cano, seed_);
        for (const Level &lvl_const : levels_) {

            Level &lvl = const_cast<Level&>(lvl_const);
            vector<uint32_t> idxs; idxs.reserve(lvl.r);
            vector<uint32_t> vals; vals.reserve(lvl.r);
            for (int i = 0; i < lvl.r; ++i) {
                uint64_t h = hash_u64_from_u64(key, lvl.seed + (uint64_t)i);
                uint32_t idx = (uint32_t)(h % (uint64_t)lvl.w);
                idxs.push_back(idx);
                vals.push_back(lvl.table[(size_t)i * (size_t)lvl.w + idx]);
            }
            uint32_t m = *min_element(vals.begin(), vals.end());
            for (int i = 0; i < lvl.r; ++i) {
                if (vals[i] == m) {

                    lvl.table[(size_t)i * (size_t)lvl.w + idxs[i]]++;
                }
            }
        }
    }


    uint64_t estimate(const string &kmer) const {
        string cano = canonical(kmer);
        uint64_t key = hash_u64(cano, seed_);
        uint64_t best = UINT64_MAX;
        for (const Level &lvl : levels_) {
            uint64_t m = UINT64_MAX;
            for (int i = 0; i < lvl.r; ++i) {
                uint64_t h = hash_u64_from_u64(key, lvl.seed + (uint64_t)i);
                uint32_t idx = (uint32_t)(h % (uint64_t)lvl.w);
                uint64_t v = lvl.table[(size_t)i * (size_t)lvl.w + idx];
                if (v < m) m = v;
            }
            if (m < best) best = m;
        }
        return best == UINT64_MAX ? 0 : best;
    }

    size_t memory_bytes() const {
        size_t s = sizeof(*this);
        for (const Level &lvl : levels_) s += lvl.table.size() * sizeof(uint32_t);
        return s;
    }

private:
    struct Level {
        int r, w;
        uint64_t seed;
        vector<uint32_t> table;
    };
    uint64_t seed_;
    vector<Level> levels_;
};

vector<pair<string,uint64_t>> load_counts_and_N(const string &path, uint64_t &N) {
    return read_counts(path, N); 
}

struct EvalResult {
    string method;
    size_t mem_bytes;
    double mean_abs_err;
    double median_abs_err;
    double max_abs_err;
    uint64_t tp, fp, fn;
};

EvalResult evaluate_countsketch(const vector<pair<string,uint64_t>> &counts, uint64_t N,
                                int d, int w, double phi, uint64_t seed=12345ULL) {
    CountSketch cs(d,w,seed);
    for (const auto &p : counts) {
        const string &k = p.first;
        uint64_t c = p.second;
        for (uint64_t i = 0; i < c; ++i) cs.insert(k);
    }

    vector<uint64_t> abs_errs; abs_errs.reserve(counts.size());
    double sum_abs = 0.0;
    uint64_t max_abs = 0;
    uint64_t phi_thresh = (uint64_t)ceil(phi * (double)N);

    uint64_t tp=0, fp=0, fn=0;
    for (const auto &p : counts) {
        uint64_t truec = p.second;
        int64_t est = cs.estimate(p.first);
        if (est < 0) est = 0;
        uint64_t aerr = (uint64_t) llabs((int64_t)truec - est);
        abs_errs.push_back(aerr);
        sum_abs += (double)aerr;
        if (aerr > max_abs) max_abs = aerr;

        bool is_hh = (truec >= phi_thresh);
        bool pred = ((uint64_t)est >= phi_thresh);
        if (is_hh) { if (pred) tp++; else fn++; } else { if (pred) fp++; }
    }
    sort(abs_errs.begin(), abs_errs.end());
    double median = abs_errs.empty() ? 0.0 : (double)abs_errs[abs_errs.size()/2];
    double mean = counts.empty() ? 0.0 : sum_abs / (double)counts.size();

    EvalResult res;
    res.method = "CountSketch";
    res.mem_bytes = cs.memory_bytes();
    res.mean_abs_err = mean;
    res.median_abs_err = median;
    res.max_abs_err = (double)max_abs;
    
    res.tp = tp; res.fp = fp; res.fn = fn;
    return res;
}

EvalResult evaluate_towersketch(const vector<pair<string,uint64_t>> &counts, uint64_t N,
                                const vector<pair<int,int>> &levels, double phi, uint64_t seed=54321ULL) {
    TowerSketch ts(levels, seed);

    for (const auto &p : counts) {
        const string &k = p.first;
        uint64_t c = p.second;
        for (uint64_t i = 0; i < c; ++i) ts.insert(k);
    }

    vector<uint64_t> abs_errs; abs_errs.reserve(counts.size());
    double sum_abs = 0.0;
    uint64_t max_abs = 0;
    uint64_t phi_thresh = (uint64_t)ceil(phi * (double)N);

    uint64_t tp=0, fp=0, fn=0;
    for (const auto &p : counts) {
        uint64_t truec = p.second;
        uint64_t est = ts.estimate(p.first);
        uint64_t aerr = (uint64_t) llabs((int64_t)truec - (int64_t)est);
        abs_errs.push_back(aerr);
        sum_abs += (double)aerr;
        if (aerr > max_abs) max_abs = aerr;

        bool is_hh = (truec >= phi_thresh);
        bool pred = (est >= phi_thresh);
        if (is_hh) { if (pred) tp++; else fn++; } else { if (pred) fp++; }
    }
    sort(abs_errs.begin(), abs_errs.end());
    double median = abs_errs.empty() ? 0.0 : (double)abs_errs[abs_errs.size()/2];
    double mean = counts.empty() ? 0.0 : sum_abs / (double)counts.size();

    EvalResult res;
    res.method = "TowerSketch-CMCU";

    TowerSketch tmp(levels,seed);
    res.mem_bytes = tmp.memory_bytes(); 
    res.mean_abs_err = mean;
    res.median_abs_err = median;
    res.max_abs_err = (double)max_abs;
    res.tp = tp; res.fp = fp; res.fn = fn;
    return res;
}

int main(int argc, char **argv) {
    if (argc < 3) {
        cerr << "Uso: " << argv[0] << " counts.csv kmers.csv [phi]\n";
        return 1;
    }
    string gt = argv[1];
    string out = argv[2];
    double phi = 1e-6;
    if (argc >= 4) phi = stod(argv[3]);

    uint64_t N;
    auto counts = load_counts_and_N(gt, N);
    if (counts.empty()) {
        cerr << "Archivo vacio o no existe: " << gt << "\n";
        return 1;
    }
    cout << "Archivo: keys="<<counts.size()<<" total_kmers="<<N<<" phi="<<phi<<"\n";

    //CountSketch
    vector<int> ds = {3,5,7};
    vector<int> ws = {1<<7, 1<<8, 1<<9};

    //TowerSketch
    vector<vector<pair<int,int>>> ts_configs;
    ts_configs.push_back({{2,1<<7},{2,1<<8}});
    ts_configs.push_back({{2,1<<9},{2,1<<10}});
    ts_configs.push_back({{3,1<<7},{2,1<<9}});
    ts_configs.push_back({{2,1<<9}});

    ofstream fout(out);
    fout << "method,params,mem_bytes,mean_abs_err,median_abs_err,max_abs_err.tp,fp,fn\n";

    // Evaluar CountSketch
    for (int d : ds) {
        for (int w : ws) {
            cerr << "[CS] d="<<d<<" w="<<w<<" ...\n";
            EvalResult r = evaluate_countsketch(counts, N, d, w, phi);
            fout << r.method << ",\"d=" << d << ",w=" << w << "\","
                 << r.mem_bytes << "," << r.mean_abs_err << "," << r.median_abs_err << ","
                 << r.max_abs_err << "," << r.tp << "," << r.fp << "," << r.fn << "\n";
            fout.flush();
        }
    }

    // Evaluar TowerSketch
    for (const auto &cfg : ts_configs) {
        stringstream ss;
        ss << "levels=[";
        for (size_t i = 0; i < cfg.size(); ++i) {
            if (i) ss << ";";
            ss << cfg[i].first << "x" << cfg[i].second;
        }
        ss << "]";
        string params = ss.str();
        cerr << "[TS] " << params << " ...\n";
        EvalResult r = evaluate_towersketch(counts, N, cfg, phi);
        fout << r.method << ",\"" << params << "\"," << r.mem_bytes << "," << r.mean_abs_err << ","
             << r.median_abs_err << "," << r.max_abs_err << "," << r.tp << "," << r.fp << "," << r.fn << "\n";
        fout.flush();
    }

    fout.close();
    cout << "Listo. Resultados en " << out << "\n";
    return 0;
}
