// Compilar: 
// g++ -O3 -std=c++17 sketches_nuevo2.cpp -o sketches_nuevo2
// sketches_nuevo2 XYkmers.csv XYcounts.csv XYcounts_HH.csv results.csv 1e-6


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
#include <bits/stdc++.h>

using namespace std;

#include "utils.cpp"

static inline uint64_t splitmix64(uint64_t x) {
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    x = x ^ (x >> 31);
    return x;
}
static inline uint64_t hash_u64(const string &s, uint64_t salt=1469598103934665603ULL) {
    uint64_t h = std::hash<std::string>{}(s);
    h ^= salt + 0x9e3779b97f4a7c15ULL + (h<<6) + (h>>2);
    return splitmix64(h);
}
static inline uint64_t hash_u64_from_u64(uint64_t x, uint64_t salt) {
    return splitmix64(x ^ salt);
}

//CountSketch
class CountSketch {
public:
    CountSketch(int d=5, int w=(1<<14), uint64_t seed=12345ULL)
    : d_(d), w_(w), seed_(seed), table_((size_t)d*(size_t)w, 0) {
        if (d_<=0 || w_<=0) throw runtime_error("Invalid CS dims");
    }

    void insert_one(const string &kmer) {
        string cano = canonical(kmer);
        uint64_t key = hash_u64(cano, seed_);
        for (int i=0;i<d_;++i) {
            uint64_t h = hash_u64_from_u64(key, (uint64_t)i + seed_);
            uint32_t idx = (uint32_t)(h % (uint64_t)w_);
            int s = ((h >> 1) & 1) ? 1 : -1;
            size_t pos = (size_t)i * (size_t)w_ + idx;
            table_[pos] += (int64_t)s;
        }
    }

    void insert_repeated(const string &kmer, uint64_t count) {
        for (uint64_t i = 0; i < count; ++i) insert_one(kmer);
    }

    int64_t estimate(const string &kmer) const {
        string cano = canonical(kmer);
        uint64_t key = hash_u64(cano, seed_);
        vector<int64_t> est; est.reserve(d_);
        for (int i=0;i<d_;++i) {
            uint64_t h = hash_u64_from_u64(key, (uint64_t)i + seed_);
            uint32_t idx = (uint32_t)(h % (uint64_t)w_);
            int s = ((h >> 1) & 1) ? 1 : -1;
            size_t pos = (size_t)i * (size_t)w_ + idx;
            est.push_back((int64_t)s * table_[pos]);
        }
        size_t m = est.size()/2;
        nth_element(est.begin(), est.begin()+m, est.end());
        return est[m];
    }

    size_t memory_bytes() const {
        return table_.size() * sizeof(int64_t) + sizeof(*this);
    }

private:
    int d_, w_;
    uint64_t seed_;
    vector<int64_t> table_;
};

//TowerSketch (CMCU)
class TowerSketch {
public:
    struct Level {
        int r;
        int w;
        uint64_t seed;
        vector<uint64_t> table;
        vector<uint32_t> idxs;
        vector<uint64_t> vals;
    };

    TowerSketch(const vector<pair<int,int>> &levels_cfg, uint64_t seed=54321ULL)
    : seed_(seed) {
        if (levels_cfg.empty()) throw runtime_error("TowerSketch needs >=1 level");
        for (size_t i=0;i<levels_cfg.size();++i) {
            int r = levels_cfg[i].first;
            int w = levels_cfg[i].second;
            if (r<=0 || w<=0) throw runtime_error("invalid level cfg");
            Level L;
            L.r = r; L.w = w;
            L.seed = seed_ + (uint64_t)i * 0x9e3779b97f4a7c15ULL;
            L.table.assign((size_t)r * (size_t)w, 0ULL);
            L.idxs.assign(r,0);
            L.vals.assign(r,0);
            levels_.push_back(std::move(L));
        }
    }

    void insert_one(const string &kmer) {
        string cano = canonical(kmer);
        uint64_t key = hash_u64(cano, seed_);
        for (Level &L : levels_) {
            int r = L.r;
            int w = L.w;

            for (int i=0;i<r;++i) {
                uint64_t h = hash_u64_from_u64(key, L.seed + (uint64_t)i);
                uint32_t idx = (uint32_t)(h % (uint64_t)w);
                L.idxs[i] = idx;
                L.vals[i] = L.table[(size_t)i * (size_t)w + idx];
            }

            uint64_t m = UINT64_MAX;
            for (int i=0;i<r;++i) if (L.vals[i] < m) m = L.vals[i];

            for (int i=0;i<r;++i) {
                if (L.vals[i] == m) {
                    size_t pos = (size_t)i * (size_t)w + L.idxs[i];
                    if (L.table[pos] != UINT64_MAX) L.table[pos] += 1;
                }
            }
        }
    }

    void insert_batch(const string &kmer, uint64_t delta) {
        if (delta == 0) return;
        string cano = canonical(kmer);
        uint64_t key = hash_u64(cano, seed_);
        for (Level &L : levels_) {
            int r = L.r;
            int w = L.w;
            for (int i=0;i<r;++i) {
                uint64_t h = hash_u64_from_u64(key, L.seed + (uint64_t)i);
                uint32_t idx = (uint32_t)(h % (uint64_t)w);
                L.idxs[i] = idx;
                L.vals[i] = L.table[(size_t)i * (size_t)w + idx];
            }
            uint64_t m = UINT64_MAX;
            for (int i=0;i<r;++i) if (L.vals[i] < m) m = L.vals[i];
            for (int i=0;i<r;++i) {
                if (L.vals[i] == m) {
                    size_t pos = (size_t)i * (size_t)w + L.idxs[i];
                    if (L.table[pos] > numeric_limits<uint64_t>::max() - delta) L.table[pos] = numeric_limits<uint64_t>::max();
                    else L.table[pos] += delta;
                }
            }
        }
    }

    uint64_t estimate(const string &kmer) const {
        string cano = canonical(kmer);
        uint64_t key = hash_u64(cano, seed_);
        uint64_t best = UINT64_MAX;
        for (const Level &L : levels_) {
            uint64_t m = UINT64_MAX;
            for (int i=0;i<L.r;++i) {
                uint64_t h = hash_u64_from_u64(key, L.seed + (uint64_t)i);
                uint32_t idx = (uint32_t)(h % (uint64_t)L.w);
                uint64_t v = L.table[(size_t)i * (size_t)L.w + idx];
                if (v < m) m = v;
            }
            if (m < best) best = m;
        }
        return best == UINT64_MAX ? 0 : best;
    }

    size_t memory_bytes() const {
        size_t s = sizeof(*this);
        for (const Level &L : levels_) s += L.table.size() * sizeof(uint64_t);
        return s;
    }

private:
    uint64_t seed_;
    vector<Level> levels_;
};

static vector<string> read_stream_kmers(const string &path) {
    vector<string> stream;
    ifstream in(path);
    if (!in.is_open()) throw runtime_error("cannot open stream file: "+path);
    string line;
    while (getline(in, line)) {
        if (line.empty()) continue;
        // trim
        while (!line.empty() && isspace((unsigned char)line.back())) line.pop_back();
        size_t comma = line.find(',');
        string token = (comma==string::npos) ? line : line.substr(0,comma);
        if (token.empty()) continue;
        stream.push_back(token);
    }
    in.close();
    return stream;
}

static vector<pair<string,uint64_t>> read_counts_csv(const string &path, uint64_t &N) {
    vector<pair<string,uint64_t>> counts;
    ifstream in(path);
    if (!in.is_open()) throw runtime_error("cannot open counts file: "+path);
    string line;

    if (!getline(in,line)) return counts;
    bool has_header = false;
    string tmp = line;
    
    string low = tmp;
    transform(low.begin(), low.end(), low.begin(), ::tolower);
    if (low.find("kmer")!=string::npos || low.find("count")!=string::npos || low.find(',')==string::npos) {
        // if there's a comma but header may include, check next char
        if (low.find("kmer")!=string::npos || low.find("count")!=string::npos) has_header = true;
    }
    if (!has_header) {
        istringstream ss(line);
        string a,b;
        if (getline(ss,a,',')) {
            if (ss >> b) {
                string k = a;
                uint64_t c = stoull(b);
                counts.emplace_back(k,c);
            }
        }
    }
    while (getline(in,line)) {
        if (line.empty()) continue;
        istringstream ss(line);
        string k; uint64_t c;
        if (getline(ss,k,',')) {
            ss >> c;
            counts.emplace_back(k,c);
        }
    }
    in.close();
    N = 0;
    for (auto &p : counts) {
        p.first = canonical(p.first);
        N += p.second;
    }
    return counts;
}

static unordered_set<string> read_hh_csv(const string &path) {
    unordered_set<string> hh;
    ifstream in(path);
    if (!in.is_open()) return hh;
    string line;

    getline(in,line);
    while (getline(in,line)) {
        if (line.empty()) continue;
        istringstream ss(line);
        string k;
        if (getline(ss,k,',')) {
            hh.insert(canonical(k));
        }
    }
    in.close();
    return hh;
}

struct EvalResult {
    string method;
    string params;
    size_t mem_bytes;
    double mean_abs;
    double median_abs;
    double max_abs;
    double precision;
    double recall;
    double f1;
    uint64_t tp, fp, fn;
};

static EvalResult evaluate_sketch_from_counts(const vector<pair<string,uint64_t>> &counts,
                                              uint64_t N,
                                              function<uint64_t(const string&)> estimator,
                                              const unordered_set<string> &true_hh,
                                              double phi)
{
    vector<uint64_t> abs_errs; abs_errs.reserve(counts.size());
    double sum_abs = 0.0;
    uint64_t max_abs = 0;
    uint64_t phi_thresh = (uint64_t)ceil(phi * (double)N);

    uint64_t tp=0, fp=0, fn=0;
    for (const auto &p : counts) {
        const string &k = p.first;
        uint64_t truec = p.second;
        uint64_t est = estimator(k);
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
    double precision = (tp + fp) == 0 ? 0.0 : (double)tp / (double)(tp + fp);
    double recall = (tp + fn) == 0 ? 0.0 : (double)tp / (double)(tp + fn);
    double f1 = (precision + recall) == 0 ? 0.0 : 2.0 * precision * recall / (precision + recall);

    EvalResult r;
    r.mean_abs = mean;
    r.median_abs = median;
    r.max_abs = (double)max_abs;
    r.precision = precision; r.recall = recall; r.f1 = f1;
    r.tp = tp; r.fp = fp; r.fn = fn;
    return r;
}

int main(int argc, char **argv) {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    if (argc < 5) {
        cerr << "Usage: " << argv[0] << " <stream_kmers.txt> <counts.csv> <hh.csv> <out_results.csv> [phi]\n";
        return 1;
    }
    string stream_file = argv[1];
    string counts_file = argv[2];
    string hh_file = argv[3];
    string out_file = argv[4];
    double phi = 1e-3;
    if (argc >= 6) phi = stod(argv[5]);

    cerr << "[1/4] Cargando ground-truth counts...\n";
    uint64_t N;
    auto counts = read_counts_csv(counts_file, N);
    if (counts.empty()) { cerr<<"counts vacio\n"; return 1; }
    unordered_map<string,uint64_t> counts_map;
    for (auto &p : counts) counts_map[p.first] = p.second;

    cerr << "[2/4] Cargando stream de k-mers (may ser grande)...\n";
    vector<string> stream = read_stream_kmers(stream_file);
    if (stream.empty()) { cerr << "stream vacio\n"; return 1; }

    cerr << "[3/4] Cargando heavy-hitters (si existe)...\n";
    unordered_set<string> hh_set = read_hh_csv(hh_file);

    cerr << "[4/4] Preparando configuraciones y calibrando...\n";
    // CountSketch configs
    vector<pair<int,int>> cs_configs = {
        {3, 128}, {3, 256}, {3, 512},
        {5, 1024}, {5, 2048}, {5, 4096},
        {7, 8192}
    };

    // TowerSketch configs
    vector<vector<pair<int,int>>> ts_configs;
    ts_configs.push_back({{2,128},{2,512}});   // pequeño + mediano
    ts_configs.push_back({{2,256},{2,1024}});
    ts_configs.push_back({{3,128},{2,512}});
    ts_configs.push_back({{2,512}});

    ofstream out(out_file);
    if (!out.is_open()) { cerr<<"No puedo abrir "<<out_file<<"\n"; return 1; }
    out << "method,params,mem_bytes,mean_abs,median_abs,max_abs,precision,recall,f1,tp,fp,fn\n";

    // Evaluate CountSketch
    for (auto &cfg : cs_configs) {
        int d = cfg.first;
        int w = cfg.second;
        cerr << "[CS] d="<<d<< " w="<<w << " ... feeding stream\n";
        CountSketch cs(d,w);
       
        for (const string &k : stream) cs.insert_one(k);
        cerr << "  estimating ...\n";
        auto est = [&cs](const string &k){ return (uint64_t)cs.estimate(k); };
        EvalResult r = evaluate_sketch_from_counts(counts, N, est, hh_set, phi);
        r.method = "CountSketch";
        {
            stringstream ss; ss << "d="<<d<<",w="<<w;
            r.params = ss.str();
        }
        r.mem_bytes = cs.memory_bytes();
        out << r.method << "," << "\"" << r.params << "\"," << r.mem_bytes << "," 
            << r.mean_abs << "," << r.median_abs << "," << r.max_abs << ","
            << r.precision << "," << r.recall << "," << r.f1 << ","
            << r.tp << "," << r.fp << "," << r.fn << "\n";
        out.flush();
    }

    for (auto &cfg : ts_configs) {

        cerr << "[TS] levels=[";
        for (size_t i=0;i<cfg.size();++i) {
            if (i) cerr << ";";
            cerr << cfg[i].first << "x" << cfg[i].second;
        }
        cerr << "] ... feeding stream\n";
        TowerSketch ts(cfg);

        for (const string &k : stream) ts.insert_one(k);
        cerr << "  estimating ...\n";
        auto est = [&ts](const string &k){ return ts.estimate(k); };
        EvalResult r = evaluate_sketch_from_counts(counts, N, est, hh_set, phi);
        r.method = "TowerSketch-CMCU";
        {
            stringstream ss; ss << "levels=[";
            for (size_t i=0;i<cfg.size();++i) {
                if (i) ss << ";";
                ss << cfg[i].first << "x" << cfg[i].second;
            }
            ss << "]";
            r.params = ss.str();
        }
        r.mem_bytes = ts.memory_bytes();
        out << r.method << "," << "\"" << r.params << "\"," << r.mem_bytes << "," 
            << r.mean_abs << "," << r.median_abs << "," << r.max_abs << ","
            << r.precision << "," << r.recall << "," << r.f1 << ","
            << r.tp << "," << r.fp << "," << r.fn << "\n";
        out.flush();
    }

    out.close();
    cerr << "Calibracion completada. Resultados en " << out_file << "\n";
    return 0;
}
