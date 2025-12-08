// holt_ga_modified.cpp
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <cctype>
#include <sstream>
#include <iomanip>
#include <random>
#include <unordered_map>

using namespace std;
using vec = vector<double>;

struct DateVal {
    int y, m, d;
    string datestr;
    double val;
};

/* ======================= CSV HELPERS ======================= */

static inline string trim(const string &s){
    size_t a = 0, b = s.size();
    while (a < b && isspace((unsigned char)s[a])) a++;
    while (b > a && isspace((unsigned char)s[b-1])) b--;
    return s.substr(a, b - a);
}

static vector<string> split_tokens(const string &line, char sep=','){
    vector<string> cols;
    string cur;
    bool inq = false;
    for (char c : line) {
        if (c == '"') { inq = !inq; continue; }
        if (c == sep && !inq) {
            cols.push_back(cur);
            cur.clear();
        } else cur.push_back(c);
    }
    cols.push_back(cur);
    return cols;
}

// parse a few common date formats (MM/DD/YYYY, M/D/YYYY, YYYY-MM-DD, etc.)
static bool parse_date_mdY(const string &s, int &y, int &m, int &d){
    string t;
    for (char c : s) if (c != '\"') t.push_back(c);
    t = trim(t);
    if (t.empty()) return false;
    for (char &c : t) if (c == '-') c = '/';

    vector<string> parts;
    string cur;
    for (char c : t){
        if (c == '/') { parts.push_back(cur); cur.clear(); }
        else cur.push_back(c);
    }
    parts.push_back(cur);

    if (parts.size() == 3) {
        int yi = -1;
        for (int i = 0; i < 3; ++i) {
            string p = trim(parts[i]);
            if (p.size() == 4) { yi = i; break; }
            try {
                int v = stoi(p);
                if (v > 31) { yi = i; break; }
            } catch (...) {}
        }
        if (yi == -1) yi = 2;
        try {
            y = stoi(trim(parts[yi]));
            vector<int> other;
            for (int i = 0; i < 3; ++i) if (i != yi) other.push_back(stoi(trim(parts[i])));
            if (yi == 0) { m = other[0]; d = other[1]; }
            else { m = other[0]; d = other[1]; }
            return true;
        } catch (...) { return false; }
    }
    return false;
}

vector<DateVal> load_csv_dates_values(const string &filename) {
    ifstream f(filename);
    if (!f.is_open()) {
        cerr << "Cannot open file: " << filename << "\n";
        exit(1);
    }

    string line;
    vector<string> header;
    vector<vector<string>> rows;
    bool first = true;

    while (getline(f, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (first) {
            header = split_tokens(line);
            first = false;
            continue;
        }
        rows.push_back(split_tokens(line));
    }
    f.close();

    auto lower = [&](string s){
        for (char &c : s) c = (char)tolower(c);
        return s;
    };

    int dateIdx = -1, valIdx = -1;
    for (int i = 0; i < (int)header.size(); ++i) {
        string h = lower(trim(header[i]));
        if (h.find("date") != string::npos) dateIdx = i;
        if (h.find("close") != string::npos ||
            h.find("price") != string::npos ||
            h.find("last") != string::npos) valIdx = i;
    }
    if (dateIdx == -1) dateIdx = 0;
    if (valIdx == -1) valIdx = 1;

    vector<DateVal> out;
    for (auto &r : rows) {
        if (dateIdx >= (int)r.size() || valIdx >= (int)r.size()) continue;
        string sdate = trim(r[dateIdx]);
        string sval  = r[valIdx];
        sval.erase(remove(sval.begin(), sval.end(), ','), sval.end());
        if (sval.empty()) continue;

        double v;
        try { v = stod(sval); } catch (...) { continue; }

        int y=0, m=0, d=0;
        if (!parse_date_mdY(sdate, y, m, d)) continue;

        out.push_back({y,m,d,sdate,v});
    }

    sort(out.begin(), out.end(), [](const DateVal &a, const DateVal &b){
        if (a.y != b.y) return a.y < b.y;
        if (a.m != b.m) return a.m < b.m;
        return a.d < b.d;
    });

    return out;
}

/* ======================= HOLT-WINTERS ======================= */

struct HWResult {
    vec fitted;
    vec forecast;
    double mse;
};

HWResult holt_winters_additive(const vec& data, int m, int horizon,
                               double alpha, double beta, double gamma,
                               bool verbose=false)
{
    int n = (int)data.size();
    if (n < 2 * m) return {{}, {}, 1e18};

    // Initial level: mean of first season
    double s1 = accumulate(data.begin(), data.begin() + m, 0.0) / m;
    // Next season mean
    double s2 = accumulate(data.begin() + m, data.begin() + 2*m, 0.0) / m;
    double level = s1;
    double trend = (s2 - s1) / static_cast<double>(m);

    int seasons = n / m;

    vec season(m, 0.0);
    vector<double> season_avg(seasons, 0.0);
    for (int s = 0; s < seasons; ++s) {
        double ss = 0;
        for (int i = 0; i < m; ++i) ss += data[s*m + i];
        season_avg[s] = ss / m;
    }
    for (int i = 0; i < m; ++i) {
        double sum = 0;
        int cnt = 0;
        for (int s = 0; s < seasons; ++s) {
            int idx = s*m + i;
            if (idx < n) {
                sum += data[idx] - season_avg[s];
                cnt++;
            }
        }
        season[i] = (cnt>0 ? sum / cnt : 0.0);
    }

    if (verbose) {
        cout << fixed << setprecision(6);
        cout << "\n[Additive] Initial level (L0) = " << level << "\n";
        cout << "[Additive] Initial trend (T0) = " << trend << "  (note: used (s2 - s1)/m)\n";
        cout << "[Additive] Initial seasonals S[0.." << m-1 << "] =";
        for (int i = 0; i < m; ++i) cout << " " << season[i];
        cout << "\n\n";
    }

    vec fitted(n);
    vec L{level}, T{trend};

    for (int t = 0; t < n; ++t) {
        int si = t % m;
        double forecast_t = L.back() + T.back() + season[si];
        fitted[t] = forecast_t;

        double newL = alpha * (data[t] - season[si]) + (1.0 - alpha) * (L.back() + T.back());
        double newT = beta * (newL - L.back()) + (1.0 - beta) * T.back();
        double newS = gamma * (data[t] - newL) + (1.0 - gamma) * season[si];

        L.push_back(newL);
        T.push_back(newT);
        season[si] = newS;
    }

    if (verbose) {
        cout << "[Additive] Final level (L_n) = " << L.back() << "\n";
        cout << "[Additive] Final trend (T_n) = " << T.back() << "\n";
        cout << "[Additive] Final seasonals S[0.." << m-1 << "] =";
        for (int i = 0; i < m; ++i) cout << " " << season[i];
        cout << "\n\n";
    }

    vec forecast(horizon);
    for (int k = 1; k <= horizon; ++k)
        forecast[k-1] = L.back() + k * T.back() + season[(n + k - 1) % m];

    double sse = 0;
    for (int i = 0; i < n; ++i) {
        double e = data[i] - fitted[i];
        sse += e*e;
    }
    return {fitted, forecast, sse / n};
}

HWResult holt_winters_multiplicative(const vec& data, int m, int horizon,
                                     double alpha, double beta, double gamma,
                                     bool verbose=false)
{
    int n = (int)data.size();
    if (n < 2 * m) return {{}, {}, 1e18};

    // Initial level: mean first season
    double s1 = accumulate(data.begin(), data.begin() + m, 0.0) / m;
    double s2 = accumulate(data.begin() + m, data.begin() + 2*m, 0.0) / m;
    double level = s1;
    double trend = (s2 - s1) / static_cast<double>(m);

    int seasons = n / m;

    vector<double> season_avg(seasons, 0.0);
    for (int s = 0; s < seasons; ++s) {
        double ss = 0;
        for (int i = 0; i < m; ++i) ss += data[s*m + i];
        season_avg[s] = ss / m;
    }

    vec season(m, 1.0);
    for (int i = 0; i < m; ++i) {
        double sum = 0;
        int cnt = 0;
        for (int s = 0; s < seasons; ++s) {
            int idx = s*m + i;
            if (idx < n) {
                sum += data[idx] / season_avg[s];
                cnt++;
            }
        }
        season[i] = (cnt>0 ? sum / cnt : 1.0);
    }

    if (verbose) {
        cout << fixed << setprecision(6);
        cout << "\n[Multiplicative] Initial level (L0) = " << level << "\n";
        cout << "[Multiplicative] Initial trend (T0) = " << trend << "  (note: used (s2 - s1)/m)\n";
        cout << "[Multiplicative] Initial seasonals S[0.." << m-1 << "] =";
        for (int i = 0; i < m; ++i) cout << " " << season[i];
        cout << "\n\n";
    }

    vec fitted(n);
    vec L{level}, T{trend};

    for (int t = 0; t < n; ++t) {
        int si = t % m;
        double forecast_t = (L.back() + T.back()) * season[si];
        fitted[t] = forecast_t;

        double denom_season = max(1e-12, season[si]);
        double newL = alpha * (data[t] / denom_season) + (1.0 - alpha) * (L.back() + T.back());
        double newT = beta * (newL - L.back()) + (1.0 - beta) * T.back();
        double denom_level = max(1e-12, newL);
        double newS = gamma * (data[t] / denom_level) + (1.0 - gamma) * season[si];

        L.push_back(newL);
        T.push_back(newT);
        season[si] = newS;
    }

    if (verbose) {
        cout << "[Multiplicative] Final level (L_n) = " << L.back() << "\n";
        cout << "[Multiplicative] Final trend (T_n) = " << T.back() << "\n";
        cout << "[Multiplicative] Final seasonals S[0.." << m-1 << "] =";
        for (int i = 0; i < m; ++i) cout << " " << season[i];
        cout << "\n\n";
    }

    vec forecast(horizon);
    for (int k = 1; k <= horizon; ++k)
        forecast[k-1] = (L.back() + k * T.back()) * season[(n + k - 1) % m];

    double sse = 0;
    for (int i = 0; i < n; ++i) {
        double e = data[i] - fitted[i];
        sse += e*e;
    }
    return {fitted, forecast, sse / n};
}

/* ======================= GENETIC ALGORITHM + CACHE ======================= */

struct HWParams { double a,b,g; };

static std::mt19937 rng_engine(123456789); // fixed seed for reproducibility

// small cache to avoid repeated evaluations of the same quantized (a,b,g)
static unordered_map<string,double> fitness_cache;

// helper to create a quantized key for params
static string params_key(const HWParams &p) {
    std::ostringstream ss;
    ss << fixed << setprecision(6) << p.a << "_" << p.b << "_" << p.g;
    return ss.str();
}

double evaluate_hw_mse(const vec &data, int season_len, const HWParams &p, bool multiplicative) {
    string key = params_key(p);
    auto it = fitness_cache.find(key);
    if (it != fitness_cache.end()) return it->second;

    HWResult res;
    if (!multiplicative) res = holt_winters_additive(data, season_len, 0, p.a, p.b, p.g, /*verbose=*/false);
    else res = holt_winters_multiplicative(data, season_len, 0, p.a, p.b, p.g, /*verbose=*/false);
    fitness_cache[key] = res.mse;
    return res.mse;
}

double uniform_real(double lo, double hi) {
    std::uniform_real_distribution<double> d(lo,hi);
    return d(rng_engine);
}
int uniform_int(int lo, int hi) {
    std::uniform_int_distribution<int> d(lo,hi);
    return d(rng_engine);
}
double clamp01(double x){ if (x<0) return 0; if (x>1) return 1; return x; }

struct Individual { HWParams p; double fitness; };

HWParams ga_optimize_hw(const vec &data, int season_len, bool multiplicative,
                        int pop_size=80, int generations=100, double crossover_rate=0.8, double mutation_rate=0.2)
{
    auto make_random_ind = [&](){
        Individual ind;
        ind.p.a = uniform_real(0.0, 1.0);
        ind.p.b = uniform_real(0.0, 1.0);
        ind.p.g = uniform_real(0.0, 1.0);
        ind.fitness = 1e18;
        return ind;
    };

    vector<Individual> pop(pop_size);
    for (int i=0;i<pop_size;++i) pop[i] = make_random_ind();
    for (int i=0;i<pop_size;++i) pop[i].fitness = evaluate_hw_mse(data, season_len, pop[i].p, multiplicative);

    for (int gen=0; gen<generations; ++gen) {
        sort(pop.begin(), pop.end(), [](const Individual &x, const Individual &y){ return x.fitness < y.fitness; });

        // optional logging: only GA summary (no HW verbose)
        if (gen % 10 == 0) {
            cout << "[GA] gen="<<gen<<" best MSE="<<pop[0].fitness
                 << " (a="<<pop[0].p.a<<" b="<<pop[0].p.b<<" g="<<pop[0].p.g<<")\n";
        }

        vector<Individual> newpop;
        // elitism: keep top 2 (if pop_size >= 2)
        if (!pop.empty()) newpop.push_back(pop[0]);
        if (pop_size > 1) newpop.push_back(pop[1]);

        while ((int)newpop.size() < pop_size) {
            // tournament selection
            auto tournament = [&](int k=3){
                Individual best = pop[uniform_int(0, pop_size-1)];
                for (int t=1;t<k;++t) {
                    int idx = uniform_int(0, pop_size-1);
                    if (pop[idx].fitness < best.fitness) best = pop[idx];
                }
                return best;
            };

            Individual p1 = tournament();
            Individual p2 = tournament();

            Individual child;
            if (uniform_real(0.0,1.0) < crossover_rate) {
                double w = uniform_real(0.0,1.0);
                child.p.a = p1.p.a * w + p2.p.a * (1-w);
                child.p.b = p1.p.b * w + p2.p.b * (1-w);
                child.p.g = p1.p.g * w + p2.p.g * (1-w);
            } else {
                child = (uniform_real(0.0,1.0) < 0.5 ? p1 : p2);
            }

            // mutations
            if (uniform_real(0.0,1.0) < mutation_rate) {
                normal_distribution<double> nd(0.0, 0.05);
                child.p.a = clamp01(child.p.a + nd(rng_engine));
            }
            if (uniform_real(0.0,1.0) < mutation_rate) {
                normal_distribution<double> nd(0.0, 0.05);
                child.p.b = clamp01(child.p.b + nd(rng_engine));
            }
            if (uniform_real(0.0,1.0) < mutation_rate) {
                normal_distribution<double> nd(0.0, 0.05);
                child.p.g = clamp01(child.p.g + nd(rng_engine));
            }

            child.fitness = evaluate_hw_mse(data, season_len, child.p, multiplicative);
            newpop.push_back(child);
        }

        pop.swap(newpop);
    }

    sort(pop.begin(), pop.end(), [](const Individual &x, const Individual &y){ return x.fitness < y.fitness; });
    cout << "[GA] finished best MSE="<<pop[0].fitness
         << " (a="<<pop[0].p.a<<" b="<<pop[0].p.b<<" g="<<pop[0].p.g<<")\n";
    return pop[0].p;
}

/* ======================= MAIN ======================= */

int main() {
    // <-- change the CSV path to your local file -->
    string csvfile = "C:\\Users\\rakkh\\OneDrive\\Desktop\\RP\\PSEi Composite Historical Data (1).csv";

    cout << "Loading CSV (date + value)...\n";
    auto dv = load_csv_dates_values(csvfile);
    cout << "Rows read: " << dv.size() << "\n";

    // split train (2010-2021) and future (>=2022)
    vector<DateVal> train_rows, future_rows;
    for (auto &x : dv) {
        if (x.y >= 2010 && x.y <= 2021) train_rows.push_back(x);
        else if (x.y >= 2022) future_rows.push_back(x);
    }
    cout << "Training rows: " << train_rows.size() << ", Future rows: " << future_rows.size() << "\n";

    if (train_rows.size() < 24) {
        cerr << "Not enough training data (need >= 24 months). Exiting.\n";
        return 1;
    }

    // build training values vector (raw prices, NO normalization)
    vec train_vals;
    for (auto &r : train_rows) train_vals.push_back(r.val);

    // GA and Holt-Winters parameters
    int season_len = 12;
    int horizon = 24;

    cout << "\nOptimizing Holt-Winters parameters using GA (on ORIGINAL price scale)...\n";
    // Reduced defaults for faster testing. Increase for final runs.
    int ga_pop = 30;
    int ga_gens = 40;

    cout << "Optimizing ADDITIVE model (GA)...\n";
    HWParams best_add = ga_optimize_hw(train_vals, season_len, false, ga_pop, ga_gens);

    cout << "Optimizing MULTIPLICATIVE model (GA)...\n";
    HWParams best_mul = ga_optimize_hw(train_vals, season_len, true, ga_pop, ga_gens);

    cout << "\nRunning Holt-Winters (Additive) with GA-found params on training (original scale)...\n";
    HWResult add = holt_winters_additive(train_vals, season_len, horizon, best_add.a, best_add.b, best_add.g, true);
    cout << "Additive in-sample MSE = " << add.mse << "\n";

    cout << "\nRunning Holt-Winters (Multiplicative) with GA-found params on training (original scale)...\n";
    HWResult mul = holt_winters_multiplicative(train_vals, season_len, horizon, best_mul.a, best_mul.b, best_mul.g, true);
    cout << "Multiplicative in-sample MSE = " << mul.mse << "\n";

    bool pick_add = (add.mse < mul.mse);
    cout << "\nChosen model on training MSE: " << (pick_add ? "Additive" : "Multiplicative") << "\n";

    vec chosen_forecast = pick_add ? add.forecast : mul.forecast;

    // save forecasts with actuals (real scale) -> forecast_REAL.csv
    ofstream fo("forecast_REAL.csv");
    fo << "Index,Forecast,Actual\n";
    DateVal last_train = train_rows.back();
    int yy = last_train.y;
    int mm = last_train.m;
    // assemble future values vector
    vec future_vals;
    for (auto &r : future_rows) future_vals.push_back(r.val);

    for (int i = 0; i < horizon; ++i) {
        mm++;
        if (mm > 12) { mm = 1; yy++; }
        double f = (i < (int)chosen_forecast.size() ? chosen_forecast[i] : NAN);
        double a = (i < (int)future_vals.size() ? future_vals[i] : NAN);
        fo << i+1 << "," << f << "," << a << "\n";
    }
    fo.close();

    // save chosen params and other model params for record
    ofstream pf("forecast_PARAMS.csv");
    pf << "model,alpha,beta,gamma,in_sample_mse\n";
    pf << "additive," << best_add.a << "," << best_add.b << "," << best_add.g << "," << add.mse << "\n";
    pf << "multiplicative," << best_mul.a << "," << best_mul.b << "," << best_mul.g << "," << mul.mse << "\n";
    pf.close();

    // evaluation on overlapping portion
    int n = min(horizon, (int)future_vals.size());
    if (n == 0) {
        cout << "No future actuals available to evaluate forecasts.\n";
        cout << "Done. Two CSV files produced: forecast_REAL.csv and forecast_PARAMS.csv\n";
        return 0;
    }

    double sse_real = 0.0, sae_real = 0.0;
    for (int i = 0; i < n; ++i) {
        double f_real = chosen_forecast[i];
        double a_real = future_vals[i];
        double e = a_real - f_real;
        sse_real += e*e;
        sae_real += fabs(e);
    }

    cout << "\nEVALUATION (overlap " << n << " months):\n";
    cout << " - REAL SCALE: MSE = " << (sse_real / n) << " , MAE = " << (sae_real / n) << "\n";

    cout << "\nFirst 8 forecasts vs actuals (real-scale):\n";
    cout << "Idx\tForecast\tActual\n";
    for (int i = 0; i < min(n, 8); ++i) {
        cout << (i+1) << "\t" << chosen_forecast[i] << "\t" << future_vals[i] << "\n";
    }

    cout << "\nDone. Two CSV files produced: forecast_REAL.csv and forecast_PARAMS.csv\n";
    cout << "Check the printed initial/final level/trend/season values above to compare with theoretical definitions.\n";

    return 0;
}
