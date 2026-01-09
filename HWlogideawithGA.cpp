// holt_2017_2021_log_vs_orig_add_vs_mul_GA_explained.cpp
//
// Goal:
// 1) Load monthly PSEi data from CSV.
// 2) Use 2017–2021 as training, 2022 as test.
// 3) Use Genetic Algorithm (GA) to find optimal Holt–Winters parameters
//    alpha, beta, gamma in [0,1] for both additive and multiplicative models.
// 4) Compare forecasts on original scale vs log scale.
// 5) Output a CSV with forecasts and errors.
//
// Key outputs:
// - GA-optimized (alpha, beta, gamma)
// - Best model type (Additive or Multiplicative) on original scale
// - CSV: forecast_2017_2021_vs_2022_add_vs_mul_log_vs_orig_GA.csv

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
#include <limits>
#include <cstdlib>

using namespace std;

// Convenience alias: vec = vector of double
using vec = vector<double>;

// Struct to store a single date and its PSEi value
struct DateVal {
    int    y;       // year (e.g., 2017)
    int    m;       // month (1–12)
    int    d;       // day   (usually 1 for monthly data)
    string datestr; // original date string from CSV
    double val;     // closing price / index value
};

/* ======================= CSV HELPERS ======================= */

// Remove leading and trailing whitespace from a string
static inline string trim(const string &s){
    size_t a = 0;
    size_t b = s.size();

    // move a forward until first non-space
    while (a < b && isspace((unsigned char)s[a])) a++;
    // move b backward until last non-space
    while (b > a && isspace((unsigned char)s[b-1])) b--;

    return s.substr(a, b - a);
}

// Split a CSV line into columns, respecting quotes
static vector<string> split_tokens(const string &line, char sep=','){
    vector<string> cols;
    string cur;
    bool inq = false;  // true if we are inside quotes

    for (char c : line) {
        if (c == '"') {
            // toggle quote state and skip the quote char itself
            inq = !inq;
            continue;
        }
        if (c == sep && !inq) {
            // separator outside quotes -> end of column
            cols.push_back(cur);
            cur.clear();
        } else {
            // normal character
            cur.push_back(c);
        }
    }
    cols.push_back(cur);  // last column
    return cols;
}

// Parse a date string of various forms into (y,m,d)
// Example formats: "2020-01-01", "01/01/2020", "2020/1/1"
static bool parse_date_mdY(const string &s, int &y, int &m, int &d){
    // Remove quotes from date
    string t;
    for (char c : s) if (c != '\"') t.push_back(c);
    t = trim(t);
    if (t.empty()) return false;

    // Normalize separators to '/'
    for (char &c : t) if (c == '-') c = '/';

    // Split by '/'
    vector<string> parts;
    string cur;
    for (char c : t){
        if (c == '/') {
            parts.push_back(cur);
            cur.clear();
        } else {
            cur.push_back(c);
        }
    }
    parts.push_back(cur);

    if (parts.size() == 3) {
        // Try to detect which part is the year
        int yi = -1;  // index of year component
        for (int i = 0; i < 3; ++i) {
            string p = trim(parts[i]);
            // Heuristic: 4-digit or >31 -> treat as year
            if (p.size() == 4) { yi = i; break; }
            try {
                int v = stoi(p);
                if (v > 31) { yi = i; break; }
            } catch(...) {}
        }
        if (yi == -1) yi = 2; // fallback: assume last is year

        try {
            y = stoi(trim(parts[yi]));
            vector<int> other;
            for (int i = 0; i < 3; ++i) {
                if (i != yi)
                    other.push_back(stoi(trim(parts[i])));
            }
            // For monthly data, we only care about year and month.
            // Here we treat the first of the remaining as month, second as day.
            if (yi == 0) { m = other[0]; d = other[1]; }
            else         { m = other[0]; d = other[1]; }
            return true;
        } catch (...) {
            return false;
        }
    }
    return false;
}

// Load CSV with at least date and close/price column
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

    // Read all lines
    while (getline(f, line)) {
        // Handle Windows CRLF line ending
        if (!line.empty() && line.back()=='\r') line.pop_back();

        if (first) {
            // First line is header
            header = split_tokens(line);
            first  = false;
            continue;
        }
        rows.push_back(split_tokens(line));
    }
    f.close();

    // helper to lowercase a string
    auto lower = [&](string s){
        for (char &c : s) c = (char)tolower(c);
        return s;
    };

    // Detect which columns correspond to date and value
    int dateIdx=-1, valIdx=-1;
    for (int i = 0; i < (int)header.size(); ++i) {
        string h = lower(trim(header[i]));
        if (h.find("date")   != string::npos) dateIdx = i;
        if (h.find("close")  != string::npos ||
            h.find("price")  != string::npos ||
            h.find("last")   != string::npos) valIdx = i;
    }
    // Fallback if detection fails
    if (dateIdx == -1) dateIdx = 0;
    if (valIdx  == -1) valIdx  = 1;

    vector<DateVal> out;
    for (auto &r : rows) {
        if (dateIdx >= (int)r.size() || valIdx >= (int)r.size()) continue;

        string sdate = trim(r[dateIdx]);  // date string
        string sval  = r[valIdx];         // price string (may have commas)

        // Remove commas from price, e.g. "6,500.23" -> "6500.23"
        sval.erase(remove(sval.begin(), sval.end(), ','), sval.end());
        if (sval.empty()) continue;

        double v;
        try {
            v = stod(sval);
        } catch (...) {
            continue;
        }

        int y=0,m=0,d=0;
        if (!parse_date_mdY(sdate,y,m,d)) continue;

        out.push_back({y,m,d,sdate,v});
    }

    // Sort by date ascending
    sort(out.begin(), out.end(), [](const DateVal &a, const DateVal &b){
        if (a.y != b.y) return a.y < b.y;
        if (a.m != b.m) return a.m < b.m;
        return a.d < b.d;
    });

    return out;
}

/* ======================= HOLT-WINTERS IMPLEMENTATIONS ======================= */

// Result struct for Holt–Winters
struct HWResult {
    vec   fitted;   // in-sample fitted values
    vec   forecast; // out-of-sample forecasts
    double mse;     // mean squared error on in-sample data
};

// Additive Holt–Winters implementation
HWResult holt_winters_additive(
    const vec& data,   // time series
    int        m,      // season length (e.g., 12 for monthly)
    int        horizon,// forecast horizon (how many steps ahead)
    double     alpha,  // smoothing for level
    double     beta,   // smoothing for trend
    double     gamma,  // smoothing for seasonality
    bool       verbose // print initial parameters?
){
    int n = (int)data.size();
    if (n < 2 * m) return {{}, {}, 1e18}; // not enough data

    // Initial level: average of first season
    double s1 = accumulate(data.begin(), data.begin() + m, 0.0) / m;
    double level = s1;

    // Initial trend: average change between first two seasons
    double trend_sum = 0;
    for (int i = 0; i < m; ++i)
        trend_sum += (data[m + i] - data[i]);
    double trend = (trend_sum / (double)m) / (double)m;

    // Initialize seasonal components (additive)
    int seasons = n / m;
    vec season(m, 0.0);
    vector<double> season_avg(seasons, 0.0);

    // Season averages
    for (int s = 0; s < seasons; ++s) {
        double ss = 0;
        for (int i = 0; i < m; ++i) ss += data[s*m + i];
        season_avg[s] = ss / m;
    }

    // Seasonal indices for each position in the season
    for (int i = 0; i < m; ++i) {
        double sum = 0;
        int    cnt = 0;
        for (int s = 0; s < seasons; ++s) {
            int idx = s*m + i;
            if (idx < n) {
                sum += data[idx] - season_avg[s];
                cnt++;
            }
        }
        season[i] = (cnt > 0 ? sum / cnt : 0.0);
    }

    if (verbose) {
        cout << fixed << setprecision(6)
             << "\n[Additive] L0=" << level
             << ", T0="           << trend
             << "\nSeasonals:";
        for (int i = 0; i < m; ++i) cout << " " << season[i];
        cout << "\n";
    }

    vec fitted(n);
    vec L{level};  // level over time
    vec T{trend};  // trend over time

    // Recursive Holt–Winters update
    for (int t = 0; t < n; ++t) {
        int si = t % m;  // season index

        // Forecast for time t (one-step ahead)
        double forecast_t = L.back() + T.back() + season[si];
        fitted[t] = forecast_t;

        // Update level, trend, seasonality
        double newL = alpha * (data[t] - season[si]) + (1.0 - alpha) * (L.back() + T.back());
        double newT = beta  * (newL - L.back())      + (1.0 - beta)  * T.back();
        double newS = gamma * (data[t] - newL)       + (1.0 - gamma) * season[si];

        L.push_back(newL);
        T.push_back(newT);
        season[si] = newS;
    }

    // Out-of-sample forecasts
    vec forecast(horizon);
    for (int k = 1; k <= horizon; ++k) {
        forecast[k-1] = L.back() + k * T.back() + season[(n + k - 1) % m];
    }

    // Compute MSE on in-sample data
    double sse = 0;
    for (int i = 0; i < n; ++i) {
        double e = data[i] - fitted[i];
        sse += e*e;
    }
    double mse = sse / n;

    return {fitted, forecast, mse};
}

// Multiplicative Holt–Winters implementation
HWResult holt_winters_multiplicative(
    const vec& data,
    int        m,
    int        horizon,
    double     alpha,
    double     beta,
    double     gamma,
    bool       verbose=false
){
    int n = (int)data.size();
    if (n < 2 * m) return {{}, {}, 1e18}; // not enough data

    // Initial level
    double s1 = accumulate(data.begin(), data.begin() + m, 0.0) / m;
    double level = s1;

    // Initial trend
    double trend_sum = 0;
    for (int i = 0; i < m; ++i)
        trend_sum += (data[m + i] - data[i]);
    double trend = (trend_sum / (double)m) / (double)m;

    // Season averages
    int seasons = n / m;
    vector<double> season_avg(seasons, 0.0);
    for (int s = 0; s < seasons; ++s) {
        double ss = 0;
        for (int i = 0; i < m; ++i) ss += data[s*m + i];
        season_avg[s] = ss / m;
    }

    // Multiplicative seasonality
    vec season(m, 1.0);
    for (int i = 0; i < m; ++i) {
        double sum = 0;
        int    cnt = 0;
        for (int s = 0; s < seasons; ++s) {
            int idx = s*m + i;
            if (idx < n) {
                sum += data[idx] / season_avg[s];
                cnt++;
            }
        }
        season[i] = (cnt > 0 ? sum / cnt : 1.0);
    }

    if (verbose) {
        cout << fixed << setprecision(6)
             << "\n[Multiplicative] L0=" << level
             << ", T0="                  << trend
             << "\nSeasonals:";
        for (int i = 0; i < m; ++i) cout << " " << season[i];
        cout << "\n";
    }

    vec fitted(n);
    vec L{level};
    vec T{trend};

    for (int t = 0; t < n; ++t) {
        int si = t % m;

        // Forecast at t
        double forecast_t = (L.back() + T.back()) * season[si];
        fitted[t] = forecast_t;

        // Update level and trend
        double denom_season = max(1e-12, season[si]);
        double newL = alpha * (data[t] / denom_season)
                    + (1.0 - alpha) * (L.back() + T.back());
        double newT = beta  * (newL - L.back())
                    + (1.0 - beta)  * T.back();

        // Update seasonality
        double denom_level = max(1e-12, newL);
        double newS = gamma * (data[t] / denom_level)
                    + (1.0 - gamma) * season[si];

        L.push_back(newL);
        T.push_back(newT);
        season[si] = newS;
    }

    vec forecast(horizon);
    for (int k = 1; k <= horizon; ++k) {
        forecast[k-1] = (L.back() + k * T.back()) * season[(n + k - 1) % m];
    }

    // MSE
    double sse = 0;
    for (int i = 0; i < n; ++i) {
        double e = data[i] - fitted[i];
        sse += e*e;
    }
    double mse = sse / n;

    return {fitted, forecast, mse};
}

/* ================== GA-BASED ALPHA/BETA/GAMMA CALCULATOR ================== */

// One GA individual (chromosome): a candidate (alpha, beta, gamma)
struct GAParams {
    double alpha;  // smoothing for level
    double beta;   // smoothing for trend
    double gamma;  // smoothing for seasonality
    double mse;    // fitness = MSE on training data
};

// Configuration of GA hyperparameters
struct GAConfig {
    int    pop_size     = 40;   // individuals per generation
    int    generations  = 80;   // number of generations
    double cx_rate      = 0.8;  // crossover probability
    double mut_rate     = 0.2;  // mutation probability
    double mut_std      = 0.05; // standard deviation of Gaussian mutation
    int    tournament_k = 3;    // tournament size for selection
};

// Clamp x to [0,1]
static inline double clamp01(double x){
    if (x < 0.0) return 0.0;
    if (x > 1.0) return 1.0;
    return x;
}

// Uniform random number in [0,1]
static double randu(){
    // RAND_MAX is typically 32767. We divide by RAND_MAX to get values in[0,1].
    return rand() / (double)RAND_MAX;
}

// Evaluate one chromosome: run HW and compute MSE
GAParams eval_chromosome(
    double alpha,
    double beta,
    double gamma,
    const vec& train_vals,
    int season_len,
    bool use_additive
){
    HWResult r;
    if (use_additive) {
        r = holt_winters_additive(
            train_vals, season_len,
            /*horizon*/ season_len,  // horizon does not matter much for training MSE
            alpha, beta, gamma, false
        );
    } else {
        r = holt_winters_multiplicative(
            train_vals, season_len,
            season_len,
            alpha, beta, gamma, false
        );
    }
    GAParams p{alpha, beta, gamma, r.mse};
    return p;
}

// Tournament selection: pick the best among k random individuals
GAParams tournament_select(const vector<GAParams>& pop, int k){
    GAParams best = pop[rand() % pop.size()];
    for (int i = 1; i < k; ++i) {
        const GAParams& cand = pop[rand() % pop.size()];
        if (cand.mse < best.mse) best = cand;
    }
    return best;
}

// Arithmetic crossover between two parents
pair<GAParams,GAParams> crossover(const GAParams& p1, const GAParams& p2){
    double w = randu();          // weight in [0,1]
    GAParams c1 = p1;
    GAParams c2 = p2;

    // child1 = weighted average of parents
    c1.alpha = clamp01(w * p1.alpha + (1-w) * p2.alpha);
    c1.beta  = clamp01(w * p1.beta  + (1-w) * p2.beta);
    c1.gamma = clamp01(w * p1.gamma + (1-w) * p2.gamma);

    // child2 = opposite weighting
    c2.alpha = clamp01(w * p2.alpha + (1-w) * p1.alpha);
    c2.beta  = clamp01(w * p2.beta  + (1-w) * p1.beta);
    c2.gamma = clamp01(w * p2.gamma + (1-w) * p1.gamma);

    return {c1, c2};
}

// Mutate GAParams using Gaussian noise on each gene with probability mut_rate
void mutate(GAParams& p, double mut_rate, double mut_std){
    auto mutate_var = [&](double &x){
        if (randu() < mut_rate) {
            // Box–Muller transform to generate N(0,1)
            double u1 = max(1e-12, randu());
            double u2 = randu();
            double z  = sqrt(-2.0 * log(u1)) * cos(2*M_PI*u2);

            // Add Gaussian noise and clamp
            x = clamp01(x + mut_std * z);
        }
    };

    mutate_var(p.alpha);
    mutate_var(p.beta);
    mutate_var(p.gamma);
}

// Run GA to optimize alpha,beta,gamma
GAParams run_ga_optimize(
    const vec& train_vals,
    int season_len,
    bool use_additive,
    const GAConfig& cfg = GAConfig()
){
    // 1) Initialize population randomly in [0,1]^3
    vector<GAParams> pop(cfg.pop_size);
    for (int i = 0; i < cfg.pop_size; ++i) {
        double a = randu();
        double b = randu();
        double g = randu();
        pop[i] = eval_chromosome(a, b, g, train_vals, season_len, use_additive);
    }

    // Track global best solution ever seen
    GAParams global_best = pop[0];
    for (auto &p : pop)
        if (p.mse < global_best.mse) global_best = p;

    // 2) Evolution loop
    for (int gen = 0; gen < cfg.generations; ++gen) {
        vector<GAParams> new_pop;
        new_pop.reserve(cfg.pop_size);

        // Elitism: keep the global best
        new_pop.push_back(global_best);

        // Fill the rest of the new population
        while ((int)new_pop.size() < cfg.pop_size) {
            // Select parents
            GAParams parent1 = tournament_select(pop, cfg.tournament_k);
            GAParams parent2 = tournament_select(pop, cfg.tournament_k);

            GAParams child1 = parent1;
            GAParams child2 = parent2;

            // With probability cx_rate, perform crossover
            if (randu() < cfg.cx_rate) {
                auto cc = crossover(parent1, parent2);
                child1 = cc.first;
                child2 = cc.second;
            }

            // Mutate children
            mutate(child1, cfg.mut_rate, cfg.mut_std);
            mutate(child2, cfg.mut_rate, cfg.mut_std);

            // Evaluate children (assign MSE)
            child1 = eval_chromosome(child1.alpha, child1.beta, child1.gamma,
                                     train_vals, season_len, use_additive);
            child2 = eval_chromosome(child2.alpha, child2.beta, child2.gamma,
                                     train_vals, season_len, use_additive);

            // Add to new population
            new_pop.push_back(child1);
            if ((int)new_pop.size() < cfg.pop_size)
                new_pop.push_back(child2);

            // Update global best
            if (child1.mse < global_best.mse) global_best = child1;
            if (child2.mse < global_best.mse) global_best = child2;
        }

        // Replace old population
        pop.swap(new_pop);
    }

    return global_best;
}

/* ======================= MAIN ======================= */

int main() {
    // Path to your PSEi CSV file
    string csvfile = "C:\\Users\\rakkh\\OneDrive\\Desktop\\RP\\PSEi Composite Historical Data (1).csv";

    cout << "Loading CSV (date + value)...\n";
    auto dv = load_csv_dates_values(csvfile);
    cout << "Rows read: " << dv.size() << "\n";

    // Split: training (2017–2021) vs future (>=2022)
    vector<DateVal> train_rows, future_rows;
    for (auto &x : dv) {
        if (x.y >= 2017 && x.y <= 2021)
            train_rows.push_back(x);
        else if (x.y >= 2022)
            future_rows.push_back(x);
    }
    cout << "Training rows (2017–2021): " << train_rows.size()
         << ", Future rows (>=2022): "   << future_rows.size() << "\n";

    if (train_rows.size() < 24) {
        cerr << "Not enough training data (need >= 24 months).\n";
        return 1;
    }

    // Extract training values (PSEi levels)
    vec train_vals;
    for (auto &r : train_rows)
        train_vals.push_back(r.val);

    // Monthly data -> season length m = 12
    int season_len = 12;
    int horizon    = 12; // forecast 12 months of 2022

    // Fix seed for reproducibility of GA
    srand(12345);

    // ---------------- GA for additive & multiplicative on original scale ----------------
    cout << "\nOptimizing additive parameters with GA...\n";
    GAParams ga_add = run_ga_optimize(train_vals, season_len, /*use_additive*/ true);
    cout << "GA additive: alpha=" << ga_add.alpha
         << ", beta="             << ga_add.beta
         << ", gamma="            << ga_add.gamma
         << ", MSE="              << ga_add.mse << "\n";

    cout << "\nOptimizing multiplicative parameters with GA...\n";
    GAParams ga_mul = run_ga_optimize(train_vals, season_len, /*use_additive*/ false);
    cout << "GA multiplicative: alpha=" << ga_mul.alpha
         << ", beta="                  << ga_mul.beta
         << ", gamma="                 << ga_mul.gamma
         << ", MSE="                   << ga_mul.mse << "\n";

    // Use GA-optimized parameters
    double add_alpha = ga_add.alpha, add_beta  = ga_add.beta, add_gamma = ga_add.gamma;
    double mul_alpha = ga_mul.alpha, mul_beta  = ga_mul.beta, mul_gamma = ga_mul.gamma;

    // ----- Original-scale Holt–Winters (Additive & Multiplicative) -----
    cout << "\nOriginal-scale Holt-Winters (Additive)...\n";
    HWResult add_orig = holt_winters_additive(
        train_vals, season_len, horizon,
        add_alpha, add_beta, add_gamma,
        false
    );
    cout << "Additive original-scale MSE = " << add_orig.mse << "\n";

    cout << "\nOriginal-scale Holt-Winters (Multiplicative)...\n";
    HWResult mul_orig = holt_winters_multiplicative(
        train_vals, season_len, horizon,
        mul_alpha, mul_beta, mul_gamma,
        false
    );
    cout << "Multiplicative original-scale MSE = " << mul_orig.mse << "\n";

    // Choose model type based on in-sample MSE on original scale
    bool use_additive = (add_orig.mse <= mul_orig.mse);
    cout << "\nChosen model (on original-scale MSE): "
         << (use_additive ? "Additive" : "Multiplicative") << "\n";

    // Forecast for 12 months on original scale
    vec forecast_orig = use_additive ? add_orig.forecast : mul_orig.forecast;

    // ----- Log-scale Holt–Winters (same GA params, but on log data) -----
    cout << "\nLOG-scale Holt-Winters...\n";

    // Transform training values to log scale
    vec train_log(train_vals.size());
    for (size_t i = 0; i < train_vals.size(); ++i)
        train_log[i] = log(train_vals[i]);

    cout << "LOG-scale Additive...\n";
    HWResult add_log = holt_winters_additive(
        train_log, season_len, horizon,
        add_alpha, add_beta, add_gamma,
        false
    );
    cout << "Additive log-scale MSE (on log values) = " << add_log.mse << "\n";

    cout << "LOG-scale Multiplicative...\n";
    HWResult mul_log = holt_winters_multiplicative(
        train_log, season_len, horizon,
        mul_alpha, mul_beta, mul_gamma,
        false
    );
    cout << "Multiplicative log-scale MSE (on log values) = " << mul_log.mse << "\n";

    // Back-transform log forecasts to original scale using exp()
    vec forecast_log_back(horizon);
    if (use_additive) {
        for (int i = 0; i < horizon; ++i) {
            double f_log = (i < (int)add_log.forecast.size()
                            ? add_log.forecast[i]
                            : numeric_limits<double>::quiet_NaN());
            forecast_log_back[i] = std::isnan(f_log) ?
                numeric_limits<double>::quiet_NaN() : exp(f_log);
        }
    } else {
        for (int i = 0; i < horizon; ++i) {
            double f_log = (i < (int)mul_log.forecast.size()
                            ? mul_log.forecast[i]
                            : numeric_limits<double>::quiet_NaN());
            forecast_log_back[i] = std::isnan(f_log) ?
                numeric_limits<double>::quiet_NaN() : exp(f_log);
        }
    }

    // ----- Collect future actual values for 2022 -----
    vec future_vals;
    for (auto &r : future_rows)
        future_vals.push_back(r.val);

    // ----- Write CSV with forecasts and errors -----
    ofstream fo("forecast_2017_2021_vs_2022_add_vs_mul_log_vs_orig_GA.csv");
    fo << "Year,Month,ModelType,Forecast_Orig,Forecast_LogBack,Actual,"
          "Error_Orig,AbsError_Orig,Error_LogBack,AbsError_LogBack\n";

    // Start from last training date (Dec 2021) and move forward month by month
    DateVal last_train = train_rows.back();
    int yy = last_train.y;
    int mm = last_train.m;

    cout << "\nYear Month Model F_orig F_logBack Actual |Err_orig| |Err_logBack|\n";
    string model_name = use_additive ? "Add" : "Mul";

    for (int i = 0; i < horizon; ++i) {
        // Advance month/year
        mm++;
        if (mm > 12) {
            mm = 1;
            yy++;
        }

        // Forecasts and actual
        double f_o = (i < (int)forecast_orig.size()
                      ? forecast_orig[i]
                      : numeric_limits<double>::quiet_NaN());
        double f_l = (i < (int)forecast_log_back.size()
                      ? forecast_log_back[i]
                      : numeric_limits<double>::quiet_NaN());
        double a   = (i < (int)future_vals.size()
                      ? future_vals[i]
                      : numeric_limits<double>::quiet_NaN());

        // Errors: Actual - Forecast
        double err_o = (!std::isnan(f_o) && !std::isnan(a)) ? (a - f_o)
                       : numeric_limits<double>::quiet_NaN();
        double err_l = (!std::isnan(f_l) && !std::isnan(a)) ? (a - f_l)
                       : numeric_limits<double>::quiet_NaN();
        double abso  = std::isnan(err_o) ? numeric_limits<double>::quiet_NaN()
                       : fabs(err_o);
        double absl  = std::isnan(err_l) ? numeric_limits<double>::quiet_NaN()
                       : fabs(err_l);

        // Write row to CSV
        fo << yy << "," << mm << "," << model_name << ","
           << (std::isnan(f_o) ? "" : to_string(f_o)) << ","
           << (std::isnan(f_l) ? "" : to_string(f_l)) << ","
           << (std::isnan(a)   ? "" : to_string(a))   << ","
           << (std::isnan(err_o)? "" : to_string(err_o)) << ","
           << (std::isnan(abso) ? "" : to_string(abso))  << ","
           << (std::isnan(err_l)? "" : to_string(err_l)) << ","
           << (std::isnan(absl) ? "" : to_string(absl))  << "\n";

        // Also print to console for quick inspection
        cout << yy << " " << setw(2) << mm << " " << model_name << " "
             << (std::isnan(f_o) ? 0 : f_o) << " "
             << (std::isnan(f_l) ? 0 : f_l) << " "
             << (std::isnan(a)   ? 0 : a)   << " | "
             << (std::isnan(abso)? 0 : abso) << " "
             << (std::isnan(absl)? 0 : absl) << "\n";
    }

    fo.close();
    cout << "\nCSV written: forecast_2017_2021_vs_2022_add_vs_mul_log_vs_orig_GA.csv\n";
    cout << "Chosen model type based on original-scale MSE: " << model_name << "\n";
    cout << "GA best additive MSE=" << ga_add.mse
         << ", GA best multiplicative MSE=" << ga_mul.mse << "\n";

    return 0;
}
