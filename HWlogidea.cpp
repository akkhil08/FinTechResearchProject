// holt_2017_2021_log_vs_orig_add_vs_mul.cpp
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

using namespace std;
using vec = vector<double>;

struct DateVal { int y,m,d; string datestr; double val; };

/* ======================= CSV HELPERS ======================= */
static inline string trim(const string &s){
    size_t a = 0, b = s.size();
    while (a < b && isspace((unsigned char)s[a])) a++;
    while (b > a && isspace((unsigned char)s[b-1])) b--;
    return s.substr(a, b - a);
}

static vector<string> split_tokens(const string &line, char sep=','){
    vector<string> cols; string cur; bool inq=false;
    for (char c : line) {
        if (c == '"') { inq = !inq; continue; }
        if (c == sep && !inq) { cols.push_back(cur); cur.clear(); }
        else cur.push_back(c);
    }
    cols.push_back(cur);
    return cols;
}

static bool parse_date_mdY(const string &s, int &y, int &m, int &d){
    string t; for (char c : s) if (c != '\"') t.push_back(c);
    t = trim(t);
    if (t.empty()) return false;
    for (char &c : t) if (c == '-') c = '/';
    vector<string> parts; string cur;
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
            } catch(...) {}
        }
        if (yi == -1) yi = 2;
        try {
            y = stoi(trim(parts[yi]));
            vector<int> other;
            for (int i = 0; i < 3; ++i)
                if (i != yi) other.push_back(stoi(trim(parts[i])));
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
    bool first=true;
    while (getline(f,line)) {
        if (!line.empty() && line.back()=='\r') line.pop_back();
        if (first) {
            header = split_tokens(line);
            first=false;
            continue;
        }
        rows.push_back(split_tokens(line));
    }
    f.close();

    auto lower = [&](string s){
        for (char &c : s) c = (char)tolower(c);
        return s;
    };

    int dateIdx=-1, valIdx=-1;
    for (int i=0;i<(int)header.size();++i) {
        string h = lower(trim(header[i]));
        if (h.find("date") != string::npos) dateIdx = i;
        if (h.find("close") != string::npos ||
            h.find("price") != string::npos ||
            h.find("last")  != string::npos) valIdx = i;
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
        int y=0,m=0,d=0;
        if (!parse_date_mdY(sdate,y,m,d)) continue;
        out.push_back({y,m,d,sdate,v});
    }
    sort(out.begin(), out.end(), [](const DateVal &a, const DateVal &b){
        if (a.y != b.y) return a.y < b.y;
        if (a.m != b.m) return a.m < b.m;
        return a.d < b.d;
    });
    return out;
}

/* ======================= HOLT-WINTERS IMPLEMENTATIONS ======================= */
struct HWResult { vec fitted; vec forecast; double mse; };

HWResult holt_winters_additive(const vec& data, int m, int horizon,
                               double alpha, double beta, double gamma,
                               bool verbose=false)
{
    int n = (int)data.size();
    if (n < 2 * m) return {{}, {}, 1e18};

    double s1 = accumulate(data.begin(), data.begin() + m, 0.0) / m;
    double level = s1;

    double trend_sum = 0;
    for (int i = 0; i < m; ++i) trend_sum += (data[m + i] - data[i]);
    double trend = (trend_sum / static_cast<double>(m)) / static_cast<double>(m);

    int seasons = n / m;
    vec season(m, 0.0);
    vector<double> season_avg(seasons, 0.0);
    for (int s = 0; s < seasons; ++s) {
        double ss = 0;
        for (int i = 0; i < m; ++i) ss += data[s*m + i];
        season_avg[s] = ss / m;
    }
    for (int i = 0; i < m; ++i) {
        double sum = 0; int cnt = 0;
        for (int s = 0; s < seasons; ++s) {
            int idx = s*m + i;
            if (idx < n) { sum += data[idx] - season_avg[s]; cnt++; }
        }
        season[i] = (cnt>0 ? sum / cnt : 0.0);
    }

    if (verbose) {
        cout<<fixed<<setprecision(6)
            <<"\n[Additive] L0="<<level<<", T0="<<trend<<"\nSeasonals:";
        for (int i=0;i<m;++i) cout<<" "<<season[i];
        cout<<"\n";
    }

    vec fitted(n); vec L{level}, T{trend};
    for (int t = 0; t < n; ++t) {
        int si = t % m;
        double forecast_t = L.back() + T.back() + season[si];
        fitted[t] = forecast_t;
        double newL = alpha * (data[t] - season[si]) + (1.0 - alpha) * (L.back() + T.back());
        double newT = beta  * (newL - L.back()) + (1.0 - beta)  * T.back();
        double newS = gamma * (data[t] - newL)    + (1.0 - gamma)* season[si];
        L.push_back(newL); T.push_back(newT); season[si] = newS;
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

    double s1 = accumulate(data.begin(), data.begin() + m, 0.0) / m;
    double level = s1;

    double trend_sum = 0;
    for (int i = 0; i < m; ++i) trend_sum += (data[m + i] - data[i]);
    double trend = (trend_sum / static_cast<double>(m)) / static_cast<double>(m);

    int seasons = n / m;
    vector<double> season_avg(seasons, 0.0);
    for (int s = 0; s < seasons; ++s) {
        double ss = 0;
        for (int i = 0; i < m; ++i) ss += data[s*m + i];
        season_avg[s] = ss / m;
    }

    vec season(m, 1.0);
    for (int i = 0; i < m; ++i) {
        double sum = 0; int cnt = 0;
        for (int s = 0; s < seasons; ++s) {
            int idx = s*m + i;
            if (idx < n) { sum += data[idx] / season_avg[s]; cnt++; }
        }
        season[i] = (cnt>0 ? sum / cnt : 1.0);
    }

    if (verbose) {
        cout<<fixed<<setprecision(6)
            <<"\n[Multiplicative] L0="<<level<<", T0="<<trend<<"\nSeasonals:";
        for (int i=0;i<m;++i) cout<<" "<<season[i];
        cout<<"\n";
    }

    vec fitted(n); vec L{level}, T{trend};
    for (int t = 0; t < n; ++t) {
        int si = t % m;
        double forecast_t = (L.back() + T.back()) * season[si];
        fitted[t] = forecast_t;
        double denom_season = max(1e-12, season[si]);
        double newL = alpha * (data[t] / denom_season) + (1.0 - alpha) * (L.back() + T.back());
        double newT = beta  * (newL - L.back())        + (1.0 - beta)  * T.back();
        double denom_level = max(1e-12, newL);
        double newS = gamma * (data[t] / denom_level)  + (1.0 - gamma) * season[si];
        L.push_back(newL); T.push_back(newT); season[si] = newS;
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

/* ======================= MAIN ======================= */

int main() {
    string csvfile = "C:\\Users\\rakkh\\OneDrive\\Desktop\\RP\\PSEi Composite Historical Data (1).csv";

    cout << "Loading CSV (date + value)...\n";
    auto dv = load_csv_dates_values(csvfile);
    cout << "Rows read: " << dv.size() << "\n";

    // Train 2017–2021, forecast 2022
    vector<DateVal> train_rows, future_rows;
    for (auto &x : dv) {
        if (x.y >= 2017 && x.y <= 2021) train_rows.push_back(x);
        else if (x.y >= 2022) future_rows.push_back(x);
    }
    cout << "Training rows (2017–2021): " << train_rows.size()
         << ", Future rows (>=2022): " << future_rows.size() << "\n";

    if (train_rows.size() < 24) {
        cerr << "Not enough training data (need >= 24 months).\n";
        return 1;
    }

    vec train_vals;
    for (auto &r : train_rows) train_vals.push_back(r.val);

    int season_len = 12;
    int horizon    = 12;   // forecast 12 months of 2022

    // GA parameters from Table I (paper)[file:1]
    double add_alpha = 0.82607, add_beta  = 0.0, add_gamma = 1.0;
    double mul_alpha = 0.49288, mul_beta  = 0.0, mul_gamma = 1.0;

    // ----- Original-scale additive & multiplicative -----
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

    // Choose model type based on original-scale MSE
    bool use_additive = (add_orig.mse <= mul_orig.mse);
    cout << "\nChosen model (on original-scale MSE): "
         << (use_additive ? "Additive" : "Multiplicative") << "\n";

    vec forecast_orig = use_additive ? add_orig.forecast : mul_orig.forecast;

    // ----- Log-scale additive & multiplicative -----
    cout << "\nLOG-scale Holt-Winters...\n";
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

    // Use the same model type as decided from original scale
    vec forecast_log_back(horizon);
    if (use_additive) {
        for (int i = 0; i < horizon; ++i) {
            double f_log = (i < (int)add_log.forecast.size()
                            ? add_log.forecast[i]
                            : numeric_limits<double>::quiet_NaN());
            forecast_log_back[i] = std::isnan(f_log) ? numeric_limits<double>::quiet_NaN() : exp(f_log);
        }
    } else {
        for (int i = 0; i < horizon; ++i) {
            double f_log = (i < (int)mul_log.forecast.size()
                            ? mul_log.forecast[i]
                            : numeric_limits<double>::quiet_NaN());
            forecast_log_back[i] = std::isnan(f_log) ? numeric_limits<double>::quiet_NaN() : exp(f_log);
        }
    }

    // ----- Future actuals (2022) -----
    vec future_vals;
    for (auto &r : future_rows) future_vals.push_back(r.val);

    // ----- CSV with errors -----
    ofstream fo("forecast_2017_2021_vs_2022_add_vs_mul_log_vs_orig.csv");
    fo << "Year,Month,ModelType,Forecast_Orig,Forecast_LogBack,Actual,"
          "Error_Orig,AbsError_Orig,Error_LogBack,AbsError_LogBack\n";

    DateVal last_train = train_rows.back();
    int yy = last_train.y;
    int mm = last_train.m;

    cout << "\nYear Month Model F_orig F_logBack Actual |Err_orig| |Err_logBack|\n";
    string model_name = use_additive ? "Add" : "Mul";

    for (int i = 0; i < horizon; ++i) {
        mm++; if (mm > 12) { mm = 1; yy++; }

        double f_o = (i < (int)forecast_orig.size()
                      ? forecast_orig[i]
                      : numeric_limits<double>::quiet_NaN());
        double f_l = (i < (int)forecast_log_back.size()
                      ? forecast_log_back[i]
                      : numeric_limits<double>::quiet_NaN());
        double a   = (i < (int)future_vals.size()
                      ? future_vals[i]
                      : numeric_limits<double>::quiet_NaN());

        double err_o = (!std::isnan(f_o) && !std::isnan(a)) ? (a - f_o) : numeric_limits<double>::quiet_NaN();
        double err_l = (!std::isnan(f_l) && !std::isnan(a)) ? (a - f_l) : numeric_limits<double>::quiet_NaN();
        double abso  = std::isnan(err_o) ? numeric_limits<double>::quiet_NaN() : fabs(err_o);
        double absl  = std::isnan(err_l) ? numeric_limits<double>::quiet_NaN() : fabs(err_l);

        fo << yy << "," << mm << "," << model_name << ","
           << (std::isnan(f_o) ? "" : to_string(f_o)) << ","
           << (std::isnan(f_l) ? "" : to_string(f_l)) << ","
           << (std::isnan(a)   ? "" : to_string(a))   << ","
           << (std::isnan(err_o)? "" : to_string(err_o)) << ","
           << (std::isnan(abso) ? "" : to_string(abso))  << ","
           << (std::isnan(err_l)? "" : to_string(err_l)) << ","
           << (std::isnan(absl) ? "" : to_string(absl))  << "\n";

        cout << yy << " " << setw(2) << mm << " " << model_name << " "
             << (std::isnan(f_o) ? 0 : f_o) << " "
             << (std::isnan(f_l) ? 0 : f_l) << " "
             << (std::isnan(a)   ? 0 : a)   << " | "
             << (std::isnan(abso)? 0 : abso) << " "
             << (std::isnan(absl)? 0 : absl) << "\n";
    }

    fo.close();
    cout << "\nCSV written: forecast_2017_2021_vs_2022_add_vs_mul_log_vs_orig.csv\n";
    cout << "Chosen model type based on original-scale MSE: " << model_name << "\n";
    return 0;
}
