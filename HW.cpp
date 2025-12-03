#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <cctype>

using namespace std;
using vec = vector<double>;

/* ================= CSV LOADER ================= */
/*
  load_csv_values:
  - Reads a CSV file line by line.
  - Handles quoted fields (very simple handling).
  - Detects header column with the word "close" (case-insensitive).
  - Returns numeric values from that column as a vector<double>.
*/
vector<double> load_csv_values(const string &filename) {
    ifstream f(filename);
    if (!f.is_open()) { cerr<<"Cannot open file: "<<filename<<"\n"; exit(1); }

    string line;
    vector<string> header;
    vector<vector<string>> rows;
    bool first = true;

    // read lines
    while (getline(f, line)) {
        if (!line.empty() && line.back()=='\r') line.pop_back(); // handle CRLF on Windows (Carriage Return and Line feed )
        vector<string> cols; string cur; bool inq=false;

        // simple CSV tokenization supporting quotes:
        for (char c: line) {
            if (c=='"') { inq=!inq; continue; }          // toggle quote state, drop quotes
            if (c==',' && !inq) { cols.push_back(cur); cur.clear(); } // separator outside quotes
            else cur.push_back(c);
        }
        cols.push_back(cur);

        if (first) { header = cols; first = false; continue; } // first line is header
        rows.push_back(cols);
    }

    int ncols = header.size();

    // lowercase helper (capture by value; returns lowered string)
    auto lower = [&](string s){ for(char &c:s) c=tolower(c); return s; };

    // try to find a header column containing "close" (case-insensitive)
    int chosen=-1;
    for (int i=0;i<ncols;i++)
        if (lower(header[i]).find("close") != string::npos) chosen=i;

    // fallback: use column index 1 if nothing found (risky but keeps behavior)
    if (chosen == -1) chosen = 1;

    vector<double> vals;
    for (auto &r : rows) {
        if (chosen >= r.size()) continue;            // skip rows shorter than header
        string s = r[chosen];
        // remove commas inside number (thousands separators)
        s.erase(remove(s.begin(), s.end(), ','), s.end());
        try { vals.push_back(stod(s)); } catch (...) {} // skip non-numeric entries
    }

    return vals;
}

/* ================= Holt-Winters ADDITIVE ================= */
/*
  holt_winters_additive:
  - data: time series values (assumes regular spacing) (Monthly closing prices vector in this case)
  - m: season length (e.g. 12 for monthly data with yearly seasonality)
  - horizon: how many future steps to forecast
  - alpha/beta/gamma: smoothing coefficients (0..1)
  Returns: fitted values (the in sample fitted value one per input time point ), forecast vector, mse 
*/
struct HWResult { vec fitted; vec forecast; double mse; };

HWResult holt_winters_additive(const vec& data,int m,int horizon,
                               double alpha,double beta,double gamma)
{
    int n=data.size(); //time series value length
    // need at least 2 full seasons to initialize s1/s2
    if (n < 2*m) return {{},{},1e18}; //code needs at least 2 full seasons to compute initial level/trend

    // initial level and trend using first two seasons' means
    double s1 = accumulate(data.begin(), data.begin()+m, 0.0)/m; //sum of first season data 0-m-1 (assume m=12 we get data[0] to data[11])
    double s2 = accumulate(data.begin()+m, data.begin()+2*m, 0.0)/m;//sum of second season data m-2m-1 (assume m=12 we get data[12] to data[23])
    double level = s1; // initial level L0
    double trend = (s2 - s1)/m; // initial average per-period change T0  ## Give an explanation for divide by m, or remove it with explanation

    int seasons = n / m; // number of complete seasons in data 
    // if (seasons < 1) seasons = 1; //(assume n=11 and m=12 we have less than 1 season, set seasons to 1 to avoid issues)

    // compute average for each season (season_avg[s] = mean of season s)
    vec season(m,0.0);
    vector<double> season_avg(seasons,0.0);
    for (int s=0;s<seasons;s++){
        double ss=0;
        for (int i=0;i<m;i++) ss+=data[s*m+i];
        season_avg[s] = ss/m; //(here suppose m=12 we get average of data[0] to data[11] for s=0, data[12] to data[23] for s=1, etc)
    }

    // initialize seasonal factors: average of (data - season_avg[s]) for each seasonal position basically measures deviation from season average example season avg is 120 for month 1, 130 for month 2, etc how much each month deviates from that yearly average
    for (int i=0;i<m;i++){
        double sum=0; int cnt=0;
        for (int s=0;s<seasons;s++){
            int idx = s*m + i;
            if (idx<n){ sum+=data[idx]-season_avg[s]; cnt++; }
        }
        season[i] = (cnt>0? sum/cnt : 0.0); // additive seasonal offsets
    }

    vec fitted(n);
    vec L{level}, T{trend}; // store history for L and T (push_back used during loop)

    // main smoothing loop computing the present fitted values
    for (int t=0;t<n;t++){
        int si = t % m; //seasonal index
        // forecast for time t using most recent level & trend (L.back(), T.back()) and seasonal factor
        double f = L.back() + T.back() + season[si];
        fitted[t] = f;

        // update equations and update of values (use observed data[t])
        double newL = alpha*(data[t]-season[si]) + (1-alpha)*(L.back()+T.back());
        double newT = beta*(newL-L.back()) + (1-beta)*T.back();
        // update seasonal factor: gamma*(obs - level) + (1-gamma)*oldSeason
        double newS = gamma*(data[t]-newL) + (1-gamma)*season[si];

        L.push_back(newL);
        T.push_back(newT);
        season[si] = newS; // update in-place for this seasonal position
    }

    // produce h-step forecast: use last level/trend and seasonal cycle like how much ahead you wish to forecast
    vec forecast(horizon);
    for (int k=1;k<=horizon;k++)
        forecast[k-1] = L.back() + k*T.back() + season[(n+k-1)%m];

    // compute MSE between observed and fitted
    double sse=0;
    for (int i=0;i<n;i++){
        double e=data[i]-fitted[i];
        sse+=e*e;
    }

    return {fitted, forecast, sse/n};
}
/*
Quick worked example (very small)

Suppose m=4 and data = {11, 9, 12, 10, 13, 11, 14, 12} (two seasons), alpha=0.5, beta=0.5, gamma=0.5.

s1 = mean(11,9,12,10) = 10.5 → level = 10.5

s2 = mean(13,11,14,12) = 12.5 → trend = (12.5-10.5)/4 = 0.5

Initial seasonal offsets would be computed for each position i=0..3 as average of (value - that season's mean) across seasons. 
Then the loop computes fitted values and updates L,T,S for each t, producing fitted outputs and forecasts for future k steps. 
*/
/* ================= Holt-Winters MULTIPLICATIVE ================= */
/*
  Same structure as additive, but seasonal factors multiply the level+trend.
  Be careful: multiplicative model requires strictly positive data (or special handling).
*/
HWResult holt_winters_multiplicative(const vec& data,int m,int horizon,
                                     double alpha,double beta,double gamma)
{
    int n=data.size();
    if (n < 2*m) return {{},{},1e18};

    int seasons = n/m;
    if (seasons < 1) seasons = 1;

    // initial level and trend (same simple approach)
    double s1 = accumulate(data.begin(), data.begin()+m, 0.0)/m;
    double s2 = accumulate(data.begin()+m, data.begin()+2*m, 0.0)/m;
    double level = s1;
    double trend = (s2 - s1)/m;

    // compute season averages
    vector<double> season_avg(seasons,0.0);
    for (int s=0;s<seasons;s++){
        double ss=0; for (int i=0;i<m;i++) ss+=data[s*m+i];
        season_avg[s] = ss/m;
    }

    // initialize multiplicative seasonal factors (ratios)
// average of (data / season_avg[s]) for each seasonal position
// basically measures the *ratio* of each month's value to that season's average.
// Example: if season average is 120 for month 1, but actual month 1 values
// across years are around 1.10 * 120 = 132, then seasonal factor ≈ 1.10.
// This means month 1 is typically 10% above its season's mean.
// If seasonal factor ≈ 0.90, that month is typically 10% below the mean.
    vec season(m,1.0);
    for (int i=0;i<m;i++){
        double sum=0; int cnt=0;
        for (int s=0;s<seasons;s++){
            int idx=s*m+i;
            if (idx<n){ sum += data[idx]/season_avg[s]; cnt++; }
        }
        season[i]= (cnt>0? sum/cnt : 1.0);
    }

    vec fitted(n);
    vec L{level}, T{trend};

    // smoothing loop
    for (int t=0;t<n;t++){
        int si = t % m;

        double f = (L.back() + T.back()) * season[si]; // multiplicative forecast
        fitted[t] = f;

        // important: avoid division by zero: newL denominator is season[si]
        // (code does not include epsilon — consider adding a small value guard)
        double newL = alpha*(data[t]/season[si]) + (1-alpha)*(L.back()+T.back());
        double newT = beta*(newL-L.back()) + (1-beta)*T.back();
        double newS = gamma*(data[t]/newL) + (1-gamma)*season[si];

        L.push_back(newL);
        T.push_back(newT);
        season[si] = newS;
    }

    vec forecast(horizon);
    for (int k=1;k<=horizon;k++)
        forecast[k-1] = (L.back()+k*T.back()) * season[(n+k-1)%m];

    double sse=0;
    for (int i=0;i<n;i++){
        double e=data[i]-fitted[i];
        sse+=e*e;
    }

    return {fitted, forecast, sse/n};
}

/* =================== MAIN ====================== */
int main() {
    // path to CSV on the original machine — change to your file path
    string csvfile = "C:\\Users\\rakkh\\OneDrive\\Desktop\\RP\\PSEi Composite Historical Data (1).csv";

    cout<<"Loading CSV...\n";
    vec data = load_csv_values(csvfile);
    cout<<"Loaded "<<data.size()<<" rows.\n\n";

    // smoothing parameters (tune these)
    double alpha = 0.82607;
    double beta  = 0.0;
    double gamma = 1.0;

    int season_len = 12; // e.g. month-over-year seasonality
    int horizon    = 24; // forecast 24 steps ahead

    cout<<"Running ADDITIVE Holt-Winters...\n";
    HWResult add = holt_winters_additive(data, season_len, horizon, alpha, beta, gamma);
    cout<<"Additive MSE = "<<add.mse<<"\n\n";

    cout<<"Running MULTIPLICATIVE Holt-Winters...\n";
    HWResult mul = holt_winters_multiplicative(data, season_len, horizon, alpha, beta, gamma);
    cout<<"Multiplicative MSE = "<<mul.mse<<"\n\n";

    // ---------- Save forecasts ----------
    ofstream fa("forecast_additive.csv");
    fa<<"Step,Forecast\n";
    for(int i=0;i<horizon;i++) fa<<(i+1)<<","<<add.forecast[i]<<"\n";
    fa.close();

    ofstream fm("forecast_multiplicative.csv");
    fm<<"Step,Forecast\n";
    for(int i=0;i<horizon;i++) fm<<(i+1)<<","<<mul.forecast[i]<<"\n";
    fm.close();

    cout<<"Forecasts saved:\n";
    cout<<" - forecast_additive.csv\n";
    cout<<" - forecast_multiplicative.csv\n\n";

    // ---------- Which model is better ----------
    if (add.mse < mul.mse)
        cout<<"➡ BEST MODEL: ADDITIVE (lower MSE)\n";
    else
        cout<<"➡ BEST MODEL: MULTIPLICATIVE (lower MSE)\n";

    return 0;
}
