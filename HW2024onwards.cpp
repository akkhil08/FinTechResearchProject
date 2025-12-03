#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <cctype>
#include <sstream>

using namespace std;
using vec = vector<double>;

struct DateVal { int y,m,d; string datestr; double val; };

/* ================= CSV LOADER (returns ordered date+value pairs) ================= */

// Trim helper
static inline string trim(const string &s){ size_t a=0,b=s.size(); while(a<b && isspace((unsigned char)s[a])) a++; while(b>a && isspace((unsigned char)s[b-1])) b--; return s.substr(a,b-a); }

// split helper
static vector<string> split_tokens(const string &line, char sep=','){
    vector<string> cols; string cur; bool inq=false;
    for (char c: line) {
        if (c=='"') { inq=!inq; continue; }
        if (c==sep && !inq) { cols.push_back(cur); cur.clear(); } else cur.push_back(c);
    }
    cols.push_back(cur);
    return cols;
}

// parse date string expecting something like MM/DD/YYYY or DD/MM/YYYY or YYYY-MM-DD
// We'll try to extract year as last numeric 4-digit token and month/day from the other tokens.
static bool parse_date_mdY(const string &s, int &y, int &m, int &d){
    string x = s;
    // remove quotes and trim
    string t; for(char c:x) if(c!='\"') t.push_back(c);
    t = trim(t);
    if (t.empty()) return false;
    // replace '-' with '/' to unify
    for (char &c : t) if (c=='-') c='/';
    vector<string> parts;
    string cur;
    for (char c: t){ if (c=='/') { parts.push_back(cur); cur.clear(); } else cur.push_back(c);} parts.push_back(cur);
    // remove empty items
    vector<string> good;
    for (auto &p: parts){ string q = trim(p); if(!q.empty()) good.push_back(q); }
    if (good.size() < 3) return false;
    // assume year is the element that has length 4 (or the last element)
    int yi=-1; for (int i=0;i<(int)good.size();++i) if ((int)good[i].size()==4){ yi=i; break; }
    if (yi==-1) yi = (int)good.size()-1;
    string sy = good[yi];
    try { y = stoi(sy); } catch(...) { return false; }
    // take other two as month/day in their order
    vector<string> others;
    for (int i=0;i<(int)good.size();++i) if (i!=yi) others.push_back(good[i]);
    // If format was YYYY/MM/DD, others will be {MM,DD} when yi==0, so handle
    if (yi == 0 && others.size()>=2) { // original was Y/M/D
        try{ m = stoi(others[0]); d = stoi(others[1]); return true; } catch(...) { return false; }
    }
    // otherwise assume order is M/D
    try{ m = stoi(others[0]); d = stoi(others[1]); } catch(...) { return false; }
    return true;
}

vector<DateVal> load_csv_dates_values(const string &filename) {
    ifstream f(filename);
    if (!f.is_open()) { cerr<<"Cannot open file: "<<filename<<"\n"; exit(1); }

    string line;
    vector<string> header;
    vector<vector<string>> rows;
    bool first = true;

    while (getline(f, line)) {
        if (!line.empty() && line.back()=='\r') line.pop_back();
        if (first){ header = split_tokens(line); first=false; continue; }
        auto cols = split_tokens(line);
        rows.push_back(cols);
    }
    f.close();

    // find date column and value/close column (case-insensitive)
    auto lower = [&](string s){ for(char &c:s) c=tolower(c); return s; };
    int dateIdx=-1, valIdx=-1;
    for (int i=0;i<(int)header.size();++i){ string h = lower(trim(header[i])); if (h.find("date")!=string::npos) dateIdx=i; if (h.find("close")!=string::npos || h.find("price")!=string::npos || h.find("last")!=string::npos) valIdx=i; }
    if (dateIdx==-1) dateIdx=0; // fallback first column
    if (valIdx==-1) valIdx=1;   // fallback second column

    vector<DateVal> out;
    for (auto &r : rows) {
        if (dateIdx >= r.size() || valIdx >= r.size()) continue;
        string sdate = trim(r[dateIdx]);
        string sval  = r[valIdx];
        // remove commas in number
        sval.erase(remove(sval.begin(), sval.end(), ','), sval.end());
        if (sval.empty()) continue;
        double v;
        try { v = stod(sval); } catch(...) { continue; }
        int y=0,m=0,d=0;
        if (!parse_date_mdY(sdate,y,m,d)) continue;
        DateVal dv{y,m,d,sdate,v};
        out.push_back(dv);
    }

    // sort by date ascending (oldest first)
    sort(out.begin(), out.end(), [](const DateVal &a,const DateVal &b){ if (a.y!=b.y) return a.y<b.y; if (a.m!=b.m) return a.m<b.m; return a.d<b.d; });

    return out;
}

/* ================= Holt-Winters (unchanged) ================= */
struct HWResult { vec fitted; vec forecast; double mse; };

HWResult holt_winters_additive(const vec& data,int m,int horizon,
                               double alpha,double beta,double gamma)
{
    int n=data.size();
    if (n < 2*m) return {{},{},1e18};

    double s1 = accumulate(data.begin(), data.begin()+m, 0.0)/m;
    double s2 = accumulate(data.begin()+m, data.begin()+2*m, 0.0)/m;
    double level = s1;
    double trend = (s2 - s1)/m;

    int seasons = n / m;
    if (seasons < 1) seasons = 1;

    vec season(m,0.0);
    vector<double> season_avg(seasons,0.0);
    for (int s=0;s<seasons;s++){
        double ss=0;
        for (int i=0;i<m;i++) ss+=data[s*m+i];
        season_avg[s] = ss/m;
    }

    for (int i=0;i<m;i++){
        double sum=0; int cnt=0;
        for (int s=0;s<seasons;s++){
            int idx = s*m + i;
            if (idx<n){ sum+=data[idx]-season_avg[s]; cnt++; }
        }
        season[i] = (cnt>0? sum/cnt : 0.0);
    }

    vec fitted(n);
    vec L{level}, T{trend};

    for (int t=0;t<n;t++){
        int si = t % m;
        double f = L.back() + T.back() + season[si];
        fitted[t] = f;

        double newL = alpha*(data[t]-season[si]) + (1-alpha)*(L.back()+T.back());
        double newT = beta*(newL-L.back()) + (1-beta)*T.back();
        double newS = gamma*(data[t]-newL) + (1-gamma)*season[si];

        L.push_back(newL);
        T.push_back(newT);
        season[si] = newS;
    }

    vec forecast(horizon);
    for (int k=1;k<=horizon;k++)
        forecast[k-1] = L.back() + k*T.back() + season[(n+k-1)%m];

    double sse=0;
    for (int i=0;i<n;i++){ double e=data[i]-fitted[i]; sse+=e*e; }

    return {fitted, forecast, sse/n};
}

HWResult holt_winters_multiplicative(const vec& data,int m,int horizon,
                                     double alpha,double beta,double gamma)
{
    int n=data.size();
    if (n < 2*m) return {{},{},1e18};

    int seasons = n/m;
    if (seasons < 1) seasons = 1;

    double s1 = accumulate(data.begin(), data.begin()+m, 0.0)/m;
    double s2 = accumulate(data.begin()+m, data.begin()+2*m, 0.0)/m;
    double level = s1;
    double trend = (s2 - s1)/m;

    vector<double> season_avg(seasons,0.0);
    for (int s=0;s<seasons;s++){
        double ss=0; for (int i=0;i<m;i++) ss+=data[s*m+i];
        season_avg[s] = ss/m;
    }

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

    for (int t=0;t<n;t++){
        int si = t % m;
        double f = (L.back() + T.back()) * season[si];
        fitted[t] = f;

        double newL = alpha*(data[t]/max(1e-12, season[si])) + (1-alpha)*(L.back()+T.back());
        double newT = beta*(newL-L.back()) + (1-beta)*T.back();
        double newS = gamma*(data[t]/max(1e-12, newL)) + (1-gamma)*season[si];

        L.push_back(newL);
        T.push_back(newT);
        season[si] = newS;
    }

    vec forecast(horizon);
    for (int k=1;k<=horizon;k++)
        forecast[k-1] = (L.back()+k*T.back()) * season[(n+k-1)%m];

    double sse=0;
    for (int i=0;i<n;i++){ double e=data[i]-fitted[i]; sse+=e*e; }

    return {fitted, forecast, sse/n};
}

/* =================== MAIN: train on data up to Dec 2023, forecast for 2024+ ====================== */
int main() {
    string csvfile = "C:\\Users\\rakkh\\OneDrive\\Desktop\\RP\\PSEi Composite Historical Data (1).csv"; // change path as needed

    cout<<"Loading CSV (date+close)...\n";
    auto dv = load_csv_dates_values(csvfile);
    cout<<"Total rows read: "<<dv.size()<<"\n";

    // split into training (<= 2023-12-31) and test (>= 2024-01-01)
    vector<DateVal> train_rows, test_rows;
    for (auto &x: dv){ if (x.y <= 2023) train_rows.push_back(x); else test_rows.push_back(x); }

    cout<<"Training rows (<= Dec 2023): "<<train_rows.size()<<"\n";
    cout<<"Test rows (>= Jan 2024): "<<test_rows.size()<<"\n\n";

    if (train_rows.size() < 24){ cerr<<"Not enough training data (need at least 2 full seasons).\n"; return 1; }

    // build training vector (values only)
    vec train_vals; train_vals.reserve(train_rows.size());
    for (auto &r: train_rows) train_vals.push_back(r.val);

    // horizon equals number of months in test set (forecast for Jan 2024 till date)
    int horizon = (int)test_rows.size();
    if (horizon == 0) horizon = 12; // default fallback if no test data found

    // smoothing parameters
    double alpha = 0.82607;
    double beta  = 0.0;
    double gamma = 1.0;
    int season_len = 12;

    cout<<"Running ADDITIVE Holt-Winters on training data...\n";
    HWResult add = holt_winters_additive(train_vals, season_len, horizon, alpha, beta, gamma);
    cout<<"Additive in-sample MSE = "<<add.mse<<"\n";

    cout<<"Running MULTIPLICATIVE Holt-Winters on training data...\n";
    HWResult mul = holt_winters_multiplicative(train_vals, season_len, horizon, alpha, beta, gamma);
    cout<<"Multiplicative in-sample MSE = "<<mul.mse<<"\n\n";

    // choose best by in-sample MSE
    bool pick_add = (add.mse < mul.mse);
    cout<<(pick_add?"➡ BEST MODEL: ADDITIVE\n":"➡ BEST MODEL: MULTIPLICATIVE\n");

    vec &chosen_forecast = (pick_add? add.forecast : mul.forecast);

    // Save forecast CSV with Date, Forecast, Actual (if available)
    ofstream fo("forecast_2024_onwards.csv");
    fo<<"Date,Forecast,Actual\n";
    // if we have test_rows, use their dates; otherwise build month sequence after last training date
    if (!test_rows.empty()){
        for (int i=0;i<horizon;i++){
            string d = test_rows[i].datestr;
            double f = chosen_forecast[i];
            double a = test_rows[i].val;
            fo<<"\""<<d<<"\","<<f<<","<<a<<"\n";
        }
    } else {
        // build future dates from last training date
        DateVal last = train_rows.back();
        int yy = last.y, mm = last.m;
        for (int i=0;i<horizon;i++){
            mm++; if (mm>12){ mm=1; yy++; }
            ostringstream ds; ds<< (mm<10?"0":"")<<mm <<"/01/"<<yy; // day set to 01
            double f = chosen_forecast[i];
            fo<<"\""<<ds.str()<<"\","<<f<<","<<""<<"\n";
        }
    }
    fo.close();
    cout<<"Saved forecasts to forecast_2024_onwards.csv\n";

    // If we have test actuals, compute forecast errors
    if (!test_rows.empty()){
        double sse=0, sae=0; int n=test_rows.size();
        double smape_num=0, smape_den=0;
        for (int i=0;i<n;i++){
            double f = chosen_forecast[i];
            double a = test_rows[i].val;
            double e = a - f;
            sse += e*e;
            sae += fabs(e);
            if (fabs(a)+fabs(f) > 1e-12) smape_num += fabs(a-f); smape_den += (fabs(a)+fabs(f))/2.0;
        }
        double mse = sse / n;
        double mae = sae / n;
        double smape = (smape_den>0 ? (smape_num / smape_den) * 100.0 : 0.0);
        cout<<"Evaluation on 2024+ actuals ("<<n<<" points):\n";
        cout<<" - MSE = "<<mse<<"\n";
        cout<<" - MAE = "<<mae<<"\n";
        cout<<" - sMAPE (%) = "<<smape<<"\n";
    } else {
        cout<<"No actual 2024+ rows found in CSV to evaluate forecasts.\n";
    }

    return 0;
}
