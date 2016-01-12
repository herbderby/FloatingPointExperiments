#include <iostream>
#include <random>
#include <tuple>
#include <limits>
#include <string>
#include <functional>

using namespace std;

class StatsAccumulator {
public:
    void Add(double v) {
        fCount += 1.0;
        double delta = v - fMean;
        fMean += delta / fCount;
        fM2 += delta * (v - fMean);
        fMin = min(fMin, v);
        fMax = max(fMax, v);
    }

    void Out(string label) {
        cout << label << ": ";
        if (fCount >= 0.0) {
            cout << " count = " << fCount;
        }
        if (fCount >= 1.0) {
            cout << " mean = " << fMean;
            cout << " min = " << fMin;
            cout << " max = " << fMax;
        }
        if (fCount >= 2.0) {
            cout << " stddef = " << sqrt(fM2 / (fCount - 1.0));
        }
        cout << endl;
    }

private:
    double fCount{0.0};
    double fMean{0.0};
    double fM2{0.0};
    double fMin{numeric_limits<double>::max()};
    double fMax{numeric_limits<double>::lowest()};
};

template<typename... Args>
class TestN {
public:
    TestN(tuple<string, function<double(Args&&...)>> ideal,
          initializer_list<tuple<string, function<float(Args&&...)>>> candidates)
        : fIdeal{ideal}, fStats{candidates.size()} {
            for (auto& c : candidates) {
                fCandidates.push_back(c);
            }
        }

    void Try(Args&... args) {
        double ideal = get<1>(fIdeal)(forward<Args>(args)...);
        double errors[fCandidates.size()];
        for (int i = 0; i < fCandidates.size(); i++) {
            double v = get<1>(fCandidates[i])(forward<Args>(args)...);
            double e = abs((ideal - v) / ideal);
            fStats[i].Add(e);
        }
    }

    void Report() {
        for (int i = 0; i < fCandidates.size(); i++) {
            string name = get<0>(fCandidates[i]);
            fStats[i].Out(name);
        }
    }

private:
    tuple<string, function<double(Args... args)>> fIdeal;
    vector<tuple<string, function<float(Args...)>>> fCandidates;
    vector<StatsAccumulator> fStats;
};

int main() {
    std::random_device rd;
    mt19937_64 mt(rd());
    uniform_real_distribution<float> randomT(0.0f, 1.0f);
    uniform_real_distribution<float> randomP(-33000000.0f, 33000000.0f);
    uniform_real_distribution<float> randomSmallP(-33000.0f, 33000.0f);

    TestN<float, float, float> tt
        {make_tuple("Ideal: t * p0 + (1 - t) * p1",
                    [](float t, float p0, float p1){
                        return (double)t * p0 + (1.0 - t) * p1;
                    }),
         {make_tuple("t * p0 + (1 - t) * p1",
                     [](float t, float p0, float p1) {
                         return t * p0 + (1.0 - t) * p1;
                     }),
          make_tuple("p1 + t * (p0 - p1)",
                     [](float t, float p0, float p1){
                         return p1 + t * (p0 - p1);
                     }),
          make_tuple("t * p0 - t * p1 + p1",
                     [](float t, float p0, float p1){
                         return t * p0 - t * p1 + p1;
                     })
         }
        };

    for (int i = 0; i < 500; i++) {
        float t = randomT(mt);
        float p0 = randomSmallP(mt);
        float p1 = randomSmallP(mt);
        tt.Try(t, p0, p1);
    }
    tt.Report();
    int best[3] = {0, 0, 0};
    int same[3] = {0, 0, 0};
    StatsAccumulator accum[3];
    for (int i = 0; i < 500; i++) {
        float t = randomT(mt);
        float p0 = randomP(mt);
        float p1 = randomP(mt);

        float v1 = t * p0 + (1.0f - t) * p1;
        float v2 = p1 + t * (p0 - p1);
        float v3 = t * p0 - t * p1 + p1;
        double dv = (double) t * p0 + (1.0 - t) * p1;
        double e1 = std::abs((v1 - dv) / dv);
        double e2 = std::abs((v2 - dv) / dv);
        double e3 = std::abs((v3 - dv) / dv);
        vector <tuple<int, double>> s{make_tuple(0, e1), make_tuple(1, e2),
        make_tuple(2, e3)};
        sort(s.begin(), s.end(), [](auto a, auto b) {return get<1>(a) < get<1>(b);});
        best[get<0>(s[0])] += 1;
        accum[0].Add(e1);
        accum[1].Add(e2);
        accum[2].Add(e3);
        if (e1 != e2) {
            cout << "winner: " << get<0>(s[0]);
            cout << " v1:e1 - " << v1 << " : " << e1;
            cout << " v2:e2 - " << v2 << " : " << e2;
            cout << " v3:e3 - " << v3 << " : " << e3;
            cout << endl;
        }
        if (e1 == e2) { same[0] += 1; }
        if (e2 == e3) { same[1] += 1; }
        if (e1 == e3) { same[2] += 1;}
    }

    for (int i = 0; i < 3; i++) {
        cout << i << " best: " << best[i] << endl;
        accum[i].Out(to_string(i));
    }
    for (int i = 0; i < 3; i++) {
        cout << i << " same: " << best[i] << endl;
    }
}
