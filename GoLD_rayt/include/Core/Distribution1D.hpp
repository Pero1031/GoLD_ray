// Distribution1D.hpp
#pragma once

#include <vector>
#include <numeric>
#include <algorithm>

namespace rayt {

    struct Distribution1D {
        std::vector<float> func; // 元の関数値（輝度など）
        std::vector<float> cdf;  // 累積分布関数 (Cumulative Distribution Function)
        float funcInt;           // 関数の積分値（全体のエネルギー合計）

        // データを渡してCDFを構築する
        Distribution1D(const float* f, int n) : func(f, f + n) {
            cdf.resize(n + 1);
            cdf[0] = 0;

            // 累積和を計算 (CDFを作る)
            for (int i = 0; i < n; ++i) {
                cdf[i + 1] = cdf[i] + func[i] / n; // 1/n は区間幅
            }

            funcInt = cdf[n]; // 合計値

            // 正規化 (CDFの最後が1.0になるように割る)
            if (funcInt == 0) {
                for (int i = 1; i < n + 1; ++i) cdf[i] = float(i) / float(n);
            }
            else {
                for (int i = 1; i < n + 1; ++i) cdf[i] /= funcInt;
            }
        }

        // 0.0〜1.0の乱数(u)を受け取り、確率に基づいたインデックスを返す
        // pdf: 選択された場所の確率密度 (Probability Density Function)
        int sampleContinuous(float u, float& pdf, int& off) const {
            // CDFの中から u に対応する場所を二分探索で探す
            auto ptr = std::upper_bound(cdf.begin(), cdf.end(), u);
            int offset = std::max(0, (int)(ptr - cdf.begin() - 1));

            off = offset; // 選ばれたインデックス

            // PDF（確率密度）を計算して返す
            // PDF = (値 / 合計) * データ数
            // ※ 確率が高い場所ほど大きな値になる
            pdf = (func[offset] / funcInt) * func.size();

            // 厳密にはオフセット内の小数を計算してリマップしますが、
            // まずは「どの画素か」が特定できればOKです
            return offset;
        }
    };

}