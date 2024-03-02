// NAIR: Network Analysis of Immune Repertoire
// Copyright (C) 2023 Li Zhang
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
//   (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.

#define ARMA_64BIT_WORD 1
#include <RcppArmadillo.h>
#include <set>
#include <vector>
#include <strings.h>
#include <algorithm>
#include <unordered_map>
#include "dropDegreeZero.h"
#include "hamDistBounded.h"
#include "levDistBounded.h"
using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]


inline int levDistBoundedMock(const std::string &a, const std::string &b, const int& k) {
    return levDistBounded(a, b, k);
}


void processWords(
    std::vector<std::string>& sortedStrings, 
    std::unordered_map<std::string, std::vector<int>>& str2idx, 
    arma::sp_umat& out, 
    std::function<int(const std::string &, const std::string &, const int&)> distBounded,
    bool reverseStrings
) {
    for (const auto &target_word : sortedStrings) {
        int half_len;
        std::string half_word, half_word_next;
        if (target_word.size() > 1) {
            if (reverseStrings) {
                half_len = std::ceil(target_word.size() / 2.);
            } else {
                half_len = int(target_word.size() / 2.);
            }
            half_word = target_word.substr(0, half_len);
            half_word_next = half_word;
            half_word_next[half_len - 1]++;
        } else if (target_word.size() == 1) {
            half_word = target_word;
            half_word_next = half_word;
            half_word_next[0]++;
        } else {
            throw std::invalid_argument("Empty string in input");
        }

        auto l = std::lower_bound(sortedStrings.begin(), sortedStrings.end(), half_word);
        auto r = std::lower_bound(sortedStrings.begin(), sortedStrings.end(), half_word_next);

        std::string target_word_processed = target_word;
        if (reverseStrings) {
            std::reverse(target_word_processed.begin(), target_word_processed.end());
        }

        for (auto it = l; it != r; ++it) {
            std::string snd_word = *it;
            int dist = distBounded(target_word, snd_word, 1);
            if (0 <= dist && dist <= 1) {
                std::string snd_word_processed = snd_word;
                if (reverseStrings) {
                    std::reverse(snd_word_processed.begin(), snd_word_processed.end());
                }
                for (int snd_idx: str2idx[snd_word_processed])
                    for (int target_idx: str2idx[target_word_processed])
                        out(target_idx, snd_idx) = 1;
            }
        }
    }
}


// [[Rcpp::export(".sortAdjacencyMatSparse")]]
arma::sp_umat sortAdjacencyMatSparse(
    const std::vector<std::string> strings,
    int cutoff,
    char metric,
    bool drop_deg_zero,
    std::string tempfile
) {
    if (cutoff != 1) {
        throw std::invalid_argument("Cutoff != 1 is not implemented for this method");
    }

    std::unordered_map<std::string, std::vector<int>> str2idx;
    arma::sp_umat out = arma::speye<arma::sp_umat>(strings.size(), strings.size());
    std::function<int(const std::string &, const std::string &, const int&)> distBounded;
    if (metric == 'H') {
        distBounded = hamDistBounded;
    } else if (metric == 'L') {
        distBounded = levDistBoundedMock;
    } else {
        throw std::invalid_argument("Choose metric param from {L, H}");
    }

    for (int i = 0; i < static_cast<int>(strings.size()); i++)
        str2idx[strings[i]].push_back(i);

    std::vector<std::string> sortedStrings;
    for (const auto& pair : str2idx) {
        sortedStrings.push_back(pair.first);
    }
    std::sort(sortedStrings.begin(), sortedStrings.end());
    std::unordered_map<std::string, std::set<std::string>> word2neighbours;
    
    processWords(sortedStrings, str2idx, out, distBounded, false);

    for (auto &word : sortedStrings) {
        std::reverse(word.begin(), word.end());
    }
    std::sort(sortedStrings.begin(), sortedStrings.end());

    processWords(sortedStrings, str2idx, out, distBounded, true); 

    dropDegreeZero(drop_deg_zero, out, tempfile);
    return out;
}