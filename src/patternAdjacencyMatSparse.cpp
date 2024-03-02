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

#include <iostream>
#include <string>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#define ARMA_64BIT_WORD 1
#include "dropDegreeZero.h"
#include <RcppArmadillo.h>
#include <Rcpp.h>
// [[Rcpp::depends(RcppArmadillo)]]

inline std::unordered_set<std::string> getHamming1Patterns(
    const std::string& str,
    std::unordered_set<std::string>* patterns = nullptr
) {
  if (patterns == nullptr)
    patterns = new std::unordered_set<std::string>();
  std::string pattern;
  for (int i = 0; i < static_cast<int>(str.length()); i++) {
    pattern = str;
    pattern[i] = '_';
    patterns->insert(pattern);
  }
  pattern = str;
  pattern.push_back('_');
  patterns->insert(pattern);
  return *patterns;
}

inline std::unordered_set<std::string> getHamming2Patterns(
    const std::string& str,
    std::unordered_set<std::string>* patterns = nullptr
) {
  if (patterns == nullptr)
    patterns = new std::unordered_set<std::string>();
  std::string pattern;
  for (int i = 0; i < static_cast<int>(str.length()); i++) {
    for (int j = i + 1; j < static_cast<int>(str.length()); j++) {
      pattern = str;
      pattern[i] = pattern[j] = '_';
      patterns->insert(pattern);
      pattern = str;
      pattern[i] = '_';
      pattern.push_back('_');
      patterns->insert(pattern);
    }
  }
  pattern = str;
  pattern.push_back('_');
  pattern.push_back('_');
  patterns->insert(pattern);
  pattern = str;
  pattern[static_cast<int>(str.length()) - 1] = '_';
  pattern.push_back('_');
  patterns->insert(pattern);
  getHamming1Patterns(str, patterns);
  return *patterns;
}

inline std::unordered_set<std::string> getLevi1Patterns(
    const std::string& str,
    std::unordered_set<std::string>* patterns = nullptr
) {
  if (patterns == nullptr)
    patterns = new std::unordered_set<std::string>();
  std::string pattern;
  for (int i = 0; i < static_cast<int>(str.length()); i++) {
    pattern = str;
    pattern[i] = '_';
    patterns->insert(pattern);

    pattern = str;
    pattern.insert(i, 1, '_');
    patterns->insert(pattern);
  }
  pattern = str;
  pattern.push_back('_');
  patterns->insert(pattern);
  return *patterns;
}

inline std::unordered_set<std::string> getLevi2Patterns(
    const std::string& str,
    std::unordered_set<std::string>* patterns = nullptr
) {
  if (patterns == nullptr)
    patterns = new std::unordered_set<std::string>();
  std::string pattern;
  for (int i = 0; i < static_cast<int>(str.length()); i++) {
    for (int j = 0; j < i; j++) {
      pattern = str;
      pattern.insert(j, 1, '_');
      pattern[i + 1] = '_';
      patterns->insert(pattern); // k + 1
    }
    for (int j = i; j < static_cast<int>(str.size()); j++) {
      if (j > i) {
        pattern = str;
        pattern[i] = '_';
        pattern[j] = '_';
        patterns->insert(pattern); // k
      }
      pattern = str;
      pattern[i] = '_';
      pattern.insert(j + 1, 1, '_');
      patterns->insert(pattern); // k + 1
      pattern = str;
      pattern.insert(i, 1, '_');
      pattern.insert(j + 1, 1, '_');
      patterns->insert(pattern); // k + 2
    }
    pattern = str;
    pattern.insert(i, 1, '_');
    pattern.push_back('_');
    patterns->insert(pattern); // k + 2
  }
  pattern = str;
  pattern.push_back('_');
  pattern.push_back('_');
  patterns->insert(pattern); // k + 2
  getLevi1Patterns(str, patterns);
  return *patterns;
}

// [[Rcpp::export(".patAdjacencyMatSparse")]]
arma::sp_umat buildG(
    const std::vector<std::string>& strings,
    int cutoff,
    char metric,
    bool drop_deg_zero,
    std::string tempfile
) {

  std::unordered_map<std::string, std::vector<std::string>> pat2str;
  std::function<std::unordered_set<std::string>(const std::string&, std::unordered_set<std::string>*)> getPatternFunc;
  std::unordered_map<std::string, std::vector<int>> str2idx;
  arma::sp_umat out = arma::speye<arma::sp_umat>(strings.size(), strings.size());
  std::unordered_set<std::string> uniqueStrings;
  std::string str, str1, str2;
  int idx1, idx2, size;

  // remember the original indices of strings
  for (int i = 0; i < static_cast<int>(strings.size()); i++)
    str2idx[strings[i]].push_back(i);

  if (cutoff == 0) {
    // for distance = 0 it is enough to have str2idx, uniqueStrings
    for (const auto &pair: str2idx) {
      str = pair.first;
      size = str2idx[str].size();
      for (int i1 = 0; i1 < size; i1++) {
        idx1 = str2idx[str][i1];
        for (int i2 = i1; i2 < size; i2++) {
          idx2 = str2idx[str][i2];
          out(idx1, idx2) = out(idx2, idx1) = 1;
        }
      }
    }
  } else {
    // choose the type of edit distance
    if (cutoff == 1 && metric == 'H') {
      getPatternFunc = getHamming1Patterns;
    } else if (cutoff == 2 && metric == 'H') {
      getPatternFunc = getHamming2Patterns;
    } else if (cutoff == 1 && metric == 'L') {
      getPatternFunc = getLevi1Patterns;
    } else if (cutoff == 2 && metric == 'L') {
      getPatternFunc = getLevi2Patterns;
    } else {
      throw std::invalid_argument("Choose metric param from {L, H} and cutoff from {0, 1, 2}");
    }

    // collect the set of unique strings
    for (const std::string& str : strings)
      uniqueStrings.insert(str);

    // calculate the pattern->strings map (for the set of unique strings)
    for (const std::string& str : uniqueStrings) {
      Rcpp::checkUserInterrupt();
      std::unordered_set<std::string> patterns = getPatternFunc(str, nullptr);
      for (const std::string& pattern : patterns) {
        pat2str[pattern].push_back(str);
      }
    }

    // calculate the adjecent matrix for unique vector of strings
    for (const auto &entry: pat2str) {
      Rcpp::checkUserInterrupt();
      if (static_cast<int>(entry.second.size()) > 1 || static_cast<int>(str2idx[entry.second[0]].size()) > 1) {
        // iterate over the pairs of unique strings corresponding to a certain pattern
        for (int i = 0; i < static_cast<int>(entry.second.size()); i++) {
          str1 = entry.second[i];
          for (int j = i; j < static_cast<int>(entry.second.size()); j++) {
            str2 = entry.second[j];
            // iterate over the pairs of all original strings
            // corresponding to the pair of unique strings
            for (const auto &idx1: str2idx[str1]) {
              for (const auto &idx2: str2idx[str2]) {
                out(idx1, idx2) = out(idx2, idx1) = 1;
              }
            }
          }
        }
      }
    }
  }

  dropDegreeZero(drop_deg_zero, out, tempfile);
  return out;
}