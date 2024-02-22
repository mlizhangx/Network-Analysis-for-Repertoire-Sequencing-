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
#include <Rcpp.h>

void dropDegreeZero(
    bool drop_deg_zero,
    arma::sp_umat& out,
    std::string tempfile
) {
    if (drop_deg_zero) {
        // sum entries columnwise
        arma::sp_umat col_sums_spmat = arma::sum(out);
        arma::urowvec col_sums(col_sums_spmat);

        // record indices of nodes with positive degree
        arma::uvec col_ids = find(col_sums > 1);

        // subset matrix to keep only network nodes
        out = out.cols(col_ids);
        out = out.t();
        out = out.cols(col_ids);

        // write indices of network nodes to file
        col_ids += 1;  // offset C++'s 0-index starting convention
        col_ids.save(tempfile, arma::raw_ascii);
    }
}
