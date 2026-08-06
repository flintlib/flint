/*
    Copyright (C) 2026 Edgar Costa

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "longlong.h"
#include "fmpz_mat.h"

/*
    Defensive iteration ceiling for the alternating row/column HNF
    loop in the iterative SNF algorithm.

    Kannan & Bachem (1979), "Polynomial Algorithms for Computing
    the Smith and Hermite Normal Forms of an Integer Matrix",
    SIAM J. Comput. 8, Theorem 5, gives a worst case of
    O(mmn^2 * log(mmn * ||A||)) HNF/LHNF passes with
    mmn = max(m, n).  We use a 4x padding and saturate at
    WORD_MAX so only true runaway loops trigger the bound.
*/
slong
_fmpz_mat_snf_iter_bound(const fmpz_mat_t A)
{
    slong mmn = FLINT_MAX(fmpz_mat_nrows(A), fmpz_mat_ncols(A));
    slong bits = FLINT_ABS(fmpz_mat_max_bits(A));
    slong log_mmn = FLINT_BIT_COUNT((ulong) mmn) + 1;
    double bound = 16.0 + 4.0 * (double) mmn * (double) mmn
                   * (double) (bits + log_mmn);

    return (bound >= (double) WORD_MAX) ? WORD_MAX : (slong) bound;
}
