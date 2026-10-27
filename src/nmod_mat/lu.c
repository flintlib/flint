/*
    Copyright (C) 2011, 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod.h"
#include "nmod_vec.h"
#include "nmod_mat.h"
#include "nmod_mat/impl.h"

/*
    Dimension min(m, n) from which nmod_mat_lu_recursive is used rather than
    nmod_mat_lu_classical_delayed, by class (AVX512-IFMA, NEON, others) and
    modulus size. Fitted on the crossovers between the two at the top level,
    on uniformly random matrices (squares 32 to 1024 and ratios 2 to 8), on
    Zen 4 and Ice Lake SP (IFMA), Apple M4 (NEON) and Cascade Lake
    (others), then raised where the recursion was slower on matrices of
    rank n / 2 (up to 29 bits: from about 96 to 192). On x86, from 31 to 52
    bits, the recursion pays off from small dimensions on.
*/
static slong
lu_rec_cutoff(flint_bitcnt_t bits)
{
#if NMOD_MAT_HAVE_MUL_U52
    return (bits <= 2) ? 96 : (bits <= 29) ? 192 : (bits <= 30) ? 64
         : (bits <= 52) ? 48 : 128;
#elif (defined(__ARM_NEON) && defined(__aarch64__)) || defined(_M_ARM64)
    return (bits <= 2) ? 96 : (bits <= 29) ? 192 : (bits <= 32) ? 64 : 192;
#else
    return (bits <= 2) ? 64 : (bits <= 29) ? 192 : (bits <= 52) ? 64 : 128;
#endif
}

slong
nmod_mat_lu(slong * P, nmod_mat_t A, int rank_check)
{
    slong nrows = A->r, ncols = A->c, n;
    n = FLINT_MIN(nrows, ncols);

    if (n <= 3 || (NMOD_BITS(A->mod) > 28 && n <= 7))
    {
        return nmod_mat_lu_classical(P, A, rank_check);
    }
    else
    {
        if (n < lu_rec_cutoff(NMOD_BITS(A->mod)))
            return nmod_mat_lu_classical_delayed(P, A, rank_check);
        else
            return nmod_mat_lu_recursive(P, A, rank_check);
    }
}

