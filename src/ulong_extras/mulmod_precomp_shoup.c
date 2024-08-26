/*
    Copyright (C) 2015 William Hart
    Copyright (C) 2015 Vladimir Glazachev
    Copyright (C) 2024 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"


/*-------------------------------------------------------------*/
/* notes on Shoup's modular multiplication with precomputation */
/*-------------------------------------------------------------*/

// Sources:
// - NTL code (in particular file FFT.cpp, consulted in NTL v11.5.1)
// - Victor Shoup, "Arithmetic Software Libraries", 2021 (https://doi.org/10.1017/9781108854207.012)
//   (chapter 9 of the book https://doi.org/10.1017/9781108854207 )

// below, B == FLINT_BITS, W == 2**B

// input: a, b, n
// output: (a*b) mod n
// constraints: a < n, and n has < B bits, i.e. 0 <= a < n < 2**(B-1)
//              (there is no restriction on b)
//
// This is intended for repeated multiplications with fixed a and n, and varying b;
// hence the initial step is seen as a precomputation (done only once, depends on a and n only)

// PRECOMPUTATION:
//
// [ ulong n_mulmod_precomp_shoup(ulong a, ulong n) ]
//
// We compute a_precomp = floor(a * W / n),
// This requires a < n, for use of udiv_qrnnd (as in n_mulmod_precomp_shoup) and
// to ensure a_precomp fits in a word (i.e. a_precomp < W).

// MAIN COMPUTATION:
//
// [ ulong n_mulmod_shoup(ulong a, ulong b, ulong a_precomp, ulong n) ]
//
// 1. a_precomp * b = p_hi * W + p_lo     (high part of double-word multiplication, p_lo will not be used)
// 2. res = a*b - p_hi*n                  (single-word multiplications)
// 3. if res >= n, return res-n, else return res
//
// Justifications:
//
// Step 1. One has p_hi = floor(a*b/n) or p_hi = floor(a*b/n) - 1.
//
// Proof: Write a*W = a_precomp * n + r, with 0 <= r < n.
// Thus a_precomp * b / W = a*b/n - r*b/(n*W)
// So p_hi = floor(a_precomp * b / W) = floor(a*b/n - r*b/(n*W))
// Clearly p_hi <= floor(a*b/n)
// And r*b < n*W (since r < n and b < W) so -r*b/(n*W) > -1, hence p_hi >= floor(a*b/n) - 1.
//
// Step 2. It follows that either res == a*b mod n or res == (a*b mod n) + n.
//
// This is where the restriction on n comes into play. If we allowed n to have
// B bits, how to detect which case we are in? A comparison such as res >= n
// would not tell. On the other hand, if n has < B bits, res >= n if and only if
// res == (a*b mod n) + n.
//
// Step 3. we detect which case we are in and correct the possible excess.

