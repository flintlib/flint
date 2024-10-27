/*
    Copyright (C) 2024 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef N_FFT_BASIC_H
#define N_FFT_BASIC_H

/*---------*/
/* helpers */
/*---------*/

/** Shoup's modular multiplication with precomputation, lazy
 * (does not perform the excess correction step)
 *  --> computes either r or r+n and store it is res, where r = (a*b) % n
 *  --> a_pr is the precomputation for n, p_hi and p_lo are temporaries
 *  --> requires nbits(n) < FLINT_BITS
 */
#define N_MULMOD_PRECOMP_LAZY(res, a, b, a_pr, n, p_hi, p_lo) \
    do {                                                      \
        umul_ppmm(p_hi, p_lo, (a_pr), (b));                   \
        res = (a) * (b) - p_hi * (n);                         \
    } while(0)

/*-------------*/
/* 2-point DFT */
/*-------------*/

/** Butterfly, node 0
 * * in [0..n) x [0..n) / out [0..2n) x [0..2n) / max < 2n
 * * In-place transform
 *                            [1  1]
 *           [a  b] <- [a  b] [1 -1]
 * * n is the modulus, tmp is a temporary
 */
#define DFT2_NODE0_LAZY12(a, b, n, tmp) \
    do {                                \
        tmp = (b);                      \
        (b) = (a) + (n) - tmp;          \
        (a) = (a) + tmp;                \
    } while(0)

/** Butterfly, node 0
 * * in [0..2n) x [0..2n) / out [0..2n) x [0..4n) / max < 4n
 * * In-place transform
 *                            [1  1]
 *           [a  b] <- [a  b] [1 -1]
 * * n2 is 2*n, tmp is a temporary
 */
#define DFT2_NODE0_LAZY24(a, b, n2, tmp) \
    do {                               \
        tmp = (b);                     \
        (b) = (a) + (n2) - tmp;        \
        (a) = (a) + tmp;               \
        if ((a) >= (n2))               \
            (a) -= (n2);               \
    } while(0)

#endif  // N_FFT_BASIC_H
