/*
    Copyright (C) 2022 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fft_small.h"

/*
    This module contains the function `sd_fft_trunc` to compute a truncated FFT,
    together with the helper functions and macros to implement it.
    See [vdH2004]_ for a high-level overview of the truncated FFT algorithm.

    For a concrete example, let `L = 3` and `itrunc = otrunc = 8`, so no truncation
    happens. The algorithm can be visualized by the following diagram:

         0   1   2   3   4   5   6   7
        [-------------]X[-------------]  // layer 0
         0   1   2   3   4   5   6   7
        [-----]X[-----] [-----]X[-----]  // layer 1
         0   1   2   3   4   5   6   7
        [-]X[-] [-]X[-] [-]X[-] [-]X[-]  // layer 2
         0   1   2   3   4   5   6   7

    We say:
    - at layer 0, there is 1 transform of length 8,
    - at layer 1, there are 2 independent contiguous transforms each of length 4,
    - at layer 2, there are 4 independent contiguous transforms each of length 2.

    Let `w` be an 8-th root of unity. If the input array represents a polynomial
    `f(x) = sum(a[i] * x^i)`, then the output array is

        [f(w^0), f(w^4), f(w^2), f(w^6), f(w^1), f(w^5), f(w^3), f(w^7)].

    Note that the order is bit-flipped, but this is fine for us. The goal
    is to compute the convolution of two sequences, so it suffices for the FFT
    result to be permuted, as long as the inverse FFT takes into account the same
    permutation. This does not affect the pointwise multiplication either.

    The bit-flipped ordering is actually more convenient. To see why, we redraw
    the diagram above with each transform labeled with the twiddle factor:

         0   1   2   3   4   5   6   7
        [------------ w^0 ------------]  // layer 0
         0   1   2   3   4   5   6   7
        [---- w^0 ----] [---- w^2 ----]  // layer 1
         0   1   2   3   4   5   6   7
        [ w^0 ] [ w^2 ] [ w^1 ] [ w^3 ]  // layer 2
         0   1   2   3   4   5   6   7

    This is exactly the ordering of the stored even-index twiddle factors `sd_fft_ctx_w2(Q, j)`.

    Here, the diagram

         a   b
        [ w^i ]
         c   d

    corresponds to the (in-place) calculation `c = a + b * w^i`, `d = a - b * w^i`.

    In the diagram above, layer 1 is shown as 2 independent contiguous transforms
    of length 4, while layer 0 is shown as 1 contiguous transform of length 8.
    However, if we consider the flat array as a 2x4 matrix in row-major order:

        [  0   1   2   3  ]
        [                 ]
        [  4   5   6   7  ]

    the layer 0 transform can be seen as 4 independent transforms operating
    on each _column_:

        [  0   1   2   3  ]
        [  |   |   |   |  ]
        [  4   5   6   7  ]

    while the layer 1 transform can be seen as 2 independent transforms operating
    on each row:

        [  0---1---2---3  ]
        [                 ]
        [  4---5---6---7  ]

    That's why "column" and "row" are used in the comments below.

    The exact transform above is performed by `sd_fft_basecase_3_1`.
    We call this **one contiguous 3-layer transform with lengths (8, 4, 2)**.
    This is done by one contiguous 1-layer transform with length 8, followed by
    two independent contiguous 2-layer transforms with lengths (4, 2).
    In general,

        A independent contiguous B-layer transforms with lengths (C, C/2, ..., C/2^(B-1))

    means

         0   1  ...   C-1   C   C+1 ...  2C-1    ........   (A-1)C (A-1)C+1 ...  AC-2   AC-1
        [--------X-------] [---------X-------] [-] ... [-] [-----------------X--------------]  // layer 0
         0   1  ...   C-1   C   C+1 ...  2C-1    ........   (A-1)C (A-1)C+1 ...  AC-2   AC-1
        [---X---] [--X---] [---X----] [--X---] [-] ... [-] [-------X--------] [-----X-------]  // layer 1
         0   1  ...   C-1   C   C+1 ...  2C-1    ........   (A-1)C (A-1)C+1 ...  AC-2   AC-1
          .
          .
          .
        [---X---] [--X---] [---X----] [--X---] [-] ... [-] [-------X--------] [-----X-------]  // layer B-1
         0   1  ...   C-1   C   C+1 ...  2C-1    ........   (A-1)C (A-1)C+1 ...  AC-2   AC-1

    Here, "contiguous" emphasizes that each independent transform operates on a contiguous
    block of memory. The "column" transforms are also independent, but not contiguous.

    When `L >= 4`, the order is less consistent: for `w` a 16-th root of unity, the
    output array is not

        [f(w^0),  f(w^8),  f(w^4),  f(w^12),
         f(w^2),  f(w^10), f(w^6),  f(w^14),
         f(w^1),  f(w^9),  f(w^5),  f(w^13),
         f(w^3),  f(w^11), f(w^7),  f(w^15)]

    but it is instead

        [f(w^0),  f(w^2),  f(w^1),  f(w^3),
         f(w^8),  f(w^10), f(w^9),  f(w^11),
         f(w^4),  f(w^6),  f(w^5),  f(w^7),
         f(w^12), f(w^14), f(w^13), f(w^15)]

    Formally, each block of 16 consecutive terms is transposed as a 4x4 matrix.

    See the comment in `sd_fft_basecase_4_1`: a final `VEC4D_TRANSPOSE` is omitted.
    This also explains why `sd_fft_ctx_trunc_index` needs to adjust the last 4 bits of `i`.

    Because the column transforms are also independent, it is not necessary
    to handle the earlier layers strictly before the later layers either.
    For example, to do a 2-layer transform with lengths (16, 8):

         0   1   2   3   4   5   6   7   8   9   10  11  12  13  14  15
        [---------------------------- w^0 ----------------------------]
         0   1   2   3   4   5   6   7   8   9   10  11  12  13  14  15
        [------------ w^0 ------------] [------------ w^4 ------------]
         0   1   2   3   4   5   6   7   8   9   10  11  12  13  14  15

    we can first do ("#" marks the transformed entries, "-" marks the
    unmodified entries):

         0   1   2   3   4   5   6   7   8   9   10  11  12  13  14  15
        [#######--------########--------########--------########------]
         0   1   2   3   4   5   6   7   8   9   10  11  12  13  14  15
        [#######--------########------] [#######--------########------]
         0   1   2   3   4   5   6   7   8   9   10  11  12  13  14  15

    followed by:

         0   1   2   3   4   5   6   7   8   9   10  11  12  13  14  15
        [-------#######---------#######---------#######---------######]
         0   1   2   3   4   5   6   7   8   9   10  11  12  13  14  15
        [-------#######---------######] [-------#######---------######]
         0   1   2   3   4   5   6   7   8   9   10  11  12  13  14  15

    We call the two steps above a

        2-layer transform with lengths (16, 8) and mask (0..2)

    and a

        2-layer transform with lengths (16, 8) and mask (2..4)

    respectively. More generally, we have a

        B-layer transform with lengths (C, C/2, ..., C/2^(B-1)) and mask (D..E)

    for 0 <= D < E <= C/2^B. The mask range is half-open: only entries
    whose index is congruent to some `i` with D <= i < E modulo C/2^B are
    transformed.

    ********************************************************************************

    N is supposed to be a good fit for the number of points to process per loop
    in the radix 4 butterflies.

        16x 4-wide AVX registers  => N = 8
        32x 2-wide NEON registers => N = 8
        32x 8-wide AVX512 registers  => N = ?
*/

#define N 8
#define VECND vec8d
#define VECNOP(op) TEMPLATE(VECND, op)

/********************* forward butterfly **************************************
    b0 = a0 + w*a1
    b1 = a0 - w*a1

    Many macros here come in pairs; each pair has a `PARAM` macro and a `MOTH` macro.
    They are meant to be used as::

        *_PARAM_*(...);
        ulong i = 0; do {
            *_MOTH_*(...);
        } while (i += N, i < BLK_SZ);

    The `PARAM` macro defines local variables used by the corresponding `MOTH` macro.

    The first argument `V` is usually `VECND` = `vec8d`. This means each `MOTH`
    macro invocation transforms `N = 8` entries from each `X*` argument,
    so `i` is increased by `N` each iteration.

    `j` is the twiddle table index. In the example above, we have:

         0   1   2   3   4   5   6   7
        [------------ j=0 ------------]  // layer 0
         0   1   2   3   4   5   6   7
        [---- j=0 ----] [---- j=1 ----]  // layer 1
         0   1   2   3   4   5   6   7
        [ j=0 ] [ j=1 ] [ j=2 ] [ j=3 ]  // layer 2
         0   1   2   3   4   5   6   7

    When `j == 0`, there is no need to multiply by the twiddle factor. This is why
    there is separate code for the `j == 0` (`J_IS_Z`) and `j != 0` cases.

    When the macros take separate `j_r` and `j_bits`, these values can be
    computed from `j` with the `SET_J_BITS_AND_J_R` macro.

    The following macros perform a 1-layer transform. The mask length of each
    transform is determined by the outer loop, in the example above, it is `BLK_SZ`.
    The spacing between X0 and X1 is half transform length.
*/

#define RADIX_2_FORWARD_PARAM_J_IS_Z(V, Q) \
    V n    = CAT(V, set_d)(Q->p); \
    V ninv = CAT(V, set_d)(Q->pinv);

#define RADIX_2_FORWARD_MOTH_J_IS_Z(V, X0, X1) \
{ \
    V x0, x1; \
    x0 = CAT(V, load)(X0); \
    x0 = CAT(V, reduce_to_pm1n)(x0, n, ninv); \
    x1 = CAT(V, load)(X1); \
    x1 = CAT(V, reduce_to_pm1n)(x1, n, ninv); \
    CAT(V, store)(X0, CAT(V, add)(x0, x1)); \
    CAT(V, store)(X1, CAT(V, sub)(x0, x1)); \
}

#define RADIX_2_FORWARD_PARAM_J_IS_NZ(V, Q, j_r, j_bits) \
    V w = CAT(V, set_d)(Q->w2tab[j_bits][j_r]); \
    V n    = CAT(V, set_d)(Q->p); \
    V ninv = CAT(V, set_d)(Q->pinv);

#define RADIX_2_FORWARD_MOTH_J_IS_NZ(V, X0, X1) \
{ \
    V x0, x1; \
    x0 = CAT(V, load)(X0); \
    x0 = CAT(V, reduce_to_pm1n)(x0, n, ninv); \
    x1 = CAT(V, load)(X1); \
    x1 = CAT(V, mulmod)(x1, w, n, ninv); \
    CAT(V, store)(X0, CAT(V, add)(x0, x1)); \
    CAT(V, store)(X1, CAT(V, sub)(x0, x1)); \
}

/**************** forward butterfly with truncation **************************/

#define RADIX_2_FORWARD_MOTH_TRUNC_2_1_J_IS_Z(V, X0, X1) \
{ \
    V x0, x1; \
    x0 = CAT(V, load)(X0); \
    x0 = CAT(V, reduce_to_pm1n)(x0, n, ninv); \
    x1 = CAT(V, load)(X1); \
    x1 = CAT(V, reduce_to_pm1n)(x1, n, ninv); \
    CAT(V, store)(X0, CAT(V, add)(x0, x1)); \
}

#define RADIX_2_FORWARD_MOTH_TRUNC_2_1_J_IS_NZ(V, X0, X1) \
{ \
    V x0, x1; \
    x0 = CAT(V, load)(X0); \
    x0 = CAT(V, reduce_to_pm1n)(x0, n, ninv); \
    x1 = CAT(V, load)(X1); \
    x1 = CAT(V, mulmod)(x1, w, n, ninv); \
    CAT(V, store)(X0, CAT(V, add)(x0, x1)); \
}

/********************* forward butterfly **************************************
    b0 = a0 + w^2*a2 +   w*(a1 + w^2*a3)
    b1 = a0 + w^2*a2 -   w*(a1 + w^2*a3)
    b2 = a0 - w^2*a2 + i*w*(a1 - w^2*a3)
    b3 = a0 - w^2*a2 - i*w*(a1 - w^2*a3)

    In other words: a 2-layer transform.
*/

#define RADIX_4_FORWARD_PARAM_J_IS_Z(V, Q) \
    V iw = CAT(V, set_d)(Q->w2tab[1][0]); \
    V n    = CAT(V, set_d)(Q->p); \
    V ninv = CAT(V, set_d)(Q->pinv);

#define RADIX_4_FORWARD_MOTH_J_IS_Z(V, X0, X1, X2, X3) \
{ \
    V x0, x1, x2, x3, y0, y1, y2, y3; \
    x0 = CAT(V, load)(X0); \
    x0 = CAT(V, reduce_to_pm1n)(x0, n, ninv); \
    x1 = CAT(V, load)(X1); \
    x2 = CAT(V, load)(X2); \
    x3 = CAT(V, load)(X3); \
    x2 = CAT(V, reduce_to_pm1n)(x2, n, ninv); \
    x3 = CAT(V, reduce_to_pm1n)(x3, n, ninv); \
    y0 = CAT(V, add)(x0, x2); \
    y1 = CAT(V, add)(x1, x3); \
    y2 = CAT(V, sub)(x0, x2); \
    y3 = CAT(V, sub)(x1, x3); \
    y1 = CAT(V, reduce_to_pm1n)(y1, n, ninv); \
    y3 = CAT(V, mulmod)(y3, iw, n, ninv); \
    x0 = CAT(V, add)(y0, y1); \
    x1 = CAT(V, sub)(y0, y1); \
    x2 = CAT(V, add)(y2, y3); \
    x3 = CAT(V, sub)(y2, y3); \
    CAT(V, store)(X0, x0); \
    CAT(V, store)(X1, x1); \
    CAT(V, store)(X2, x2); \
    CAT(V, store)(X3, x3); \
}

#define RADIX_4_FORWARD_PARAM_J_IS_NZ(V, Q, j_r, j_bits) \
    FLINT_ASSERT(j_bits > 0); \
    V w  = CAT(V, set_d)(Q->w2tab[1+j_bits][2*j_r]); \
    V w2 = CAT(V, set_d)(Q->w2tab[0+j_bits][j_r]); \
    V iw = CAT(V, set_d)(Q->w2tab[1+j_bits][2*j_r+1]); \
    V n    = CAT(V, set_d)(Q->p); \
    V ninv = CAT(V, set_d)(Q->pinv);

#define RADIX_4_FORWARD_MOTH_J_IS_NZ(V, X0, X1, X2, X3) \
{ \
    V x0, x1, x2, x3, y0, y1, y2, y3; \
    x0 = CAT(V, load)(X0); \
    x0 = CAT(V, reduce_to_pm1n)(x0, n, ninv); \
    x1 = CAT(V, load)(X1); \
    x2 = CAT(V, load)(X2); \
    x3 = CAT(V, load)(X3); \
    x2 = CAT(V, mulmod)(x2, w2, n, ninv); \
    x3 = CAT(V, mulmod)(x3, w2, n, ninv); \
    y0 = CAT(V, add)(x0, x2); \
    y1 = CAT(V, add)(x1, x3); \
    y2 = CAT(V, sub)(x0, x2); \
    y3 = CAT(V, sub)(x1, x3); \
    y1 = CAT(V, mulmod)(y1, w, n, ninv); \
    y3 = CAT(V, mulmod)(y3, iw, n, ninv); \
    x0 = CAT(V, add)(y0, y1); \
    x1 = CAT(V, sub)(y0, y1); \
    x2 = CAT(V, add)(y2, y3); \
    x3 = CAT(V, sub)(y2, y3); \
    CAT(V, store)(X0, x0); \
    CAT(V, store)(X1, x1); \
    CAT(V, store)(X2, x2); \
    CAT(V, store)(X3, x3); \
}

#define LENGTH4_ANY_J(T, x0, x1, x2, x3, n, ninv, w2, w, iw) \
{ \
    T X0 = x0, X1 = x1, X2 = x2, X3 = x3, Y0, Y1, Y2, Y3; \
    X0 = CAT(T, reduce_to_pm1n)(X0, n, ninv); \
    X2 = CAT(T, mulmod)(X2, w2, n, ninv); \
    X3 = CAT(T, mulmod)(X3, w2, n, ninv); \
    Y0 = CAT(T, add)(X0, X2); \
    Y1 = CAT(T, add)(X1, X3); \
    Y2 = CAT(T, sub)(X0, X2); \
    Y3 = CAT(T, sub)(X1, X3); \
    Y1 = CAT(T, mulmod)(Y1, w, n, ninv); \
    Y3 = CAT(T, mulmod)(Y3, iw, n, ninv); \
    x0 = CAT(T, add)(Y0, Y1); \
    x1 = CAT(T, sub)(Y0, Y1); \
    x2 = CAT(T, add)(Y2, Y3); \
    x3 = CAT(T, sub)(Y2, Y3); \
}

#define LENGTH4_ZERO_J(T, x0, x1, x2, x3, n, ninv, e14) \
{ \
    T X0 = x0, X1 = x1, X2 = x2, X3 = x3, Y0, Y1, Y2, Y3; \
    X0 = CAT(T, reduce_to_pm1n)(X0, n, ninv); \
    X2 = CAT(T, reduce_to_pm1n)(X2, n, ninv); \
    X3 = CAT(T, reduce_to_pm1n)(X3, n, ninv); \
    Y0 = CAT(T, add)(X0, X2); \
    Y1 = CAT(T, add)(X1, X3); \
    Y2 = CAT(T, sub)(X0, X2); \
    Y3 = CAT(T, sub)(X1, X3); \
    Y1 = CAT(T, reduce_to_pm1n)(Y1, n, ninv); \
    Y3 = CAT(T, mulmod)(Y3, e14, n, ninv); \
    x0 = CAT(T, add)(Y0, Y1); \
    x1 = CAT(T, sub)(Y0, Y1); \
    x2 = CAT(T, add)(Y2, Y3); \
    x3 = CAT(T, sub)(Y2, Y3); \
}

/* A 3-layer transform. */

#define LENGTH8_ANY_J(T, x0, x1, x2, x3, x4, x5, x6, x7, n, ninv, w2, w, iw, ww0, ww1, ww2, ww3) \
{ \
    T X0 = x0, X1 = x1, X2 = x2, X3 = x3, X4 = x4, X5 = x5, X6 = x6, X7 = x7; \
    T Y0, Y1, Y2, Y3, Y4, Y5, Y6, Y7, Z0, Z1, Z2, Z3, Z4, Z5, Z6, Z7; \
    X0 = CAT(T, reduce_to_pm1n)(X0, n, ninv); \
    X1 = CAT(T, reduce_to_pm1n)(X1, n, ninv); \
    X4 = CAT(T, mulmod)(X4, w2, n, ninv); \
    X5 = CAT(T, mulmod)(X5, w2, n, ninv); \
    X6 = CAT(T, mulmod)(X6, w2, n, ninv); \
    X7 = CAT(T, mulmod)(X7, w2, n, ninv); \
    Y0 = CAT(T, add)(X0, X4); \
    Y1 = CAT(T, add)(X1, X5); \
    Y2 = CAT(T, add)(X2, X6); \
    Y3 = CAT(T, add)(X3, X7); \
    Y4 = CAT(T, sub)(X0, X4); \
    Y5 = CAT(T, sub)(X1, X5); \
    Y6 = CAT(T, sub)(X2, X6); \
    Y7 = CAT(T, sub)(X3, X7); \
    Y2 = CAT(T, mulmod)(Y2, w, n, ninv); \
    Y3 = CAT(T, mulmod)(Y3, w, n, ninv); \
    Y6 = CAT(T, mulmod)(Y6, iw, n, ninv); \
    Y7 = CAT(T, mulmod)(Y7, iw, n, ninv); \
    Z0 = CAT(T, add)(Y0, Y2); \
    Z1 = CAT(T, add)(Y1, Y3); \
    Z2 = CAT(T, sub)(Y0, Y2); \
    Z3 = CAT(T, sub)(Y1, Y3); \
    Z4 = CAT(T, add)(Y4, Y6); \
    Z5 = CAT(T, add)(Y5, Y7); \
    Z6 = CAT(T, sub)(Y4, Y6); \
    Z7 = CAT(T, sub)(Y5, Y7); \
    Z0 = CAT(T, reduce_to_pm1n)(Z0, n, ninv); \
    Z1 = CAT(T, mulmod)(Z1, ww0, n, ninv); \
    Z2 = CAT(T, reduce_to_pm1n)(Z2, n, ninv); \
    Z3 = CAT(T, mulmod)(Z3, ww1, n, ninv); \
    Z4 = CAT(T, reduce_to_pm1n)(Z4, n, ninv); \
    Z5 = CAT(T, mulmod)(Z5, ww2, n, ninv); \
    Z6 = CAT(T, reduce_to_pm1n)(Z6, n, ninv); \
    Z7 = CAT(T, mulmod)(Z7, ww3, n, ninv); \
    x0 = CAT(T, add)(Z0, Z1); \
    x1 = CAT(T, sub)(Z0, Z1); \
    x2 = CAT(T, add)(Z2, Z3); \
    x3 = CAT(T, sub)(Z2, Z3); \
    x4 = CAT(T, add)(Z4, Z5); \
    x5 = CAT(T, sub)(Z4, Z5); \
    x6 = CAT(T, add)(Z6, Z7); \
    x7 = CAT(T, sub)(Z6, Z7); \
}

#define LENGTH8_ZERO_J(T, x0, x1, x2, x3, x4, x5, x6, x7, n, ninv, e14, e18, e38) \
{ \
    T X0 = x0, X1 = x1, X2 = x2, X3 = x3, X4 = x4, X5 = x5, X6 = x6, X7 = x7; \
    T Y0, Y1, Y2, Y3, Y4, Y5, Y6, Y7, Z0, Z1, Z2, Z3, Z4, Z5, Z6, Z7; \
    Y0 = CAT(T, reduce_to_pm1n)(CAT(T, add)(X0, X4), n, ninv); \
    Y1 = CAT(T, reduce_to_pm1n)(CAT(T, add)(X1, X5), n, ninv); \
    Y2 = CAT(T, reduce_to_pm1n)(CAT(T, add)(X2, X6), n, ninv); \
    Y3 = CAT(T, reduce_to_pm1n)(CAT(T, add)(X3, X7), n, ninv); \
    Y4 = CAT(T, reduce_to_pm1n)(CAT(T, sub)(X0, X4), n, ninv); \
    Y5 = CAT(T, reduce_to_pm1n)(CAT(T, sub)(X1, X5), n, ninv); \
    Y6 = CAT(T, reduce_to_pm1n)(CAT(T, sub)(X2, X6), n, ninv); \
    Y7 = CAT(T, reduce_to_pm1n)(CAT(T, sub)(X3, X7), n, ninv); \
    Z0 = CAT(T, add)(Y0, Y2); \
    Z1 = CAT(T, add)(Y1, Y3); \
    Z2 = CAT(T, sub)(Y0, Y2); \
    Z3 = CAT(T, sub)(Y1, Y3); \
    Y6 = CAT(T, mulmod)(e14, Y6, n, ninv); \
    Y7 = CAT(T, mulmod)(e14, Y7, n, ninv); \
    Z4 = CAT(T, add)(Y4, Y6); \
    Z5 = CAT(T, add)(Y5, Y7); \
    Z6 = CAT(T, sub)(Y4, Y6); \
    Z7 = CAT(T, sub)(Y5, Y7); \
    x0 = CAT(T, add)(Z0, Z1); \
    x1 = CAT(T, sub)(Z0, Z1); \
    Z3 = CAT(T, mulmod)(e14, Z3, n, ninv); \
    Z5 = CAT(T, mulmod)(e18, Z5, n, ninv); \
    Z7 = CAT(T, mulmod)(e38, Z7, n, ninv); \
    x2 = CAT(T, add)(Z2, Z3); \
    x3 = CAT(T, sub)(Z2, Z3); \
    x4 = CAT(T, add)(Z4, Z5); \
    x5 = CAT(T, sub)(Z4, Z5); \
    x6 = CAT(T, add)(Z6, Z7); \
    x7 = CAT(T, sub)(Z6, Z7); \
}


/**************** basecase transform of size 2^m **********************/
/* template<int m, bool j_is_zero> sd_fft_basecase(Q, X, j_r, j_bits)

    With notation as above: sd_fft_basecase_{m} performs a contiguous m-layer
    transform with lengths (2^m, 2^(m-1), ..., 2^1).  */

static void sd_fft_basecase_0_1(const sd_fft_ctx_t FLINT_UNUSED(Q), double* FLINT_UNUSED(X))
{
}


static void sd_fft_basecase_1_1(const sd_fft_ctx_t Q, double* X)
{
    double n    = Q->p;
    double ninv = Q->pinv;
    double x0 = vec1d_reduce_to_pm1n(X[0], n, ninv);
    double x1 = vec1d_reduce_to_pm1n(X[1], n, ninv);
    X[0] = vec1d_add(x0, x1);
    X[1] = vec1d_sub(x0, x1);
}


static void sd_fft_basecase_2_1(const sd_fft_ctx_t Q, double* X)
{
    LENGTH4_ZERO_J(vec1d, X[0], X[1], X[2], X[3], Q->p, Q->pinv, Q->w2tab[1][0])
}


static void sd_fft_basecase_3_1(const sd_fft_ctx_t Q, double* X)
{
    LENGTH8_ZERO_J(vec1d, X[0], X[1], X[2], X[3], X[4], X[5], X[6], X[7], Q->p, Q->pinv,
                   Q->w2tab[1][0], Q->w2tab[2][0], Q->w2tab[2][1])
}


/* missing final transpose gives length >= 16 a worse-than-bit-reversed order */
static void sd_fft_basecase_4_1(const sd_fft_ctx_t Q, double* X)
{
    vec4d n    = vec4d_set_d(Q->p);
    vec4d ninv = vec4d_set_d(Q->pinv);
    vec4d w, w2, iw, u, v;

    vec4d x0 = vec4d_load(X+4*0);
    vec4d x1 = vec4d_load(X+4*1);
    vec4d x2 = vec4d_load(X+4*2);
    vec4d x3 = vec4d_load(X+4*3);

    FLINT_ASSERT(SD_FFT_CTX_W2TAB_INIT >= 4); /* Q.w2tab[0] points to consecutive entries */
    iw = vec4d_set_d(Q->w2tab[0][1]);
    LENGTH4_ZERO_J(vec4d, x0,x1,x2,x3, n,ninv, iw);

    u = vec4d_load_aligned(&Q->w2tab[0][0]);
    v = vec4d_load_aligned(&Q->w2tab[0][4]);
    w2 = u;
    w  = vec4d_unpack_lo_permute_0_2_1_3(u, v);
    iw = vec4d_unpack_hi_permute_0_2_1_3(u, v);
    VEC4D_TRANSPOSE(x0,x1,x2,x3, x0,x1,x2,x3);
    LENGTH4_ANY_J(vec4d, x0,x1,x2,x3, n,ninv, w2, w,iw);

    /* VEC4D_TRANSPOSE(x0,x1,x2,x3, x0,x1,x2,x3) */ /* skipped */
    vec4d_store(X+4*0, x0);
    vec4d_store(X+4*1, x1);
    vec4d_store(X+4*2, x2);
    vec4d_store(X+4*3, x3);
}

#if !defined(__AVX512F__)
static void sd_fft_basecase_4_0(const sd_fft_ctx_t Q, double* X, ulong j_r, ulong j_bits)
{
    vec4d n    = vec4d_set_d(Q->p);
    vec4d ninv = vec4d_set_d(Q->pinv);
    vec4d w, w2, iw, u, v;

    vec4d x0 = vec4d_load(X+4*0);
    vec4d x1 = vec4d_load(X+4*1);
    vec4d x2 = vec4d_load(X+4*2);
    vec4d x3 = vec4d_load(X+4*3);

    w2 = vec4d_set_d(Q->w2tab[0+j_bits][1*j_r+0]);
    w  = vec4d_set_d(Q->w2tab[1+j_bits][2*j_r+0]);
    iw = vec4d_set_d(Q->w2tab[1+j_bits][2*j_r+1]);
    LENGTH4_ANY_J(vec4d, x0,x1,x2,x3, n,ninv, w2, w,iw);

    u  = vec4d_load_aligned(&Q->w2tab[3+j_bits][8*j_r+0]);
    v  = vec4d_load_aligned(&Q->w2tab[3+j_bits][8*j_r+4]);
    w2 = vec4d_load_aligned(&Q->w2tab[2+j_bits][4*j_r+0]);
    w  = vec4d_unpack_lo_permute_0_2_1_3(u, v);
    iw = vec4d_unpack_hi_permute_0_2_1_3(u, v);
    VEC4D_TRANSPOSE(x0,x1,x2,x3, x0,x1,x2,x3);
    LENGTH4_ANY_J(vec4d, x0,x1,x2,x3, n,ninv, w2, w,iw);

    /* VEC4D_TRANSPOSE(x0,x1,x2,x3, x0,x1,x2,x3) */ /* skipped */
    vec4d_store(X+4*0, x0);
    vec4d_store(X+4*1, x1);
    vec4d_store(X+4*2, x2);
    vec4d_store(X+4*3, x3);
}
#endif

/*
The length 32 transform can be broken up as
   (a) 8 transforms of length 4 within columns, followed by 4 transforms of length 8 in the rows, or
   (b) 4 transforms of length 8 within columns, followed by 8 transforms of length 4 in the rows
Since the length 16 basecase is missing the final 4x4 transpose, so the output
is worse than bit-reversed. If the length 32 transform used a order different from 16's,
then we will have a problem at a higher level since it would be difficult to keep track
of what basecase happened to have been used. Therefore, the length 16 and 32 basecases
should produce the same order, and this is easier with (b).
*/
static void sd_fft_basecase_5_1(const sd_fft_ctx_t Q, double* X)
{
    vec4d n    = vec4d_set_d(Q->p);
    vec4d ninv = vec4d_set_d(Q->pinv);
    vec4d u, v, w0, ww0, ww1, www2, www3;

    vec4d x0 = vec4d_load(X+4*0);
    vec4d x1 = vec4d_load(X+4*1);
    vec4d x2 = vec4d_load(X+4*2);
    vec4d x3 = vec4d_load(X+4*3);
    vec4d x4 = vec4d_load(X+4*4);
    vec4d x5 = vec4d_load(X+4*5);
    vec4d x6 = vec4d_load(X+4*6);
    vec4d x7 = vec4d_load(X+4*7);

    ww1  = vec4d_set_d(Q->w2tab[1][0]);
    www2 = vec4d_set_d(Q->w2tab[2][0]);
    www3 = vec4d_set_d(Q->w2tab[2][1]);
    LENGTH8_ZERO_J(vec4d, x0,x1,x2,x3,x4,x5,x6,x7, n,ninv, ww1, www2,www3);

    VEC4D_TRANSPOSE(x0,x1,x2,x3, x0,x1,x2,x3);
    VEC4D_TRANSPOSE(x4,x5,x6,x7, x4,x5,x6,x7);

    /* j = 0, 1, 2, 3 */
    w0  = vec4d_set_d4(Q->w2tab[0][0], Q->w2tab[0+(1)][1*(0)+0], Q->w2tab[0+(2)][1*(0)+0], Q->w2tab[0+(2)][1*(1)+0]);
    ww0 = vec4d_set_d4(Q->w2tab[0][0], Q->w2tab[1+(1)][2*(0)+0], Q->w2tab[1+(2)][2*(0)+0], Q->w2tab[1+(2)][2*(1)+0]);
    ww1 = vec4d_set_d4(Q->w2tab[1][0], Q->w2tab[1+(1)][2*(0)+1], Q->w2tab[1+(2)][2*(0)+1], Q->w2tab[1+(2)][2*(1)+1]);
    LENGTH4_ANY_J(vec4d, x0,x1,x2,x3, n,ninv, w0,ww0,ww1);

    /* j = 4, 5, 6, 7 */
    w0 = vec4d_set_d4(Q->w2tab[0+(3)][1*(0)+0], Q->w2tab[0+(3)][1*(1)+0], Q->w2tab[0+(3)][1*(2)+0], Q->w2tab[0+(3)][1*(3)+0]);
    u  = vec4d_load_aligned(&Q->w2tab[1+(3)][2*(0)+0]);
    v  = vec4d_load_aligned(&Q->w2tab[1+(3)][2*(0)+4]);
    ww0 = vec4d_unpack_lo_permute_0_2_1_3(u, v);
    ww1 = vec4d_unpack_hi_permute_0_2_1_3(u, v);
    LENGTH4_ANY_J(vec4d, x4,x5,x6,x7, n,ninv, w0,ww0,ww1);

    /* VEC4D_TRANSPOSE(x0,x1,x2,x3, x0,x1,x2,x3); */ /* skipped */
    /* VEC4D_TRANSPOSE(x4,x5,x6,x7, x4,x5,x6,x7); */
    vec4d_store(X+4*0, x0);
    vec4d_store(X+4*1, x1);
    vec4d_store(X+4*2, x2);
    vec4d_store(X+4*3, x3);
    vec4d_store(X+4*4, x4);
    vec4d_store(X+4*5, x5);
    vec4d_store(X+4*6, x6);
    vec4d_store(X+4*7, x7);
}

static void sd_fft_basecase_5_0(const sd_fft_ctx_t Q, double* X, ulong j_r, ulong j_bits)
{
    vec4d n    = vec4d_set_d(Q->p);
    vec4d ninv = vec4d_set_d(Q->pinv);
    vec4d u, v, w0, ww0, ww1, www0, www1, www2, www3;

    vec4d x0 = vec4d_load(X+4*0);
    vec4d x1 = vec4d_load(X+4*1);
    vec4d x2 = vec4d_load(X+4*2);
    vec4d x3 = vec4d_load(X+4*3);
    vec4d x4 = vec4d_load(X+4*4);
    vec4d x5 = vec4d_load(X+4*5);
    vec4d x6 = vec4d_load(X+4*6);
    vec4d x7 = vec4d_load(X+4*7);

    w0   = vec4d_set_d(Q->w2tab[0+j_bits][1*j_r+0]);
    ww0  = vec4d_set_d(Q->w2tab[1+j_bits][2*j_r+0]);
    ww1  = vec4d_set_d(Q->w2tab[1+j_bits][2*j_r+1]);
    www0 = vec4d_set_d(Q->w2tab[2+j_bits][4*j_r+0]);
    www1 = vec4d_set_d(Q->w2tab[2+j_bits][4*j_r+1]);
    www2 = vec4d_set_d(Q->w2tab[2+j_bits][4*j_r+2]);
    www3 = vec4d_set_d(Q->w2tab[2+j_bits][4*j_r+3]);
    LENGTH8_ANY_J(vec4d, x0,x1,x2,x3,x4,x5,x6,x7, n,ninv, w0, ww0,ww1, www0,www1,www2,www3);

    VEC4D_TRANSPOSE(x0,x1,x2,x3, x0,x1,x2,x3);
    VEC4D_TRANSPOSE(x4,x5,x6,x7, x4,x5,x6,x7);

    /* j = 8*j+0, 8*j+1, 8*j+2, 8*j+3 */
    w0 = vec4d_load_aligned(&Q->w2tab[0+(3+j_bits)][1*(8*j_r+0)+0]);
    u  = vec4d_load_aligned(&Q->w2tab[1+(3+j_bits)][2*(8*j_r+0)+0]);
    v  = vec4d_load_aligned(&Q->w2tab[1+(3+j_bits)][2*(8*j_r+0)+4]);
    ww0 = vec4d_unpack_lo_permute_0_2_1_3(u, v);
    ww1 = vec4d_unpack_hi_permute_0_2_1_3(u, v);
    LENGTH4_ANY_J(vec4d, x0,x1,x2,x3,n,ninv, w0, ww0, ww1);

    /* j = 8*j+4, 8*j+5, 8*j+6, 8*j+7 */
    w0 = vec4d_load_aligned(&Q->w2tab[0+(3+j_bits)][1*(8*j_r+4)+0]);
    u  = vec4d_load_aligned(&Q->w2tab[1+(3+j_bits)][2*(8*j_r+4)+0]);
    v  = vec4d_load_aligned(&Q->w2tab[1+(3+j_bits)][2*(8*j_r+4)+4]);
    ww0 = vec4d_unpack_lo_permute_0_2_1_3(u, v);
    ww1 = vec4d_unpack_hi_permute_0_2_1_3(u, v);
    LENGTH4_ANY_J(vec4d, x4,x5,x6,x7,n,ninv, w0, ww0, ww1);

    /* VEC4D_TRANSPOSE(x0,x1,x2,x3, x0,x1,x2,x3); */ /* skipped */
    /* VEC4D_TRANSPOSE(x4,x5,x6,x7, x4,x5,x6,x7); */
    vec4d_store(X+4*0, x0);
    vec4d_store(X+4*1, x1);
    vec4d_store(X+4*2, x2);
    vec4d_store(X+4*3, x3);
    vec4d_store(X+4*4, x4);
    vec4d_store(X+4*5, x5);
    vec4d_store(X+4*6, x6);
    vec4d_store(X+4*7, x7);
}


/* use with n = m-2 and m >= 6 */
#define EXTEND_BASECASE(n, m) \
static void CAT3(sd_fft_basecase, m, 1)(const sd_fft_ctx_t Q, double* X) \
{ \
    ulong l = n_pow2(m - 2); \
    RADIX_4_FORWARD_PARAM_J_IS_Z(VECND, Q) \
    ulong i = 0; do { \
        RADIX_4_FORWARD_MOTH_J_IS_Z(VECND, X+0*l+i, X+1*l+i, X+2*l+i, X+3*l+i) \
    } while (i += N, i < l); \
    CAT3(sd_fft_basecase, n, 1)(Q, X+0*l); \
    CAT3(sd_fft_basecase, n, 0)(Q, X+1*l, 0, 1); \
    CAT3(sd_fft_basecase, n, 0)(Q, X+2*l, 0, 2); \
    CAT3(sd_fft_basecase, n, 0)(Q, X+3*l, 1, 2); \
} \
static void CAT3(sd_fft_basecase, m, 0)(const sd_fft_ctx_t Q, double* X, ulong j_r, ulong j_bits) \
{ \
    ulong l = n_pow2(m - 2); \
    RADIX_4_FORWARD_PARAM_J_IS_NZ(VECND, Q, j_r, j_bits) \
    ulong i = 0; do { \
        RADIX_4_FORWARD_MOTH_J_IS_NZ(VECND, X+0*l+i, X+1*l+i, X+2*l+i, X+3*l+i) \
    } while (i += N, i < l); \
    CAT3(sd_fft_basecase, n, 0)(Q, X+0*l, 4*j_r+0, j_bits+2); \
    CAT3(sd_fft_basecase, n, 0)(Q, X+1*l, 4*j_r+1, j_bits+2); \
    CAT3(sd_fft_basecase, n, 0)(Q, X+2*l, 4*j_r+2, j_bits+2); \
    CAT3(sd_fft_basecase, n, 0)(Q, X+3*l, 4*j_r+3, j_bits+2); \
}
#if defined(__AVX512F__)
/* For input x0 = {x00,x01,x02,x03,x04,x05,x06,x07} through
   x7 = {x70,x71,x72,x73,x74,x75,x76,x77}, store
   X[0..63] = {x00,x40,x01,x41,x10,x50,x11,x51,
               x20,x60,x21,x61,x30,x70,x31,x71,
               x02,x42,x03,x43,x12,x52,x13,x53,
               x22,x62,x23,x63,x32,x72,x33,x73,
               x04,x44,x05,x45,x14,x54,x15,x55,
               x24,x64,x25,x65,x34,x74,x35,x75,
               x06,x46,x07,x47,x16,x56,x17,x57,
               x26,x66,x27,x67,x36,x76,x37,x77}. */
FLINT_FORCE_INLINE void
sd_fft_store_vec8dz_basecase_order(double* X, vec8dz x0, vec8dz x1,
                                   vec8dz x2, vec8dz x3, vec8dz x4,
                                   vec8dz x5, vec8dz x6, vec8dz x7)
{
    const vec8nz i01 = vec8nz_set_n8(0, 8, 1, 9, 0, 8, 1, 9);
    const vec8nz i23 = vec8nz_set_n8(2, 10, 3, 11, 2, 10, 3, 11);
    const vec8nz i45 = vec8nz_set_n8(4, 12, 5, 13, 4, 12, 5, 13);
    const vec8nz i67 = vec8nz_set_n8(6, 14, 7, 15, 6, 14, 7, 15);
    vec8dz a01 = vec8dz_permute2var(x0, i01, x4);
    vec8dz b01 = vec8dz_permute2var(x1, i01, x5);
    vec8dz c01 = vec8dz_permute2var(x2, i01, x6);
    vec8dz d01 = vec8dz_permute2var(x3, i01, x7);
    vec8dz a23 = vec8dz_permute2var(x0, i23, x4);
    vec8dz b23 = vec8dz_permute2var(x1, i23, x5);
    vec8dz c23 = vec8dz_permute2var(x2, i23, x6);
    vec8dz d23 = vec8dz_permute2var(x3, i23, x7);
    vec8dz a45 = vec8dz_permute2var(x0, i45, x4);
    vec8dz b45 = vec8dz_permute2var(x1, i45, x5);
    vec8dz c45 = vec8dz_permute2var(x2, i45, x6);
    vec8dz d45 = vec8dz_permute2var(x3, i45, x7);
    vec8dz a67 = vec8dz_permute2var(x0, i67, x4);
    vec8dz b67 = vec8dz_permute2var(x1, i67, x5);
    vec8dz c67 = vec8dz_permute2var(x2, i67, x6);
    vec8dz d67 = vec8dz_permute2var(x3, i67, x7);

    vec8dz_store(X +  0, vec8dz_permute2_128_0_1_4_5(a01, b01));
    vec8dz_store(X +  8, vec8dz_permute2_128_0_1_4_5(c01, d01));
    vec8dz_store(X + 16, vec8dz_permute2_128_0_1_4_5(a23, b23));
    vec8dz_store(X + 24, vec8dz_permute2_128_0_1_4_5(c23, d23));
    vec8dz_store(X + 32, vec8dz_permute2_128_0_1_4_5(a45, b45));
    vec8dz_store(X + 40, vec8dz_permute2_128_0_1_4_5(c45, d45));
    vec8dz_store(X + 48, vec8dz_permute2_128_0_1_4_5(a67, b67));
    vec8dz_store(X + 56, vec8dz_permute2_128_0_1_4_5(c67, d67));
}

static void sd_fft_basecase_6_1(const sd_fft_ctx_t Q, double* X)
{
    vec8dz n    = vec8dz_set_d(Q->p);
    vec8dz ninv = vec8dz_set_d(Q->pinv);
    vec8dz w0, ww0, ww1, www0, www1, www2, www3;

    vec8dz x0 = vec8dz_load(X+8*0);
    vec8dz x1 = vec8dz_load(X+8*1);
    vec8dz x2 = vec8dz_load(X+8*2);
    vec8dz x3 = vec8dz_load(X+8*3);
    vec8dz x4 = vec8dz_load(X+8*4);
    vec8dz x5 = vec8dz_load(X+8*5);
    vec8dz x6 = vec8dz_load(X+8*6);
    vec8dz x7 = vec8dz_load(X+8*7);

    ww1  = vec8dz_set_d(Q->w2tab[1][0]);
    www2 = vec8dz_set_d(Q->w2tab[2][0]);
    www3 = vec8dz_set_d(Q->w2tab[2][1]);
    LENGTH8_ZERO_J(vec8dz, x0,x1,x2,x3,x4,x5,x6,x7, n,ninv, ww1, www2,www3);

    vec8dz_transpose(&x0,&x1,&x2,&x3,&x4,&x5,&x6,&x7, x0,x1,x2,x3,x4,x5,x6,x7);

    w0 = vec8dz_load_aligned(&Q->w2tab[0][0]);
    vec8dz_load_deinterleave2_aligned(&Q->w2tab[0][0], &ww0, &ww1);
    vec8dz_load_deinterleave4_aligned(&Q->w2tab[0][0], &www0, &www1, &www2, &www3);
    LENGTH8_ANY_J(vec8dz, x0,x1,x2,x3,x4,x5,x6,x7, n,ninv, w0,ww0,ww1, www0,www1,www2,www3);

    sd_fft_store_vec8dz_basecase_order(X, x0,x1,x2,x3,x4,x5,x6,x7);
}

static void sd_fft_basecase_6_0(const sd_fft_ctx_t Q, double* X, ulong j_r, ulong j_bits)
{
    vec8dz n    = vec8dz_set_d(Q->p);
    vec8dz ninv = vec8dz_set_d(Q->pinv);
    vec8dz w0, ww0, ww1, www0, www1, www2, www3;
    ulong j = 8*j_r;

    vec8dz x0 = vec8dz_load(X+8*0);
    vec8dz x1 = vec8dz_load(X+8*1);
    vec8dz x2 = vec8dz_load(X+8*2);
    vec8dz x3 = vec8dz_load(X+8*3);
    vec8dz x4 = vec8dz_load(X+8*4);
    vec8dz x5 = vec8dz_load(X+8*5);
    vec8dz x6 = vec8dz_load(X+8*6);
    vec8dz x7 = vec8dz_load(X+8*7);

    w0   = vec8dz_set_d(Q->w2tab[0+j_bits][1*j_r+0]);
    ww0  = vec8dz_set_d(Q->w2tab[1+j_bits][2*j_r+0]);
    ww1  = vec8dz_set_d(Q->w2tab[1+j_bits][2*j_r+1]);
    www0 = vec8dz_set_d(Q->w2tab[2+j_bits][4*j_r+0]);
    www1 = vec8dz_set_d(Q->w2tab[2+j_bits][4*j_r+1]);
    www2 = vec8dz_set_d(Q->w2tab[2+j_bits][4*j_r+2]);
    www3 = vec8dz_set_d(Q->w2tab[2+j_bits][4*j_r+3]);
    LENGTH8_ANY_J(vec8dz, x0,x1,x2,x3,x4,x5,x6,x7, n,ninv, w0,ww0,ww1, www0,www1,www2,www3);

    vec8dz_transpose(&x0,&x1,&x2,&x3,&x4,&x5,&x6,&x7, x0,x1,x2,x3,x4,x5,x6,x7);

    w0 = vec8dz_load_aligned(&Q->w2tab[0+(3+j_bits)][j]);
    vec8dz_load_deinterleave2_aligned(&Q->w2tab[1+(3+j_bits)][2*j], &ww0, &ww1);
    vec8dz_load_deinterleave4_aligned(&Q->w2tab[2+(3+j_bits)][4*j], &www0, &www1, &www2, &www3);
    LENGTH8_ANY_J(vec8dz, x0,x1,x2,x3,x4,x5,x6,x7, n,ninv, w0,ww0,ww1, www0,www1,www2,www3);

    sd_fft_store_vec8dz_basecase_order(X, x0,x1,x2,x3,x4,x5,x6,x7);
}
#else
EXTEND_BASECASE(4, 6)
#endif
EXTEND_BASECASE(5, 7)
EXTEND_BASECASE(6, 8)
EXTEND_BASECASE(7, 9)
#undef EXTEND_BASECASE

/* The `sd_fft_base_{m}_*` functions take `j`, unlike `sd_fft_basecase_{m}_*`
   which takes `j_r` and `j_bits` (or nothing, if `j` is zero).  */

/* parameter 1: j can be zero */
static void sd_fft_base_8_1(const sd_fft_ctx_t Q, double* x, ulong j)
{
    ulong j_bits, j_r;

    FLINT_ASSERT(8 == LG_BLK_SZ);

    SET_J_BITS_AND_J_R(j_bits, j_r, j);

    if (FLINT_UNLIKELY(j == 0))
        sd_fft_basecase_8_1(Q, x);
    else
        sd_fft_basecase_8_0(Q, x, j_r, j_bits);
}

/* parameter 0: j cannot be zero */
static void sd_fft_base_8_0(const sd_fft_ctx_t Q, double* x, ulong j)
{
    ulong j_bits, j_r;

    FLINT_ASSERT(j != 0);
    FLINT_ASSERT(8 == LG_BLK_SZ);

    SET_J_BITS_AND_J_R(j_bits, j_r, j);

    sd_fft_basecase_8_0(Q, x, j_r, j_bits);
}

static void sd_fft_base_9_1(const sd_fft_ctx_t Q, double* x, ulong j)
{
    ulong j_bits, j_r;

    FLINT_ASSERT(8 == LG_BLK_SZ);

    SET_J_BITS_AND_J_R(j_bits, j_r, j);

    if (FLINT_UNLIKELY(j == 0))
        sd_fft_basecase_9_1(Q, x);
    else
        sd_fft_basecase_9_0(Q, x, j_r, j_bits);
}


/**************** forward butterfly with truncation **************************/

/*
    Let `D = X1-X0`, measured in `double` entries. Assume `D = X2-X1 = X3-X2` and `D >= BLK_SZ`.
    `sd_fft_moth_trunc_block_{itrunc}_{otrunc}_{j_is_zero}` computes a possibly
    truncated 2-layer transform with lengths (4*D, 2*D) and mask (0..BLK_SZ).
    Only the first `itrunc` blocks are read (`2 <= itrunc <= 4`),
    missing input blocks are treated as zero, and only the first `otrunc` output
    blocks are written (`1 <= otrunc <= 4`).
*/

/* third parameter is j == 0 */
#define DEFINE_IT(itrunc, otrunc) \
static void CAT4(sd_fft_moth_trunc_block, itrunc, otrunc, 1)( \
    const sd_fft_ctx_t Q, \
    ulong FLINT_UNUSED(j_r), ulong FLINT_UNUSED(j_bits), \
    double* X0, double* X1, double* X2, double* X3) \
{ \
    RADIX_4_FORWARD_PARAM_J_IS_Z(VECND, Q); \
    ulong i = 0; do { \
        VECND x0, x1, x2, x3, y0, y1, y2, y3; \
        x0 = x1 = x2 = x3 = VECNOP(zero)(); \
        if (0 < itrunc) x0 = VECNOP(load)(X0+i); \
        if (0 < itrunc) x0 = VECNOP(reduce_to_pm1n)(x0, n, ninv); \
        if (1 < itrunc) x1 = VECNOP(load)(X1+i); \
        if (2 < itrunc) x2 = VECNOP(load)(X2+i); \
        if (2 < itrunc) x2 = VECNOP(reduce_to_pm1n)(x2, n, ninv); \
        if (3 < itrunc) x3 = VECNOP(load)(X3+i); \
        if (3 < itrunc) x3 = VECNOP(reduce_to_pm1n)(x3, n, ninv); \
        y0 = (2 < itrunc) ? VECNOP(add)(x0, x2) : x0; \
        y1 = (3 < itrunc) ? VECNOP(add)(x1, x3) : x1; \
        y2 = (2 < itrunc) ? VECNOP(sub)(x0, x2) : x0; \
        y3 = (3 < itrunc) ? VECNOP(sub)(x1, x3) : x1; \
        y1 = VECNOP(reduce_to_pm1n)(y1, n, ninv); \
        y3 = VECNOP(mulmod)(y3, iw, n, ninv); \
        x0 = VECNOP(add)(y0, y1); \
        x1 = VECNOP(sub)(y0, y1); \
        x2 = VECNOP(add)(y2, y3); \
        x3 = VECNOP(sub)(y2, y3); \
        if (0 < otrunc) VECNOP(store)(X0+i, x0); \
        if (1 < otrunc) VECNOP(store)(X1+i, x1); \
        if (2 < otrunc) VECNOP(store)(X2+i, x2); \
        if (3 < otrunc) VECNOP(store)(X3+i, x3); \
    } while (i += N, i < BLK_SZ); \
    FLINT_ASSERT(i == BLK_SZ); \
} \
static void CAT4(sd_fft_moth_trunc_block, itrunc, otrunc, 0)( \
    const sd_fft_ctx_t Q, \
    ulong j_r, ulong j_bits, \
    double* X0, double* X1, double* X2, double* X3) \
{ \
    RADIX_4_FORWARD_PARAM_J_IS_NZ(VECND, Q, j_r, j_bits); \
    ulong i = 0; do { \
        VECND x0, x1, x2, x3, y0, y1, y2, y3; \
        x0 = x1 = x2 = x3 = VECNOP(zero)(); \
        if (0 < itrunc) x0 = VECNOP(load)(X0+i); \
        if (0 < itrunc) x0 = VECNOP(reduce_to_pm1n)(x0, n, ninv); \
        if (1 < itrunc) x1 = VECNOP(load)(X1+i); \
        if (2 < itrunc) x2 = VECNOP(load)(X2+i); \
        if (2 < itrunc) x2 = VECNOP(mulmod)(x2, w2, n, ninv); \
        if (3 < itrunc) x3 = VECNOP(load)(X3+i); \
        if (3 < itrunc) x3 = VECNOP(mulmod)(x3, w2, n, ninv); \
        y0 = (2 < itrunc) ? VECNOP(add)(x0, x2) : x0; \
        y1 = (3 < itrunc) ? VECNOP(add)(x1, x3) : x1; \
        y2 = (2 < itrunc) ? VECNOP(sub)(x0, x2) : x0; \
        y3 = (3 < itrunc) ? VECNOP(sub)(x1, x3) : x1; \
        y1 = VECNOP(mulmod)(y1, w, n, ninv); \
        y3 = VECNOP(mulmod)(y3, iw, n, ninv); \
        x0 = VECNOP(add)(y0, y1); \
        x1 = VECNOP(sub)(y0, y1); \
        x2 = VECNOP(add)(y2, y3); \
        x3 = VECNOP(sub)(y2, y3); \
        if (0 < otrunc) VECNOP(store)(X0+i, x0); \
        if (1 < otrunc) VECNOP(store)(X1+i, x1); \
        if (2 < otrunc) VECNOP(store)(X2+i, x2); \
        if (3 < otrunc) VECNOP(store)(X3+i, x3); \
    } while (i += N, i < BLK_SZ); \
    FLINT_ASSERT(i == BLK_SZ); \
}

DEFINE_IT(2, 1)
DEFINE_IT(2, 2)
DEFINE_IT(2, 3)
DEFINE_IT(2, 4)
DEFINE_IT(3, 1)
DEFINE_IT(3, 2)
DEFINE_IT(3, 3)
DEFINE_IT(3, 4)
DEFINE_IT(4, 1)
DEFINE_IT(4, 2)
DEFINE_IT(4, 3)
DEFINE_IT(4, 4)
#undef DEFINE_IT

/************************ the recursive stuff ********************************/

/*
    Compute an untruncated k-layer transform with lengths
    (BLK_SZ*S*2^k, BLK_SZ*S*2^(k-1), ..., BLK_SZ*S*2) and mask (0..BLK_SZ).
*/
static void sd_fft_no_trunc_block(
    const sd_fft_ctx_t Q,
    double* x,
    ulong S, /* stride */
    ulong k, /* BLK_SZ transforms each of length 2^k */
    ulong j)
{
    ulong j_bits, j_r;

    if (k > 4)
    {
        ulong k1 = k/2;
        ulong k2 = k - k1;

        ulong l2 = n_pow2(k2);
        ulong a = 0; do {
            sd_fft_no_trunc_block(Q, x + BLK_SZ*(a*S), S<<k2, k1, j);
        } while (a++, a < l2);

        /* row ffts */
        ulong l1 = n_pow2(k1);
        ulong b = 0; do {
            sd_fft_no_trunc_block(Q, x + BLK_SZ*((b<<k2)*S), S, k2, (j<<k1) + b);
        } while (b++, b < l1);

        return;
    }

    SET_J_BITS_AND_J_R(j_bits, j_r, j);

    if (k >= 2)
    {
        ulong k1 = 2;
        ulong k2 = k - k1;
        ulong l2 = n_pow2(k2);

        /* column ffts */
        if (FLINT_UNLIKELY(j_bits == 0))
        {
            RADIX_4_FORWARD_PARAM_J_IS_Z(VECND, Q)
            ulong a = 0; do {
                double* X0 = x + BLK_SZ*(a*S + (S<<k2)*0);
                double* X1 = x + BLK_SZ*(a*S + (S<<k2)*1);
                double* X2 = x + BLK_SZ*(a*S + (S<<k2)*2);
                double* X3 = x + BLK_SZ*(a*S + (S<<k2)*3);
                ulong i = 0; do {
                    RADIX_4_FORWARD_MOTH_J_IS_Z(VECND, X0+i, X1+i, X2+i, X3+i);
                } while (i += N, i < BLK_SZ);
            } while (a++, a < l2);
        }
        else
        {
            RADIX_4_FORWARD_PARAM_J_IS_NZ(VECND, Q, j_r, j_bits)
            ulong a = 0; do {
                double* X0 = x + BLK_SZ*(a*S + (S<<k2)*0);
                double* X1 = x + BLK_SZ*(a*S + (S<<k2)*1);
                double* X2 = x + BLK_SZ*(a*S + (S<<k2)*2);
                double* X3 = x + BLK_SZ*(a*S + (S<<k2)*3);
                ulong i = 0; do {
                    RADIX_4_FORWARD_MOTH_J_IS_NZ(VECND, X0+i, X1+i, X2+i, X3+i);
                } while (i += N, i < BLK_SZ);
            } while (a++, a < l2);
        }

        if (l2 == 1)
            return;

        /* row ffts */
        ulong l1 = n_pow2(k1);
        ulong b = 0; do {
            sd_fft_no_trunc_block(Q, x + BLK_SZ*((b<<k2)*S), S, k2, (j<<k1) + b);
        } while (b++, b < l1);
    }
    else if (k == 1)
    {
        double* X0 = x + BLK_SZ*(S*0);
        double* X1 = x + BLK_SZ*(S*1);
        if (FLINT_UNLIKELY(j_bits == 0))
        {
            RADIX_2_FORWARD_PARAM_J_IS_Z(VECND, Q)
            ulong i = 0; do {
                RADIX_2_FORWARD_MOTH_J_IS_Z(VECND, X0+i, X1+i);
            } while (i += N, i < BLK_SZ);
        }
        else
        {
            RADIX_2_FORWARD_PARAM_J_IS_NZ(VECND, Q, j_r, j_bits)
            ulong i = 0; do {
                RADIX_2_FORWARD_MOTH_J_IS_NZ(VECND, X0+i, X1+i);
            } while (i += N, i < BLK_SZ);
        }
    }
}

/*
    Computes an untruncated (LG_BLK_SZ + k)-layer contiguous transform with
    lengths (BLK_SZ*2^k, BLK_SZ*2^(k-1), ..., BLK_SZ, BLK_SZ/2, ..., 2).
*/
static void sd_fft_no_trunc_internal(
    const sd_fft_ctx_t Q,
    double* x,
    ulong k,    /* 1 transform of length BLK_SZ*2^k */
    ulong j)    /* twist param */
{
    if (k > 2)
    {
        ulong k1 = k/2;
        ulong k2 = k - k1;

        /* column ffts */
        ulong l2 = n_pow2(k2);
        ulong a = 0; do {
            sd_fft_no_trunc_block(Q, x + BLK_SZ*a, n_pow2(k2), k1, j);
        } while (a++, a < l2);

        /* row ffts */
        ulong l1 = n_pow2(k1);
        ulong b = 0; do {
            sd_fft_no_trunc_internal(Q, x + BLK_SZ*(b<<k2), k2, (j<<k1) + b);
        } while (b++, b < l1);

        return;
    }

    if (k == 2)
    {
        /* k1 = 2; k2 = 0 */
        sd_fft_no_trunc_block(Q, x, 1, 2, j);
        sd_fft_base_8_1(Q, x + BLK_SZ*0, 4*j + 0);
        sd_fft_base_8_0(Q, x + BLK_SZ*1, 4*j + 1);
        sd_fft_base_8_0(Q, x + BLK_SZ*2, 4*j + 2);
        sd_fft_base_8_0(Q, x + BLK_SZ*3, 4*j + 3);
    }
    else if (k == 1)
    {
        sd_fft_base_9_1(Q, x, j);
    }
    else
    {
        /* currently unreachable because all ffts are called with k > 0 */
        sd_fft_base_8_1(Q, x, j);
    }
}


static void sd_fft_trunc_block(
    const sd_fft_ctx_t Q,
    double* x,
    ulong S,
    ulong k, /* transform length 2^k */
    ulong j,
    ulong itrunc,
    ulong otrunc)
{
    ulong j_bits, j_r;

    FLINT_ASSERT(itrunc <= n_pow2(k));
    FLINT_ASSERT(otrunc <= n_pow2(k));

    if (otrunc < 1)
        return;

    if (itrunc <= 1)
    {
        if (itrunc < 1)
        {
            for (ulong a = 0; a < otrunc; a++)
            {
                double* X0 = x + BLK_SZ*(S*a);
                VECND z = VECNOP(zero)();
                ulong i = 0; do {
                    VECNOP(store)(X0+i, z);
                } while (i += N, i < BLK_SZ);
            }
        }
        else
        {
            double* X0 = x + BLK_SZ*(S*0);
            for (ulong a = 1; a < otrunc; a++)
            {
                double* X1 = x + BLK_SZ*(S*a);
                ulong i = 0; do {
                    VECND u = VECNOP(load)(X0+i);
                    VECNOP(store)(X1+i, u);
                } while (i += N, i < BLK_SZ);
            }
        }

        return;
    }

    if (itrunc == otrunc && otrunc == n_pow2(k))
    {
        sd_fft_no_trunc_block(Q, x, S, k, j);
        return;
    }

    if (k > 2)
    {
        ulong k1 = k/2;
        ulong k2 = k - k1;

        ulong l2 = n_pow2(k2);
        ulong n1 = otrunc >> k2;
        ulong n2 = otrunc & (l2 - 1);
        ulong z1 = itrunc >> k2;
        ulong z2 = itrunc & (l2 - 1);
        ulong n1p = n1 + (n2 != 0);
        ulong z2p = n_min(l2, itrunc);

        /* columns */
        for (ulong a = 0; a < z2p; a++)
            sd_fft_trunc_block(Q, x + BLK_SZ*(a*S), S << k2, k1, j, z1 + (a < z2), n1p);

        /* full rows */
        for (ulong b = 0; b < n1; b++)
            sd_fft_trunc_block(Q, x + BLK_SZ*(b*(S << k2)), S, k2, (j << k1) + b, z2p, l2);

        /* last partial row */
        if (n2 > 0)
            sd_fft_trunc_block(Q, x + BLK_SZ*(n1*(S << k2)), S, k2, (j << k1) + n1, z2p, n2);

        return;
    }

    SET_J_BITS_AND_J_R(j_bits, j_r, j);

    if (k == 2)
    {
#define IT(ii, oo) CAT4(sd_fft_moth_trunc_block, ii, oo, 0), \
                   CAT4(sd_fft_moth_trunc_block, ii, oo, 1)
#define LOOKUP_IT(ii, oo, j_is_zero) tab[(j_is_zero) + 2*((oo)-1 + 4*((ii)-2))]
        static void (*tab[3*4*2])(const sd_fft_ctx_t, ulong, ulong, double*, double*, double*, double*) =
                        {IT(2,1), IT(2,2), IT(2,3), IT(2,4),
                         IT(3,1), IT(3,2), IT(3,3), IT(3,4),
                         IT(4,1), IT(4,2), IT(4,3), IT(4,4)};

        double* X0 = x + BLK_SZ*(S*0);
        double* X1 = x + BLK_SZ*(S*1);
        double* X2 = x + BLK_SZ*(S*2);
        double* X3 = x + BLK_SZ*(S*3);
        LOOKUP_IT(itrunc, otrunc, j == 0)(Q, j_r, j_bits, X0, X1, X2, X3);
#undef LOOKUP_IT
#undef IT
    }
    else if (k == 1)
    {
        double* X0 = x + BLK_SZ*(S*0);
        double* X1 = x + BLK_SZ*(S*1);
        FLINT_ASSERT(itrunc == 2);
        FLINT_ASSERT(otrunc == 1);
        if (FLINT_UNLIKELY(j_bits == 0))
        {
            RADIX_2_FORWARD_PARAM_J_IS_Z(VECND, Q)
            ulong i = 0; do {
                RADIX_2_FORWARD_MOTH_TRUNC_2_1_J_IS_Z(VECND, X0 + i, X1 + i);
            } while (i += N, i < BLK_SZ);
        }
        else
        {
            RADIX_2_FORWARD_PARAM_J_IS_NZ(VECND, Q, j_r, j_bits)
            ulong i = 0; do {
                RADIX_2_FORWARD_MOTH_TRUNC_2_1_J_IS_NZ(VECND, X0 + i, X1 + i);
            } while (i += N, i < BLK_SZ);
        }
    }
}


static void sd_fft_trunc_internal(
    const sd_fft_ctx_t Q,
    double* x,      /* x = data + BLK_SZ*I  where I = starting index */
    ulong k,        /* transform length BLK_SZ*2^k */
    ulong j,
    ulong itrunc,   /* actual trunc is BLK_SZ*itrunc */
    ulong otrunc)   /* actual trunc is BLK_SZ*otrunc */
{
    if (otrunc < 1)
        return;

    if (itrunc < 1)
    {
        for (ulong a = 0; a < otrunc; a++)
        {
            double* X0 = x + BLK_SZ*a;
            VECND z = VECNOP(zero)();
            ulong i = 0; do {
                VECNOP(store)(X0 + i, z);
            } while (i += N, i < BLK_SZ);
        }

        return;
    }

    if (itrunc == otrunc && otrunc == n_pow2(k))
    {
        sd_fft_no_trunc_internal(Q, x, k, j);
        return;
    }

    if (k > 2)
    {
        ulong k1 = k/2;
        ulong k2 = k - k1;

        ulong l2 = n_pow2(k2);
        ulong n1 = otrunc >> k2;
        ulong n2 = otrunc & (l2 - 1);
        ulong z1 = itrunc >> k2;
        ulong z2 = itrunc & (l2 - 1);
        ulong n1p = n1 + (n2 != 0);
        ulong z2p = n_min(l2, itrunc);

        /* columns */
        for (ulong a = 0; a < z2p; a++)
            sd_fft_trunc_block(Q, x + BLK_SZ*a, n_pow2(k2), k1, j, z1 + (a < z2), n1p);

        /* full rows */
        for (ulong b = 0; b < n1; b++)
            sd_fft_trunc_internal(Q, x + BLK_SZ*(b << k2), k2, (j << k1) + b, z2p, l2);

        /* last partial row */
        if (n2 > 0)
            sd_fft_trunc_internal(Q, x + BLK_SZ*(n1 << k2), k2, (j << k1) + n1, z2p, n2);

        return;
    }

    if (k == 2)
    {
        sd_fft_trunc_block(Q, x, 1, 2, j, itrunc, otrunc);
                        sd_fft_base_8_1(Q, x + BLK_SZ*0, 4*j + 0);
        if (otrunc > 1) sd_fft_base_8_0(Q, x + BLK_SZ*1, 4*j + 1);
        if (otrunc > 2) sd_fft_base_8_0(Q, x + BLK_SZ*2, 4*j + 2);
        if (otrunc > 3) sd_fft_base_8_0(Q, x + BLK_SZ*3, 4*j + 3);
    }
    else if (k == 1)
    {
        sd_fft_trunc_block(Q, x, 1, 1, j, itrunc, otrunc);
                        sd_fft_base_8_1(Q, x + BLK_SZ*0, 2*j + 0);
        if (otrunc > 1) sd_fft_base_8_0(Q, x + BLK_SZ*1, 2*j + 1);
    }
    else
    {
        /* currently unreachable for the same reason */
        sd_fft_base_8_1(Q, x, j);
    }
}

/********************* interface functions ***********************/

/*
Compute a truncated FFT in place in `d`, assuming all terms after the first `itrunc`
are zero.

The output satisfies
    eval_poly(in_data, sd_fft_ctx_w(Q, i)) = out_data[sd_fft_ctx_trunc_index(L, i)]
for all `0 <= i < otrunc`. This invariant is tested in `test/t-sd_fft.c`.
Usually, it only makes sense to have `otrunc >= itrunc` and `n_max(itrunc, otrunc) >= 2^(L-1)`.
The array `d` needs to have size at least `2^L`.
*/
void sd_fft_trunc(
    sd_fft_ctx_t Q,
    double* d,
    ulong L,    /* convolution length 2^L */
    ulong itrunc, ulong otrunc)
{
    FLINT_ASSERT(itrunc <= n_pow2(L));
    FLINT_ASSERT(otrunc <= n_pow2(L));

    sd_fft_ctx_fit_depth(Q, L);

    if (L > LG_BLK_SZ)
    {
        ulong new_itrunc, new_otrunc;

        new_itrunc = n_cdiv(itrunc, BLK_SZ);
        new_otrunc = n_cdiv(otrunc, BLK_SZ);
        /* this isn't very clever */
        for (int i = 0; i < ((-(int)itrunc)&(BLK_SZ-1)); i++)
            d[itrunc+i] = 0.0;

        sd_fft_trunc_internal(Q, d, L - LG_BLK_SZ, 0, new_itrunc, new_otrunc);
        return;
    }

    /* neither is this */
    for (int i = itrunc; i < (1<<L); i++)
        d[i] = 0;

    /* L=8 reads from w2tab[7] */
    FLINT_ASSERT(LG_BLK_SZ <= SD_FFT_CTX_W2TAB_INIT);

    switch (L) {
        case 0: sd_fft_basecase_0_1(Q, d); break;
        case 1: sd_fft_basecase_1_1(Q, d); break;
        case 2: sd_fft_basecase_2_1(Q, d); break;
        case 3: sd_fft_basecase_3_1(Q, d); break;
        case 4: sd_fft_basecase_4_1(Q, d); break;
        case 5: sd_fft_basecase_5_1(Q, d); break;
        case 6: sd_fft_basecase_6_1(Q, d); break;
        case 7: sd_fft_basecase_7_1(Q, d); break;
        case 8: sd_fft_basecase_8_1(Q, d); break;
        default: FLINT_ASSERT(0);
    }
}
