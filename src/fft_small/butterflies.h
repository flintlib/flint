/*
    Copyright (C) 2022 Daniel Schultz
    Copyright (C) 2026 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef FFT_SMALL_BUTTERFLIES_H
#define FFT_SMALL_BUTTERFLIES_H

/*
    Butterflies shared by the four transforms sd_fft.c, sd_ifft.c,
    sd_fft_tr.c and sd_ifft_tr.c. They come in two shapes:

    - DIT, the shape of the forward transform of sd_fft.c:
          (a0, a1)  |-->  (a0 + w*a1, a0 - w*a1);
    - DIF, the shape of the inverse transform of sd_ifft.c:
          (a0, a1)  |-->  (a0 + a1, W*(a1 - a0)),
      where sd_ifft.c uses W = -w^-1.

    The transpose of a DIT butterfly with twiddle w is the DIF butterfly
    with W = -w, and the transpose of a DIF butterfly with twiddle W is the
    DIT butterfly with w = -W. So sd_fft_tr.c uses the DIF butterflies, and
    sd_ifft_tr.c the DIT ones, with negated twiddles; the radix-4
    butterflies and the register basecases below are compositions of these
    and transpose in the same way.

    The negation is folded into the multiplications: every butterfly takes
    the multiplication to use as an argument, MULF (and MUL, where an
    input may exceed the range of mulmod_fast), which is mulmod_fast (and
    mulmod) in sd_fft.c and sd_ifft.c, and nmulmod_fast (and nmulmod) in
    the transposed transforms. Since nmulmod_fast(a, b) is exactly
    -mulmod_fast(a, b), and likewise for nmulmod, the transposed transforms
    compute bit for bit what they would with negated twiddles, and every
    intermediate value has the magnitude it has in the transform whose
    butterflies they share: the bound analyses of the reverse and forward
    butterflies in mulmod_satisfies_bounds.c apply unchanged.

    The twiddles are explicit arguments. The butterflies acting on blocks
    (MOTH, on vectors of type V loaded from and stored to X0, X1, ...) use
    the variables n and ninv of the caller, which its PARAM macros define;
    the register basecases (LENGTH, on variables of type T) take them as
    arguments. See sd_fft.c for how these macros are used.
*/

/************************** DIT butterflies ***********************************
    b0 = a0 + w*a1
    b1 = a0 - w*a1
*/

#define SD_DIT2_MOTH_J_IS_Z(V, X0, X1) \
{ \
    V x0, x1; \
    x0 = CAT(V, load)(X0); \
    x0 = CAT(V, reduce_to_pm1n)(x0, n, ninv); \
    x1 = CAT(V, load)(X1); \
    x1 = CAT(V, reduce_to_pm1n)(x1, n, ninv); \
    CAT(V, store)(X0, CAT(V, add)(x0, x1)); \
    CAT(V, store)(X1, CAT(V, sub)(x0, x1)); \
}

#define SD_DIT2_MOTH_J_IS_NZ(V, MULF, w, X0, X1) \
{ \
    V x0, x1; \
    x0 = CAT(V, load)(X0); \
    x0 = CAT(V, reduce_to_pm1n)(x0, n, ninv); \
    x1 = CAT(V, load)(X1); \
    x1 = CAT(V, MULF)(x1, w, n, ninv); \
    CAT(V, store)(X0, CAT(V, add)(x0, x1)); \
    CAT(V, store)(X1, CAT(V, sub)(x0, x1)); \
}

/*
    b0 = a0 + w^2*a2 +   w*(a1 + w^2*a3)
    b1 = a0 + w^2*a2 -   w*(a1 + w^2*a3)
    b2 = a0 - w^2*a2 + i*w*(a1 - w^2*a3)
    b3 = a0 - w^2*a2 - i*w*(a1 - w^2*a3)
    with (w2, w, iw) = (w^2, w, i*w); for j == 0, w = 1 and iw = i.
*/

#define SD_DIT4_MOTH_J_IS_Z(V, MULF, iw, X0, X1, X2, X3) \
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
    y3 = CAT(V, MULF)(y3, iw, n, ninv); \
    x0 = CAT(V, add)(y0, y1); \
    x1 = CAT(V, sub)(y0, y1); \
    x2 = CAT(V, add)(y2, y3); \
    x3 = CAT(V, sub)(y2, y3); \
    CAT(V, store)(X0, x0); \
    CAT(V, store)(X1, x1); \
    CAT(V, store)(X2, x2); \
    CAT(V, store)(X3, x3); \
}

#define SD_DIT4_MOTH_J_IS_NZ(V, MULF, w2, w, iw, X0, X1, X2, X3) \
{ \
    V x0, x1, x2, x3, y0, y1, y2, y3; \
    x0 = CAT(V, load)(X0); \
    x0 = CAT(V, reduce_to_pm1n)(x0, n, ninv); \
    x1 = CAT(V, load)(X1); \
    x2 = CAT(V, load)(X2); \
    x3 = CAT(V, load)(X3); \
    x2 = CAT(V, MULF)(x2, w2, n, ninv); \
    x3 = CAT(V, MULF)(x3, w2, n, ninv); \
    y0 = CAT(V, add)(x0, x2); \
    y1 = CAT(V, add)(x1, x3); \
    y2 = CAT(V, sub)(x0, x2); \
    y3 = CAT(V, sub)(x1, x3); \
    y1 = CAT(V, MULF)(y1, w, n, ninv); \
    y3 = CAT(V, MULF)(y3, iw, n, ninv); \
    x0 = CAT(V, add)(y0, y1); \
    x1 = CAT(V, sub)(y0, y1); \
    x2 = CAT(V, add)(y2, y3); \
    x3 = CAT(V, sub)(y2, y3); \
    CAT(V, store)(X0, x0); \
    CAT(V, store)(X1, x1); \
    CAT(V, store)(X2, x2); \
    CAT(V, store)(X3, x3); \
}

/* DIT register basecases: 1, 2 and 3 layers */

#define SD_DIT_LENGTH2_ZERO_J(T, x0, x1, n, ninv) \
{ \
    T Z0 = CAT(T, reduce_to_pm1n)(x0, n, ninv); \
    T Z1 = CAT(T, reduce_to_pm1n)(x1, n, ninv); \
    x0 = CAT(T, add)(Z0, Z1); \
    x1 = CAT(T, sub)(Z0, Z1); \
}

#define SD_DIT_LENGTH2_ANY_J(T, MULF, x0, x1, n, ninv, w0) \
{ \
    T Z0 = CAT(T, reduce_to_pm1n)(x0, n, ninv); \
    T Z1 = CAT(T, MULF)(x1, w0, n, ninv); \
    x0 = CAT(T, add)(Z0, Z1); \
    x1 = CAT(T, sub)(Z0, Z1); \
}

#define SD_DIT_LENGTH4_ANY_J(T, MULF, x0, x1, x2, x3, n, ninv, w2, w, iw) \
{ \
    T X0 = x0, X1 = x1, X2 = x2, X3 = x3, Y0, Y1, Y2, Y3; \
    X0 = CAT(T, reduce_to_pm1n)(X0, n, ninv); \
    X2 = CAT(T, MULF)(X2, w2, n, ninv); \
    X3 = CAT(T, MULF)(X3, w2, n, ninv); \
    Y0 = CAT(T, add)(X0, X2); \
    Y1 = CAT(T, add)(X1, X3); \
    Y2 = CAT(T, sub)(X0, X2); \
    Y3 = CAT(T, sub)(X1, X3); \
    Y1 = CAT(T, MULF)(Y1, w, n, ninv); \
    Y3 = CAT(T, MULF)(Y3, iw, n, ninv); \
    x0 = CAT(T, add)(Y0, Y1); \
    x1 = CAT(T, sub)(Y0, Y1); \
    x2 = CAT(T, add)(Y2, Y3); \
    x3 = CAT(T, sub)(Y2, Y3); \
}

#define SD_DIT_LENGTH4_ZERO_J(T, MULF, x0, x1, x2, x3, n, ninv, e14) \
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
    Y3 = CAT(T, MULF)(Y3, e14, n, ninv); \
    x0 = CAT(T, add)(Y0, Y1); \
    x1 = CAT(T, sub)(Y0, Y1); \
    x2 = CAT(T, add)(Y2, Y3); \
    x3 = CAT(T, sub)(Y2, Y3); \
}

#define SD_DIT_LENGTH8_ANY_J(T, MULF, x0, x1, x2, x3, x4, x5, x6, x7, n, ninv, w2, w, iw, ww0, ww1, ww2, ww3) \
{ \
    T X0 = x0, X1 = x1, X2 = x2, X3 = x3, X4 = x4, X5 = x5, X6 = x6, X7 = x7; \
    T Y0, Y1, Y2, Y3, Y4, Y5, Y6, Y7, Z0, Z1, Z2, Z3, Z4, Z5, Z6, Z7; \
    X0 = CAT(T, reduce_to_pm1n)(X0, n, ninv); \
    X1 = CAT(T, reduce_to_pm1n)(X1, n, ninv); \
    X4 = CAT(T, MULF)(X4, w2, n, ninv); \
    X5 = CAT(T, MULF)(X5, w2, n, ninv); \
    X6 = CAT(T, MULF)(X6, w2, n, ninv); \
    X7 = CAT(T, MULF)(X7, w2, n, ninv); \
    Y0 = CAT(T, add)(X0, X4); \
    Y1 = CAT(T, add)(X1, X5); \
    Y2 = CAT(T, add)(X2, X6); \
    Y3 = CAT(T, add)(X3, X7); \
    Y4 = CAT(T, sub)(X0, X4); \
    Y5 = CAT(T, sub)(X1, X5); \
    Y6 = CAT(T, sub)(X2, X6); \
    Y7 = CAT(T, sub)(X3, X7); \
    Y2 = CAT(T, MULF)(Y2, w, n, ninv); \
    Y3 = CAT(T, MULF)(Y3, w, n, ninv); \
    Y6 = CAT(T, MULF)(Y6, iw, n, ninv); \
    Y7 = CAT(T, MULF)(Y7, iw, n, ninv); \
    Z0 = CAT(T, add)(Y0, Y2); \
    Z1 = CAT(T, add)(Y1, Y3); \
    Z2 = CAT(T, sub)(Y0, Y2); \
    Z3 = CAT(T, sub)(Y1, Y3); \
    Z4 = CAT(T, add)(Y4, Y6); \
    Z5 = CAT(T, add)(Y5, Y7); \
    Z6 = CAT(T, sub)(Y4, Y6); \
    Z7 = CAT(T, sub)(Y5, Y7); \
    Z0 = CAT(T, reduce_to_pm1n)(Z0, n, ninv); \
    Z1 = CAT(T, MULF)(Z1, ww0, n, ninv); \
    Z2 = CAT(T, reduce_to_pm1n)(Z2, n, ninv); \
    Z3 = CAT(T, MULF)(Z3, ww1, n, ninv); \
    Z4 = CAT(T, reduce_to_pm1n)(Z4, n, ninv); \
    Z5 = CAT(T, MULF)(Z5, ww2, n, ninv); \
    Z6 = CAT(T, reduce_to_pm1n)(Z6, n, ninv); \
    Z7 = CAT(T, MULF)(Z7, ww3, n, ninv); \
    x0 = CAT(T, add)(Z0, Z1); \
    x1 = CAT(T, sub)(Z0, Z1); \
    x2 = CAT(T, add)(Z2, Z3); \
    x3 = CAT(T, sub)(Z2, Z3); \
    x4 = CAT(T, add)(Z4, Z5); \
    x5 = CAT(T, sub)(Z4, Z5); \
    x6 = CAT(T, add)(Z6, Z7); \
    x7 = CAT(T, sub)(Z6, Z7); \
}

/************************** DIF butterflies ***********************************
    b0 =    (a0 + a1)
    b1 = W*(a1 - a0)
*/

#define SD_DIF2_MOTH_J_IS_Z(V, X0, X1) \
{ \
    V x0, x1, y0, y1; \
    x0 = CAT(V, load)(X0); \
    x1 = CAT(V, load)(X1); \
    y0 = CAT(V, add)(x0, x1); \
    y1 = CAT(V, sub)(x0, x1); \
    y0 = CAT(V, reduce_to_pm1n)(y0, n, ninv); \
    y1 = CAT(V, reduce_to_pm1n)(y1, n, ninv); \
    CAT(V, store)(X0, y0); \
    CAT(V, store)(X1, y1); \
}

#define SD_DIF2_MOTH_J_IS_NZ(V, MULF, W, X0, X1) \
{ \
    V x0, x1, y0, y1; \
    x0 = CAT(V, load)(X0); \
    x1 = CAT(V, load)(X1); \
    y0 = CAT(V, add)(x0, x1); \
    y1 = CAT(V, sub)(x1, x0); \
    y0 = CAT(V, reduce_to_pm1n)(y0, n, ninv); \
    y1 = CAT(V, MULF)(y1, W, n, ninv); \
    CAT(V, store)(X0, y0); \
    CAT(V, store)(X1, y1); \
}

/*
    b0 =      (a0 + a1) + (a2 + a3)
    b1 =      IW*(a3 - a2) - W*(a0 - a1)
    b2 = W2*( (a2 + a3) - (a0 + a1) )
    b3 = W2*( IW*(a3 - a2) + W*(a0 - a1) )
    with, in sd_ifft.c, (W2, W, IW) = (-w^-2, -w^-1, i*w^-1). For j == 0
    only IW is used: b0, b2 = (a0 + a1) +- (a2 + a3) and
    b1, b3 = (a0 - a1) -+ IW*(a2 - a3).
*/

#define SD_DIF4_MOTH_J_IS_Z(V, MULF, IW, X0, X1, X2, X3) \
{ \
    V x0, x1, x2, x3, y0, y1, y2, y3; \
    x0 = CAT(V, load)(X0); \
    x1 = CAT(V, load)(X1); \
    x2 = CAT(V, load)(X2); \
    x3 = CAT(V, load)(X3); \
    y0 = CAT(V, add)(x0, x1); \
    y1 = CAT(V, add)(x2, x3); \
    y2 = CAT(V, sub)(x0, x1); \
    y3 = CAT(V, sub)(x2, x3); \
    y0 = CAT(V, reduce_to_pm1n)(y0, n, ninv); \
    y1 = CAT(V, reduce_to_pm1n)(y1, n, ninv); \
    y2 = CAT(V, reduce_to_pm1n)(y2, n, ninv); \
    y3 = CAT(V, MULF)(y3, IW, n, ninv); \
    CAT(V, store)(X0, CAT(V, add)(y0, y1)); \
    CAT(V, store)(X2, CAT(V, sub)(y0, y1)); \
    CAT(V, store)(X1, CAT(V, sub)(y2, y3)); \
    CAT(V, store)(X3, CAT(V, add)(y2, y3)); \
}

#define SD_DIF4_MOTH_J_IS_NZ(V, MULF, MUL, W2, W, IW, X0, X1, X2, X3) \
{ \
    V x0, x1, x2, x3, y0, y1, y2, y3; \
    x0 = CAT(V, load)(X0); \
    x1 = CAT(V, load)(X1); \
    x2 = CAT(V, load)(X2); \
    x3 = CAT(V, load)(X3); \
    y0 = CAT(V, add)(x0, x1); \
    y1 = CAT(V, add)(x2, x3); \
    y2 = CAT(V, sub)(x0, x1); \
    y3 = CAT(V, sub)(x3, x2); \
    y2 = CAT(V, MULF)(y2, W, n, ninv); \
    y3 = CAT(V, MULF)(y3, IW, n, ninv); \
    x0 = CAT(V, add)(y0, y1); \
    x1 = CAT(V, sub)(y3, y2); \
    CAT(V, store)(X1, x1); \
    x2 = CAT(V, sub)(y1, y0); \
    x3 = CAT(V, add)(y3, y2); \
    x0 = CAT(V, reduce_to_pm1n)(x0, n, ninv); \
    /* |x2| < 8n: may exceed the mulmod_fast range */ \
    x2 = CAT(V, MUL)(x2, W2, n, ninv); \
    x3 = CAT(V, MULF)(x3, W2, n, ninv); \
    CAT(V, store)(X0, x0); \
    CAT(V, store)(X2, x2); \
    CAT(V, store)(X3, x3); \
}

/* DIF register basecases: 1, 2 and 3 layers */

#define SD_DIF_LENGTH2_ZERO_J(T, x0, x1, n, ninv) \
{ \
   T Z0 = x0, Z1 = x1; \
   x0 = CAT(T, reduce_to_pm1n)(CAT(T, add)(Z0, Z1), n, ninv); \
   x1 = CAT(T, reduce_to_pm1n)(CAT(T, sub)(Z0, Z1), n, ninv); \
}

#define SD_DIF_LENGTH2_ANY_J(T, MULF, x0, x1, n, ninv, w0) \
{ \
   T Z0 = x0, Z1 = x1; \
   x0 = CAT(T, reduce_to_pm1n)(CAT(T, add)(Z0, Z1), n, ninv); \
   x1 = CAT(T, MULF)(CAT(T, sub)(Z1, Z0), w0, n, ninv); \
}

#define SD_DIF_LENGTH4_ANY_J(T, MULF, MUL, x0, x1, x2, x3, n, ninv, w0, ww0, ww1) \
{ \
    T X0 = x0, X1 = x1, X2 = x2, X3 = x3, Y0, Y1, Y2, Y3, Z0, Z1, Z2, Z3; \
    Y0 = CAT(T, add)(X0, X1); \
    Y1 = CAT(T, add)(X2, X3); \
    Y2 = CAT(T, sub)(X0, X1); \
    Y3 = CAT(T, sub)(X3, X2); \
    Y2 = CAT(T, MULF)(Y2, ww0, n, ninv); \
    Y3 = CAT(T, MULF)(Y3, ww1, n, ninv); \
    Z0 = CAT(T, add)(Y0, Y1); \
    Z1 = CAT(T, sub)(Y3, Y2); \
    Z2 = CAT(T, sub)(Y1, Y0); \
    Z3 = CAT(T, add)(Y3, Y2); \
    x0 = CAT(T, reduce_to_pm1n)(Z0, n, ninv); \
    x1 = Z1; \
    /* |Z2| < 8n: may exceed the mulmod_fast range */ \
    x2 = CAT(T, MUL)(Z2, w0, n, ninv); \
    x3 = CAT(T, MULF)(Z3, w0, n, ninv); \
}

#define SD_DIF_LENGTH4_ZERO_J(T, MULF, x0, x1, x2, x3, n, ninv, ww1) \
{ \
    T X0 = x0, X1 = x1, X2 = x2, X3 = x3, Y0, Y1, Y2, Y3; \
    Y0 = CAT(T, add)(X0, X1); \
    Y1 = CAT(T, add)(X2, X3); \
    Y2 = CAT(T, sub)(X0, X1); \
    Y3 = CAT(T, sub)(X2, X3); \
    Y0 = CAT(T, reduce_to_pm1n)(Y0, n, ninv); \
    Y1 = CAT(T, reduce_to_pm1n)(Y1, n, ninv); \
    Y2 = CAT(T, reduce_to_pm1n)(Y2, n, ninv); \
    Y3 = CAT(T, MULF)(Y3, ww1, n, ninv); \
    x0 = CAT(T, add)(Y0, Y1); \
    x2 = CAT(T, sub)(Y0, Y1); \
    x1 = CAT(T, sub)(Y2, Y3); \
    x3 = CAT(T, add)(Y2, Y3); \
}

#define SD_DIF_LENGTH8_ANY_J(T, MULF, MUL, x0, x1, x2, x3, x4, x5, x6, x7, n, ninv, w0, ww0, ww1, www0, www1, www2, www3) \
{ \
    T A0 = x0, A1 = x1, A2 = x2, A3 = x3, A4 = x4, A5 = x5, A6 = x6, A7 = x7; \
    SD_DIF_LENGTH2_ANY_J(T, MULF, A0,A1, n,ninv, www0) \
    SD_DIF_LENGTH2_ANY_J(T, MULF, A2,A3, n,ninv, www1) \
    SD_DIF_LENGTH2_ANY_J(T, MULF, A4,A5, n,ninv, www2) \
    SD_DIF_LENGTH2_ANY_J(T, MULF, A6,A7, n,ninv, www3) \
    SD_DIF_LENGTH4_ANY_J(T, MULF, MUL, A0,A2,A4,A6, n,ninv, w0,ww0,ww1) \
    SD_DIF_LENGTH4_ANY_J(T, MULF, MUL, A1,A3,A5,A7, n,ninv, w0,ww0,ww1) \
    x0 = A0; \
    x1 = A1; \
    x2 = A2; \
    x3 = A3; \
    x4 = A4; \
    x5 = A5; \
    x6 = A6; \
    x7 = A7; \
}

#define SD_DIF_LENGTH8_ZERO_J(T, MULF, x0, x1, x2, x3, x4, x5, x6, x7, n, ninv, ww1, www2, www3) \
{ \
    T A0 = x0, A1 = x1, A2 = x2, A3 = x3, A4 = x4, A5 = x5, A6 = x6, A7 = x7; \
    SD_DIF_LENGTH2_ZERO_J(T, A0,A1, n,ninv) \
    SD_DIF_LENGTH2_ANY_J(T, MULF, A2,A3, n,ninv, ww1) \
    SD_DIF_LENGTH2_ANY_J(T, MULF, A4,A5, n,ninv, www2) \
    SD_DIF_LENGTH2_ANY_J(T, MULF, A6,A7, n,ninv, www3) \
    SD_DIF_LENGTH4_ZERO_J(T, MULF, A0,A2,A4,A6, n,ninv, ww1) \
    SD_DIF_LENGTH4_ZERO_J(T, MULF, A1,A3,A5,A7, n,ninv, ww1) \
    x0 = A0; \
    x1 = A1; \
    x2 = A2; \
    x3 = A3; \
    x4 = A4; \
    x5 = A5; \
    x6 = A6; \
    x7 = A7; \
}

#endif
