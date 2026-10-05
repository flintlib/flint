/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fp_contract.h"
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <string.h>
#include "dfloat.h"
#include "gr.h"

/* The error of a double product by Dekker's algorithm (Veltkamp
   splitting), exact when nothing overflows or underflows; a reference
   for fma(a, b, -a b) that does not use fma. */
static double
_dfloat_dekker_err(double a, double b, double p)
{
    double c, ah, al, bh, bl;
    c = 134217729.0 * a;
    ah = c - (c - a);
    al = a - ah;
    c = 134217729.0 * b;
    bh = c - (c - b);
    bl = b - bh;
    return ((ah * bh - p) + ah * bl + al * bh) + al * bl;
}

/* Whether fma() rounds once, as far as a few checks can tell. Some C
   libraries implement it in long double arithmetic or worse (the
   mingw-w64 runtime, as used without -mfma, got 10% of the products
   fma(a, b, -a b) wrong in a test), which breaks the error-free
   products; hardware FMA is always right. The inputs go through
   volatile variables so that the compiler cannot evaluate the calls
   at compile time. */
static int
_dfloat_fma_is_correct(void)
{
    volatile double va, vb, vc;
    double a, b, p, e;
    uint64_t s = UINT64_C(0x9e3779b97f4a7c15), u;
    int i;

    /* double rounding through a 64-bit significand: the exact value
       1 + 2^-52 + 2^-53 - 2^-64 rounds to 1 + 2^-52 */
    va = 0x1.0000000000001p0; vb = 1.0; vc = 0x1.ffcp-54;
    if (fma(va, vb, vc) != 0x1.0000000000001p0)
        return 0;
    /* (1 + 2^-52)^2 - (1 + 2^-51) = 2^-104 */
    va = 0x1.0000000000001p0; vb = 0x1.0000000000001p0; vc = -0x1.0000000000002p0;
    if (fma(va, vb, vc) != 0x1p-104)
        return 0;

    for (i = 0; i < 256; i++)
    {
        s ^= s << 13; s ^= s >> 7; s ^= s << 17;
        /* a significand from s, an exponent in [-64, 63] */
        u = (s & UINT64_C(0x800fffffffffffff)) | ((uint64_t) (1023 - 64 + (int) ((s >> 52) & 127)) << 52);
        memcpy(&a, &u, sizeof(a));
        s ^= s << 13; s ^= s >> 7; s ^= s << 17;
        u = (s & UINT64_C(0x800fffffffffffff)) | ((uint64_t) (1023 - 64 + (int) ((s >> 52) & 127)) << 52);
        memcpy(&b, &u, sizeof(b));
        va = a; vb = b;
        p = a * b;
        vc = -p;
        e = fma(va, vb, vc);
        if (e != _dfloat_dekker_err(a, b, p))
            return 0;
    }

    return 1;
}

/* The error-free transformations need every double operation rounded
   once to double precision (with x87 excess precision, FLT_EVAL_METHOD
   2, they fail; this is evaluated with the flags the module is compiled
   with), no contraction of multiplications and additions into fma (see
   fp_contract.h) and a correctly rounded fma(). */
int
dfloat_is_supported(void)
{
#if !DFLOAT_NO_CONTRACTION
    return 0;
#elif defined(FLT_EVAL_METHOD) && FLT_EVAL_METHOD != 0
    return 0;
#else
    /* 0: not yet known, 1: supported, -1: not supported (the race of
       two threads computing it is harmless) */
    static volatile int supported = 0;

    if (supported == 0)
        supported = (FLT_RADIX == 2 && DBL_MANT_DIG == 53
                        && _dfloat_fma_is_correct()) ? 1 : -1;

    return supported == 1;
#endif
}

extern gr_method_tab_input _d1_gr_methods_input[], _d2_gr_methods_input[],
    _d3_gr_methods_input[], _d4_gr_methods_input[];
extern gr_method_tab_input _d1_gr_strong_methods_input[], _d2_gr_strong_methods_input[],
    _d3_gr_strong_methods_input[], _d4_gr_strong_methods_input[];
extern gr_method_tab_input _d1b_gr_methods_input[], _d2b_gr_methods_input[],
    _d3b_gr_methods_input[], _d4b_gr_methods_input[];
extern gr_method_tab_input _d1b_gr_strong_methods_input[], _d2b_gr_strong_methods_input[],
    _d3b_gr_strong_methods_input[], _d4b_gr_strong_methods_input[];
extern gr_method_tab_input _d1c_gr_methods_input[], _d2c_gr_methods_input[],
    _d3c_gr_methods_input[], _d4c_gr_methods_input[];
extern gr_method_tab_input _d1c_gr_strong_methods_input[], _d2c_gr_strong_methods_input[],
    _d3c_gr_strong_methods_input[], _d4c_gr_strong_methods_input[];
extern gr_method_tab_input _d1cb_gr_methods_input[], _d2cb_gr_methods_input[],
    _d3cb_gr_methods_input[], _d4cb_gr_methods_input[];
extern gr_method_tab_input _d1cb_gr_strong_methods_input[], _d2cb_gr_strong_methods_input[],
    _d3cb_gr_strong_methods_input[], _d4cb_gr_strong_methods_input[];
extern gr_method_tab_input _d1b_gr_fast_methods_input[], _d2b_gr_fast_methods_input[],
    _d3b_gr_fast_methods_input[], _d4b_gr_fast_methods_input[];
extern gr_method_tab_input _d1cb_gr_fast_methods_input[], _d2cb_gr_fast_methods_input[],
    _d3cb_gr_fast_methods_input[], _d4cb_gr_fast_methods_input[];

/* the generic methods of the complex rings (complex.c), extended by the
   formats' own */
extern gr_method_tab_input _dfloat_complex_gr_methods_input[];
extern gr_method_tab_input _dfloat_complex_plain_gr_methods_input[];
extern gr_method_tab_input _dfloat_complex_ball_gr_methods_input[];

/* [kind][n]: kind = ball + 2 complex */
static gr_method_tab_input * const _dfloat_tabs[4][DFLOAT_MAX_N + 1] = {
    {NULL, _d1_gr_methods_input, _d2_gr_methods_input, _d3_gr_methods_input, _d4_gr_methods_input},
    {NULL, _d1b_gr_methods_input, _d2b_gr_methods_input, _d3b_gr_methods_input, _d4b_gr_methods_input},
    {NULL, _d1c_gr_methods_input, _d2c_gr_methods_input, _d3c_gr_methods_input, _d4c_gr_methods_input},
    {NULL, _d1cb_gr_methods_input, _d2cb_gr_methods_input, _d3cb_gr_methods_input, _d4cb_gr_methods_input},
};
static gr_method_tab_input * const _dfloat_strong_tabs[4][DFLOAT_MAX_N + 1] = {
    {NULL, _d1_gr_strong_methods_input, _d2_gr_strong_methods_input, _d3_gr_strong_methods_input, _d4_gr_strong_methods_input},
    {NULL, _d1b_gr_strong_methods_input, _d2b_gr_strong_methods_input, _d3b_gr_strong_methods_input, _d4b_gr_strong_methods_input},
    {NULL, _d1c_gr_strong_methods_input, _d2c_gr_strong_methods_input, _d3c_gr_strong_methods_input, _d4c_gr_strong_methods_input},
    {NULL, _d1cb_gr_strong_methods_input, _d2cb_gr_strong_methods_input, _d3cb_gr_strong_methods_input, _d4cb_gr_strong_methods_input},
};
static gr_method_tab_input * const _dfloat_fast_tabs[4][DFLOAT_MAX_N + 1] = {
    {NULL, NULL, NULL, NULL, NULL},
    {NULL, _d1b_gr_fast_methods_input, _d2b_gr_fast_methods_input, _d3b_gr_fast_methods_input, _d4b_gr_fast_methods_input},
    {NULL, NULL, NULL, NULL, NULL},
    {NULL, _d1cb_gr_fast_methods_input, _d2cb_gr_fast_methods_input, _d3cb_gr_fast_methods_input, _d4cb_gr_fast_methods_input},
};

/* [kind][strong][fast][n] */
static gr_static_method_table _dfloat_methods[4][2][2][DFLOAT_MAX_N + 1];
static int _dfloat_methods_initialized[4][2][2][DFLOAT_MAX_N + 1];

int
gr_ctx_init_dfloat(gr_ctx_t ctx, int n, int flags)
{
    int ball, strong, fast, cplx, kind;

    if (n < 1 || n > DFLOAT_MAX_N || !dfloat_is_supported())
        return GR_UNABLE;

    ball = (flags & DFLOAT_BALL) != 0;
    strong = (flags & DFLOAT_STRONG) != 0;
    fast = ball && (flags & DFLOAT_FAST) != 0;
    cplx = (flags & DFLOAT_COMPLEX) != 0;
    kind = ball + 2 * cplx;

    ctx->which_ring = cplx ? (ball ? GR_CTX_DFLOAT_COMPLEX_BALL : GR_CTX_DFLOAT_COMPLEX)
                           : (ball ? GR_CTX_DFLOAT_BALL : GR_CTX_DFLOAT);
    ctx->sizeof_elem = sizeof(double) * (n + ball) * (1 + cplx);
    ctx->size_limit = WORD_MAX;
    DFLOAT_CTX_N(ctx) = n;
    DFLOAT_CTX_BALL(ctx) = ball;
    DFLOAT_CTX_STRONG(ctx) = strong;
    DFLOAT_CTX_FAST(ctx) = fast;
    DFLOAT_CTX_COMPLEX(ctx) = cplx;

    if (!_dfloat_methods_initialized[kind][strong][fast][n])
    {
        gr_funcptr * methods = _dfloat_methods[kind][strong][fast][n];

        if (cplx)
        {
            gr_method_tab_init(methods, _dfloat_complex_gr_methods_input);
            gr_method_tab_extend(methods, ball ? _dfloat_complex_ball_gr_methods_input
                                               : _dfloat_complex_plain_gr_methods_input);
            gr_method_tab_extend(methods, _dfloat_tabs[kind][n]);
        }
        else
            gr_method_tab_init(methods, _dfloat_tabs[kind][n]);
        if (fast)
            gr_method_tab_extend(methods, _dfloat_fast_tabs[kind][n]);
        /* the strong wrappers consult the fast flag themselves */
        if (strong)
            gr_method_tab_extend(methods, _dfloat_strong_tabs[kind][n]);

        _dfloat_methods_initialized[kind][strong][fast][n] = 1;
    }

    ctx->methods = _dfloat_methods[kind][strong][fast][n];
    return GR_SUCCESS;
}
