/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Regression tests for extreme exponents, indices and the exact precision
   (integer overflows that random tests rarely or never reach). */

#include <string.h>
#include "test_helpers.h"
#include "decimal.h"
#include "fmpq.h"
#include "fmpzi.h"
#include "arb.h"
#include "qqbar.h"
#include "gr.h"

static void
_xcheck(int ok, const char * what)
{
    if (!ok)
    {
        flint_printf("FAIL: %s\n", what);
        flint_abort();
    }
}

static void
_xcheck_str(char * s, const char * expected, const char * what)
{
    if (strcmp(s, expected))
    {
        flint_printf("FAIL: %s\ngot %s, expected %s\n", what, s, expected);
        flint_abort();
    }
    flint_free(s);
}

#define XCHECK(expr) _xcheck((expr), #expr)
#define XCHECK_STR(s, t) _xcheck_str((s), (t), #s)

TEST_FUNCTION_START(decimal_extreme, state)
{
    gr_ctx_t F, B, BX, CBX, ZI;
    decfloat_t x, y;
    decball_t X, Y;
    deccball_t Z, W, V;
    decmag_t u, v;
    fmpz_t m, t;
    fmpq_t q;
    arb_t a;
    arf_t r;
    qqbar_t s;
    fmpzi_t g;
    slong e, k;
    int c;

    gr_ctx_init_decfloat(F, 30, 0);
    gr_ctx_init_decball(B, 30, 0);
    gr_ctx_init_decball(BX, DECIMAL_PREC_EXACT, 0);
    gr_ctx_init_deccball(CBX, DECIMAL_PREC_EXACT, 0);
    gr_ctx_init_fmpzi(ZI);

    decfloat_init(x, F);
    decfloat_init(y, F);
    fmpz_init(m);
    fmpz_init(t);
    fmpq_init(q);
    arb_init(a);
    arf_init(r);
    qqbar_init(s);
    fmpzi_init(g);

    e = DECIMAL_CTX_E(F);

    /* conversion of m 2^t with t = WORD_MIN (FLINT_ABS(t) overflowed and
       GMP was asked for 5^(huge)) */
    decball_init(X, B);
    arf_set_ui_2exp_si(arb_midref(a), 3, WORD_MIN);
    mag_zero(arb_radref(a));
    XCHECK(decball_set_arb(X, a, B) == GR_SUCCESS);
    XCHECK(decfloat_set_round_arf(x, arb_midref(a), 30, DECIMAL_RND_NEAR, F) == GR_SUCCESS);
    XCHECK(!DECFLOAT_IS_SPECIAL(x));
    mag_set_ui_2exp_si(arb_radref(a), 1, WORD_MIN / 2);
    XCHECK(decball_set_arb(X, a, B) == GR_SUCCESS);
    decball_clear(X, B);

    /* exact conversion to arf of 3 * 10^(-k e) where k e (about
       11 2^(FLINT_BITS - 4)) wraps around to a negative slong, which was
       used as a table index in the dyadic pre-check (an out-of-bounds read
       far from the table) */
    k = (slong) ((UWORD(11) << (FLINT_BITS - 4)) / (ulong) e);
    fmpz_set_ui(m, 3);
    fmpz_set_si(t, -k);
    fmpz_mul_ui(t, t, e);
    XCHECK(decfloat_set_fmpz_10exp_fmpz(x, m, t, F) == GR_SUCCESS);
    XCHECK(decfloat_get_arf(r, x, ARF_PREC_EXACT, ARF_RND_DOWN, F) == GR_UNABLE);

    /* limb and digit access at extreme indices */
    XCHECK(decfloat_set_si(x, 1, F) == GR_SUCCESS);
    XCHECK(decfloat_set_limb_si(y, x, WORD_MIN, 1, F) == GR_UNABLE);
    XCHECK(decfloat_set_limb_si(y, x, WORD_MAX, 1, F) == GR_UNABLE);
    XCHECK(decfloat_set_digit_si(y, x, WORD_MIN, 1, F) == GR_UNABLE);
    XCHECK(decfloat_set_digit_si(y, x, WORD_MAX, 1, F) == GR_UNABLE);
    XCHECK(decfloat_get_digit_si(x, WORD_MIN, F) == 0);
    XCHECK(decfloat_get_digit_si(x, WORD_MIN + 5, F) == 0);
    XCHECK(decfloat_get_digit_si(x, WORD_MAX, F) == 0);
    XCHECK(decfloat_get_limb_si(x, WORD_MIN, F) == 0);
    XCHECK(decfloat_get_limb_si(x, WORD_MAX, F) == 0);
    XCHECK(decfloat_get_digit_si(x, 0, F) == 1);
    XCHECK(decfloat_set_digit_si(y, x, -3, 7, F) == GR_SUCCESS);
    XCHECK_STR(decfloat_get_str(y, F), "1.007");

    /* decmag exponents near WORD_MIN / WORD_MAX */
    _decmag_init(u, B);
    _decmag_init(v, B);
    _decmag_one(u, B);
    XCHECK(_decmag_cmp_10exp_si(u, WORD_MIN, B) > 0);
    XCHECK(_decmag_cmp_10exp_si(u, WORD_MAX, B) < 0);
    _decmag_mul_10exp_si(v, u, WORD_MAX, B);
    XCHECK(_decmag_cmp_10exp_si(v, WORD_MAX - 1, B) > 0);
    _decmag_mul_10exp_si(v, u, WORD_MIN, B);
    XCHECK(_decmag_cmp_10exp_si(v, WORD_MIN + 1, B) < 0);
    _decmag_clear(u, B);
    _decmag_clear(v, B);

    decball_init(X, B);
    decball_zero(X, B);
    XCHECK(decball_add_error_10exp_si(X, WORD_MIN, B) == GR_SUCCESS);
    XCHECK(_decmag_cmp_10exp_si(&X->rad, -1000, B) < 0);
    decball_clear(X, B);

    /* containment of a rational: huge and tiny radii (10^|exp| must not
       be formed), and moderate radii decided exactly */
    decball_init(X, B);
    fmpq_set_si(q, 1, 3);
    decball_zero(X, B);
    XCHECK(decball_add_error_10exp_si(X, WORD_MAX / 2, B) == GR_SUCCESS);
    XCHECK(_decball_contains_fmpq(X, q, B) == 1);
    XCHECK(decball_set_si(X, 1, B) == GR_SUCCESS);
    XCHECK(decball_add_error_10exp_si(X, -(WORD_MAX / 2), B) == GR_SUCCESS);
    XCHECK(_decball_contains_fmpq(X, q, B) == 0);
    XCHECK(decball_set_str(X, "[0.3 +/- 0.04]", B) == GR_SUCCESS);
    XCHECK(_decball_contains_fmpq(X, q, B) == 1);
    XCHECK(decball_set_str(X, "[0.3 +/- 0.03]", B) == GR_SUCCESS);
    XCHECK(_decball_contains_fmpq(X, q, B) == 0);
    decball_clear(X, B);

    /* working precision caps saturate instead of overflowing */
    XCHECK(_decimal_ziv_wp_max(100, 0) == 64 * 100 + 1000);
    XCHECK(_decimal_ziv_wp_max(DECIMAL_PREC_MAX, WORD_MAX) == WORD_MAX / 4);
    XCHECK(_decimal_wp_cap(DECIMAL_PREC_MAX, 100, 1000) == WORD_MAX / 4);

    /* exact ball contexts: no guard digits may be added to
       DECIMAL_PREC_EXACT */
    decball_init(X, BX);
    decball_init(Y, BX);
    XCHECK(decball_set_si(X, 4, BX) == GR_SUCCESS);
    XCHECK(decball_rsqrt(Y, X, BX) == GR_SUCCESS);
    XCHECK_STR(decball_get_str(Y, BX), "0.5");
    XCHECK(decball_exp(Y, X, BX) == GR_UNABLE);
    XCHECK(decball_set_si(X, 3, BX) == GR_SUCCESS);
    fmpz_set_ui(t, 10);
    XCHECK(decball_mul_2exp_fmpz(Y, X, t, BX) == GR_SUCCESS);
    XCHECK_STR(decball_get_str(Y, BX), "3072");
    fmpz_set_ui(t, 100000);
    XCHECK(decball_mul_2exp_fmpz(Y, X, t, BX) == GR_UNABLE);
    qqbar_set_ui(s, 2);
    qqbar_sqrt(s, s);
    XCHECK(decball_set_qqbar(Y, s, BX) == GR_UNABLE);
    qqbar_set_si(s, 1);
    qqbar_div_ui(s, s, 4);
    XCHECK(decball_set_qqbar(Y, s, BX) == GR_SUCCESS);
    XCHECK_STR(decball_get_str(Y, BX), "0.25");
    decball_clear(X, BX);
    decball_clear(Y, BX);

    deccball_init(Z, CBX);
    deccball_init(W, CBX);
    deccball_init(V, CBX);
    XCHECK(deccball_set_si(Z, 2, CBX) == GR_SUCCESS);
    XCHECK(decball_set_si(&W->re, 1, CBX) == GR_SUCCESS);
    XCHECK(decball_set_si(&W->im, 1, CBX) == GR_SUCCESS);
    XCHECK(deccball_div(V, Z, W, CBX) == GR_SUCCESS);
    XCHECK_STR(deccball_get_str(V, CBX), "(1 - 1*I)");
    XCHECK(deccball_cmpabs(&c, Z, W, CBX) == GR_SUCCESS);
    XCHECK(c == 1);
    XCHECK(deccball_exp(V, W, CBX) == GR_UNABLE);
    /* a Gaussian integer via the generic acb path keeps a tight enclosure */
    fmpz_set_str(fmpzi_realref(g), "1000000000000000000000000000001", 10);
    fmpz_set_si(fmpzi_imagref(g), 7);
    XCHECK(deccball_set_other(V, g, ZI, CBX) == GR_SUCCESS);
    XCHECK(_decmag_cmp_10exp_si(&V->re.rad, 10, CBX) < 0);
    deccball_clear(Z, CBX);
    deccball_clear(W, CBX);
    deccball_clear(V, CBX);

#if FLINT_BITS == 64
    /* exponent limits: scientific exponents beyond about WORD_MAX / 2 are
       compared exactly against emin and emax (they were declared
       overflows or underflows regardless of the limits) */
    {
        gr_ctx_t FL, BL;
        decball_t P, Q;

        gr_ctx_init_decfloat(FL, 20, 0);
        decimal_ctx_set_exp_limits(FL, -100, WORD_MAX);
        XCHECK(decfloat_set_str(x, "1e5000000000000000000", FL) == GR_SUCCESS);
        XCHECK_STR(decfloat_get_str(x, FL), "1e5000000000000000000");
        decimal_ctx_set_exp_limits(FL, WORD_MIN, 100);
        XCHECK(decfloat_set_str(x, "1e-5000000000000000000", FL) == GR_SUCCESS);
        decimal_ctx_set_exp_limits(FL, -WORD(9000000000000000000), WORD(9000000000000000000));
        XCHECK(decfloat_set_str(x, "1e5000000000000000000", FL) == GR_SUCCESS);
        XCHECK(decfloat_set_str(x, "1e-5000000000000000000", FL) == GR_SUCCESS);
        gr_ctx_clear(FL);

        /* a spurious flush to zero gave a ball not containing the result */
        gr_ctx_init_decball(BL, 20, DECIMAL_ALLOW_UNDERFLOW);
        decimal_ctx_set_exp_limits(BL, -WORD(9000000000000000000), 100);
        decball_init(P, BL);
        decball_init(Q, BL);
        XCHECK(decball_set_str(P, "1e-3000000000000000000", BL) == GR_SUCCESS);
        XCHECK(decball_mul(Q, P, P, BL) == GR_SUCCESS);
        XCHECK_STR(decball_get_str(Q, BL), "1e-6000000000000000000");
        decball_clear(P, BL);
        decball_clear(Q, BL);
        gr_ctx_clear(BL);
    }
#endif

    decfloat_clear(x, F);
    decfloat_clear(y, F);
    fmpz_clear(m);
    fmpz_clear(t);
    fmpq_clear(q);
    arb_clear(a);
    arf_clear(r);
    qqbar_clear(s);
    fmpzi_clear(g);

    gr_ctx_clear(F);
    gr_ctx_clear(B);
    gr_ctx_clear(BX);
    gr_ctx_clear(CBX);
    gr_ctx_clear(ZI);

    TEST_FUNCTION_END(state);
}
