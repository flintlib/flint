/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq.h"
#include "gr.h"
#include "gr_special.h"
#include "gr_tower.h"
#include "gr_tower_lazy.h"

/*
    Exponentials whose arguments have rational parts in pi i: x = exp(a C +
    b pi i), y = exp(c C + d pi i) with a/c = p/q are related by x^q =
    y^p exp((b q - d p) pi i). The relation is found by the zero test and
    must be expressed through exp(C / L) and a root of unity of the order
    of the arguments (degree phi(60) = 16 for the first pair), not through
    a radical of degree about q over the other exponential.
*/

static void
_et_check(const char * C, slong a1, slong a2, slong b1, slong b2,
                          slong c1, slong c2, slong d1, slong d2, slong maxdeg, gr_ctx_t K)
{
    gr_ptr x, y, u, v, w;
    fmpq_t a, b, c, d, r, s;
    slong p, q, deg, level;
    gr_tower_struct * T;
    truth_t t;

    GR_TMP_INIT5(x, y, u, v, w, K);
    fmpq_init(a); fmpq_init(b); fmpq_init(c); fmpq_init(d); fmpq_init(r); fmpq_init(s);
    fmpq_set_si(a, a1, a2); fmpq_set_si(b, b1, b2);
    fmpq_set_si(c, c1, c2); fmpq_set_si(d, d1, d2);

    GR_MUST_SUCCEED(gr_set_str(u, C, K));
    GR_MUST_SUCCEED(gr_set_str(v, "pi*i", K));

    GR_MUST_SUCCEED(gr_mul_fmpq(x, u, a, K));
    GR_MUST_SUCCEED(gr_mul_fmpq(w, v, b, K));
    GR_MUST_SUCCEED(gr_add(x, x, w, K));
    GR_MUST_SUCCEED(gr_exp(x, x, K));
    GR_MUST_SUCCEED(gr_mul_fmpq(y, u, c, K));
    GR_MUST_SUCCEED(gr_mul_fmpq(w, v, d, K));
    GR_MUST_SUCCEED(gr_add(y, y, w, K));
    GR_MUST_SUCCEED(gr_exp(y, y, K));

    /* a/c = p/q; x^q = y^p exp((b q - d p) pi i) */
    fmpq_div(r, a, c);
    p = fmpz_get_si(fmpq_numref(r));
    q = fmpz_get_si(fmpq_denref(r));
    fmpq_mul_si(r, b, q);
    fmpq_mul_si(s, d, p);
    fmpq_sub(r, r, s);
    GR_MUST_SUCCEED(gr_mul_fmpq(w, v, r, K));
    GR_MUST_SUCCEED(gr_exp(w, w, K));
    GR_MUST_SUCCEED(gr_pow_si(u, y, p, K));
    GR_MUST_SUCCEED(gr_mul(u, u, w, K));
    GR_MUST_SUCCEED(gr_pow_si(w, x, q, K));

    t = gr_equal(w, u, K);
    if (t != T_TRUE)
    {
        flint_printf("FAIL: relation (%d), C = %s\n", (int) t, C);
        gr_println(x, K);
        gr_println(y, K);
        flint_abort();
    }

    /* the field of x y */
    GR_MUST_SUCCEED(gr_mul(u, x, y, K));
    T = gr_tower_lazy_get_tower(&level, u, K);
    deg = gr_tower_degree(T);
    if (deg > maxdeg)
    {
        flint_printf("FAIL: degree %wd > %wd, C = %s\n", deg, maxdeg, C);
        gr_println(u, K);
        flint_abort();
    }

    GR_TMP_CLEAR5(x, y, u, v, w, K);
    fmpq_clear(a); fmpq_clear(b); fmpq_clear(c); fmpq_clear(d); fmpq_clear(r); fmpq_clear(s);
}

/*
    exp(a + r pi i) for a with pi in a denominator (so that the rational
    multiple of pi i is not split off), r = 1/2, -1/2, 1: found at the
    adjunction to be i exp(a), -i exp(a), -exp(a) (a modulus of degree
    one, not the reducible X^2 + exp(a)^2)
*/
static void
_et_quarter(gr_ctx_t K)
{
    static const char * r[3] = { "pi*i/2", "-pi*i/2", "pi*i" };
    static const char * q[3] = { "i", "-i", "-1" };
    gr_ptr a, x, y, u;
    slong k, level;
    gr_tower_struct * T;

    GR_TMP_INIT4(a, x, y, u, K);
    for (k = 0; k < 3; k++)
    {
        GR_MUST_SUCCEED(gr_set_str(a, "pi^2/(pi + 1)", K));
        GR_MUST_SUCCEED(gr_exp(x, a, K));
        GR_MUST_SUCCEED(gr_set_str(u, r[k], K));
        GR_MUST_SUCCEED(gr_add(a, a, u, K));
        GR_MUST_SUCCEED(gr_exp(y, a, K));
        GR_MUST_SUCCEED(gr_mul(u, x, y, K));
        T = gr_tower_lazy_get_tower(&level, u, K);
        if (gr_tower_degree(T) > 2)
        {
            flint_printf("FAIL: exp(a + %s) not linear over exp(a)\n", r[k]);
            gr_println(u, K);
            flint_abort();
        }
        GR_MUST_SUCCEED(gr_div(u, y, x, K));
        GR_MUST_SUCCEED(gr_set_str(a, q[k], K));
        if (gr_equal(u, a, K) != T_TRUE)
        {
            flint_printf("FAIL: exp(a + %s) / exp(a)\n", r[k]);
            flint_abort();
        }
    }
    GR_TMP_CLEAR4(a, x, y, u, K);
}

TEST_FUNCTION_START(gr_tower_exp_torsion, state)
{
    gr_ctx_t QQ, K;
    const char * C[4] = { "1", "pi", "sqrt(2)", "pi + sqrt(3)" };
    slong i;

    gr_ctx_init_fmpq(QQ);

    for (i = 0; i < 4; i++)
    {
        gr_ctx_init_tower_lazy(K, QQ, GR_TOWER_MERGE_EXPRESS);
        /* exp(16 C/225 + 2 pi i/15), exp(C/12 + 7 pi i/10): zeta_60 (and
           sqrt(2) or sqrt(3), whose step may stay of degree 2) */
        _et_check(C[i], 16, 225, 2, 15, 1, 12, 7, 10, 16 * 2, K);
        gr_ctx_clear(K);

        gr_ctx_init_tower_lazy(K, QQ, GR_TOWER_MERGE_EXPRESS);
        /* exp(C/7 + pi i/5), exp(C/9 + pi i/3): zeta_30 and i, so zeta_60 */
        _et_check(C[i], 1, 7, 1, 5, 1, 9, 1, 3, 16 * 2, K);
        gr_ctx_clear(K);
    }

    gr_ctx_init_tower_lazy(K, QQ, GR_TOWER_MERGE_EXPRESS);
    _et_quarter(K);
    gr_ctx_clear(K);

    gr_ctx_clear(QQ);
    TEST_FUNCTION_END(state);
}
