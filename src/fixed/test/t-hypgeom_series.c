/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "arb.h"
#include "fixed.h"

/* the series in signed-mpn form, pointing into the fmpz coefficients
   (limb arrays of the absolute values, stored in d) */
typedef struct
{
    fixed_hypgeom_series_struct s;
    fixed_hypgeom_int_struct * c;
    nn_ptr d;
}
test_series_t;

static void
_set_int(fixed_hypgeom_int_struct * x, nn_ptr * d, const fmpz_t v)
{
    fmpz_t a;
    slong n = fmpz_size(v);

    fmpz_init(a);
    fmpz_abs(a, v);
    if (n > 0)
        fmpz_get_ui_array(*d, n, a);
    x->d = *d;
    x->n = n;
    x->neg = (fmpz_sgn(v) < 0);
    *d += FLINT_MAX(n, 1);
    fmpz_clear(a);
}

static void
test_series_init(test_series_t * t, int power, const fmpz_t cP,
    const fmpz_t cQ, const fmpz_t cD, const fmpz * P, slong Plen,
    const fmpz * Q, slong Qlen, const fmpz * R, slong Rlen)
{
    slong i, tot = 3, limbs = 3;
    nn_ptr d;

    for (i = 0; i < Plen; i++)
        limbs += FLINT_MAX(fmpz_size(P + i), 1);
    for (i = 0; i < Qlen; i++)
        limbs += FLINT_MAX(fmpz_size(Q + i), 1);
    for (i = 0; i < Rlen; i++)
        limbs += FLINT_MAX(fmpz_size(R + i), 1);
    limbs += fmpz_size(cP) + fmpz_size(cQ) + fmpz_size(cD);
    tot += Plen + Qlen + Rlen;

    t->c = flint_malloc(tot * sizeof(fixed_hypgeom_int_struct));
    t->d = d = flint_malloc(limbs * sizeof(ulong));

    for (i = 0; i < Plen; i++)
        _set_int(t->c + i, &d, P + i);
    for (i = 0; i < Qlen; i++)
        _set_int(t->c + Plen + i, &d, Q + i);
    for (i = 0; i < Rlen; i++)
        _set_int(t->c + Plen + Qlen + i, &d, R + i);
    _set_int(&t->s.coefP, &d, cP);
    _set_int(&t->s.coefQ, &d, cQ);
    _set_int(&t->s.coefD, &d, cD);

    t->s.power = power;
    t->s.P = t->c;
    t->s.Plen = Plen;
    t->s.Q = t->c + Plen;
    t->s.Qlen = Qlen;
    t->s.R = t->c + Plen + Qlen;
    t->s.Rlen = Rlen;
}

static void
test_series_clear(test_series_t * t)
{
    flint_free(t->c);
    flint_free(t->d);
}

/* reference: M terms summed in arb, plus a radius 2^-(M-1) |coefP| for
   the tail (|term_k| <= 2^-(k-1) by construction) */
static void
_reference(arb_t res, int power, const fmpz_t cP, const fmpz_t cQ,
    const fmpz_t cD, const fmpz * P, slong Plen, const fmpz * Q, slong Qlen,
    const fmpz * R, slong Rlen, slong M, slong prec)
{
    arb_t sum, prod, t, u;
    fmpz_t pk, qk, rk, kk;
    mag_t e;
    slong k;

    arb_init(sum); arb_init(prod); arb_init(t); arb_init(u);
    fmpz_init(pk); fmpz_init(qk); fmpz_init(rk); fmpz_init(kk);
    mag_init(e);

    arb_one(prod);
    for (k = 1; k <= M; k++)
    {
        fmpz_set_ui(kk, k);
        _fmpz_poly_evaluate_fmpz(pk, P, Plen, kk);
        _fmpz_poly_evaluate_fmpz(qk, Q, Qlen, kk);
        _fmpz_poly_evaluate_fmpz(rk, R, Rlen, kk);
        arb_mul_fmpz(t, prod, pk, prec);
        arb_div_fmpz(t, t, qk, prec);
        arb_add(sum, sum, t, prec);
        arb_mul_fmpz(prod, prod, rk, prec);
        arb_div_fmpz(prod, prod, qk, prec);
    }

    mag_set_fmpz(e, cP);
    mag_mul_2exp_si(e, e, -(M - 1));
    arb_mul_fmpz(sum, sum, cP, prec);
    arb_add_error_mag(sum, e);
    arb_add_fmpz(sum, sum, cQ, prec);
    arb_set_fmpz(u, cD);
    if (power == 1)
        arb_div(res, sum, u, prec);
    else
        arb_div(res, u, sum, prec);

    arb_clear(sum); arb_clear(prod); arb_clear(t); arb_clear(u);
    fmpz_clear(pk); fmpz_clear(qk); fmpz_clear(rk); fmpz_clear(kk);
    mag_clear(e);
}

static void
_randbig(fmpz_t x, flint_rand_t state, slong bits)
{
    fmpz_randbits(x, state, 1 + n_randint(state, bits));
    fmpz_abs(x, x);
}

TEST_FUNCTION_START(fixed_hypgeom_series, state)
{
    slong iter;

    /* pi (Chudnovsky, the y-cruncher formula) and log 2 through the
       int64 interface */
    for (iter = 0; iter < 30 * flint_test_multiplier(); iter++)
    {
        flint_set_num_threads(1 + n_randint(state, 4));
        static const int64_t piP[] = {-67957045, -2100495856,
            INT64_C(23608573992), INT64_C(-57896553024),
            INT64_C(39250089648)};
        static const int64_t piQ[] = {0, 0, 0, INT64_C(-10939058860032000)};
        static const int64_t piR[] = {-5, 46, -108, 72};
        static const int64_t l2P[] = {0, -1497, 1200, 3588};
        static const int64_t l2Q[] = {1080, 7776, 7776};
        static const int64_t l2R[] = {0, -1, 2};
        slong n = 1 + n_randint(state, (iter % 5 == 0) ? 500 : 50);
        fball_t x, y;
        arb_t a, b;

        fball_init(x); fball_init(y);
        arb_init(a); arb_init(b);

        if (iter % 2 == 0)
        {
            /* pi = S^-1 / sqrt(10005) */
            fball_hypgeom_series_int64(x, -1, 1, 13591409,
                INT64_C(4270934400), piP, 5, piQ, 4, piR, 4, n);
            fball_rsqrt_ui(y, 10005, n);
            fball_mul(x, x, y, n);
            arb_const_pi(b, FLINT_BITS * n + 64);
        }
        else
        {
            fball_hypgeom_series_int64(x, 1, 1, 1497, 2160, l2P, 4, l2Q, 3,
                l2R, 3, n);
            arb_const_log2(b, FLINT_BITS * n + 64);
        }

        fball_get_arb(a, x);
        if (!arb_overlaps(a, b)
            || arb_rel_accuracy_bits(a) < FLINT_BITS * (n - 1) - 32)
        {
            flint_printf("FAIL: constant %wd, n = %wd\n", iter % 2, n);
            arb_printd(a, 50); flint_printf("\n");
            arb_printd(b, 50); flint_printf("\n");
            flint_abort();
        }

        fball_clear(x); fball_clear(y);
        arb_clear(a); arb_clear(b);
    }

    /* atan(p/q), atanh(p/q) = (p/q) (1 + sum_k (-+p^2/q^2)^k / (2k+1)):
       P = R = -+p^2 (2k - 1), Q = q^2 (2k + 1), coefP = coefQ = p,
       coefD = q (content mode for large q) */
    for (iter = 0; iter < 200 * flint_test_multiplier(); iter++)
    {
        flint_set_num_threads(1 + n_randint(state, 4));
        slong n = 1 + n_randint(state, (iter % 10 == 0) ? 300 : 40);
        int hyp = n_randint(state, 2);
        fmpz * P, * Q;
        fmpz_t p, q, p2, q2;
        test_series_t t;
        fball_t x;
        arb_t a, b;

        P = _fmpz_vec_init(2);
        Q = _fmpz_vec_init(2);
        fmpz_init(p); fmpz_init(q); fmpz_init(p2); fmpz_init(q2);
        fball_init(x);
        arb_init(a); arb_init(b);

        /* 0 < p/q <= 0.98 */
        for (;;)
        {
            fmpz_randtest_unsigned(q, state, 2 + n_randint(state, 200));
            fmpz_add_ui(q, q, 2);
            if (n_randint(state, 2))
                fmpz_one(p);
            else
                fmpz_randm(p, state, q);
            fmpz_mul_ui(p2, p, 100);
            fmpz_mul_ui(q2, q, 98);
            if (!fmpz_is_zero(p) && fmpz_cmp(p2, q2) <= 0)
                break;
        }

        fmpz_mul(p2, p, p);
        fmpz_mul(q2, q, q);
        /* -+p^2 (2k - 1) = +-p^2 -+ 2 p^2 k */
        fmpz_set(P + 0, p2);
        fmpz_mul_si(P + 1, p2, -2);
        if (hyp)
            _fmpz_vec_neg(P, P, 2);
        fmpz_set(Q + 0, q2);
        fmpz_mul_ui(Q + 1, q2, 2);

        test_series_init(&t, 1, p, p, q, P, 2, Q, 2, P, 2);
        fball_hypgeom_series(x, &t.s, n);
        test_series_clear(&t);

        fball_get_arb(a, x);
        arb_set_fmpz(b, p);
        arb_div_fmpz(b, b, q, FLINT_BITS * n + 64);
        if (hyp)
            arb_atanh(b, b, FLINT_BITS * n + 64);
        else
            arb_atan(b, b, FLINT_BITS * n + 64);

        if (!arb_overlaps(a, b)
            || arb_rel_accuracy_bits(a) < FLINT_BITS * (n - 1) - 32)
        {
            flint_printf("FAIL: atan, hyp = %d, n = %wd\np = ", hyp, n);
            fmpz_print(p); flint_printf("\nq = "); fmpz_print(q);
            flint_printf("\n");
            arb_printd(a, 50); flint_printf("\n");
            arb_printd(b, 50); flint_printf("\n");
            flint_abort();
        }

        _fmpz_vec_clear(P, 2);
        _fmpz_vec_clear(Q, 2);
        fmpz_clear(p); fmpz_clear(q); fmpz_clear(p2); fmpz_clear(q2);
        fball_clear(x);
        arb_clear(a); arb_clear(b);
    }

    /* random series with coefficients of any size, constructed so that
       |P(k)/Q(k)| <= 1 and |R(k)/Q(k)| <= 1/2 (Q with positive
       coefficients, |r_i| <= q_i / 2, |p_i| <= q_i): sometimes R | P,
       a large content in Q or R, a terminating R, coefficients beyond
       the double range, power -1 */
    for (iter = 0; iter < 300 * flint_test_multiplier(); iter++)
    {
        flint_set_num_threads(1 + n_randint(state, 4));
        slong n = 1 + n_randint(state, (iter % 10 == 0) ? 60 : 12);
        slong d = n_randint(state, 4), Plen, Qlen = d + 1, Rlen, i, M;
        slong bits = (iter % 3 == 0) ? 300 : 60;
        int power = n_randint(state, 3) == 0 ? -1 : 1;
        int kind = n_randint(state, 6);
        fmpz * P, * Q, * R;
        fmpz_t cP, cQ, cD, c, u;
        test_series_t t;
        fball_t x;
        arb_t a, b;

        P = _fmpz_vec_init(Qlen);
        Q = _fmpz_vec_init(Qlen);
        R = _fmpz_vec_init(Qlen);
        fmpz_init(cP); fmpz_init(cQ); fmpz_init(cD);
        fmpz_init(c); fmpz_init(u);
        fball_init(x);
        arb_init(a); arb_init(b);

        /* Q: positive coefficients, a nonzero constant term (Q(k) >= 1) */
        for (i = 0; i < Qlen; i++)
        {
            _randbig(Q + i, state, bits);
            fmpz_add_ui(Q + i, Q + i, 2);
        }

        /* R: |r_i| <= q_i / 2, degree <= d */
        Rlen = 1 + n_randint(state, Qlen);
        for (i = 0; i < Rlen; i++)
        {
            fmpz_fdiv_q_2exp(u, Q + i, 1);
            fmpz_add_ui(u, u, 1);
            fmpz_randm(R + i, state, u);
            if (n_randint(state, 2))
                fmpz_neg(R + i, R + i);
        }
        if (fmpz_is_zero(R + Rlen - 1))
            fmpz_one(R + Rlen - 1);

        if (kind == 0 && d >= 1)
        {
            /* terminating: R = r1 (k - j0) with q1 >= 2 |r1|,
               q0 >= 2 |r1| j0 */
            slong j0 = 1 + n_randint(state, 30);
            fmpz_fdiv_q_2exp(u, Q + 1, 1);
            fmpz_add_ui(u, u, 1);
            fmpz_randm(R + 1, state, u);
            if (fmpz_is_zero(R + 1))
                fmpz_one(R + 1);
            fmpz_mul_ui(u, R + 1, 2 * j0);
            if (fmpz_cmp(Q + 0, u) < 0)
                fmpz_set(Q + 0, u);
            fmpz_mul_si(R + 0, R + 1, -j0);
            Rlen = 2;
        }

        /* P */
        if (kind == 1)
        {
            /* P = a R, |a| <= 2 */
            slong s = n_randint(state, 5) - 2;
            Plen = Rlen;
            _fmpz_vec_scalar_mul_si(P, R, Rlen, s ? s : 1);
        }
        else
        {
            Plen = 1 + n_randint(state, Qlen);
            for (i = 0; i < Plen; i++)
            {
                fmpz_add_ui(u, Q + i, 1);
                fmpz_randm(P + i, state, u);
                if (n_randint(state, 2))
                    fmpz_neg(P + i, P + i);
            }
        }

        if (kind == 2)
        {
            /* a large content in Q (and R, keeping |R| <= |Q|/2) */
            fmpz_randbits(c, state, 100 + n_randint(state, 500));
            fmpz_abs(c, c);
            fmpz_add_ui(c, c, 2);
            _fmpz_vec_scalar_mul_fmpz(Q, Q, Qlen, c);
            _fmpz_vec_scalar_mul_fmpz(P, P, Plen, c);
            if (n_randint(state, 2))
            {
                /* |R/Q| = |R'| / (c^0 |Q'|): keep R' c / c^0 */
                _fmpz_vec_scalar_mul_fmpz(R, R, Rlen, c);
            }
        }

        if (kind == 3)
        {
            /* huge coefficients (beyond the double range): a common
               factor of P, Q, R leaves the series unchanged */
            fmpz_randbits(c, state, 1000 + n_randint(state, 1500));
            fmpz_abs(c, c);
            fmpz_add_ui(c, c, 1);
            _fmpz_vec_scalar_mul_fmpz(Q, Q, Qlen, c);
            _fmpz_vec_scalar_mul_fmpz(P, P, Plen, c);
            _fmpz_vec_scalar_mul_fmpz(R, R, Rlen, c);
        }

        _randbig(cP, state, (iter % 2) ? 200 : 20);
        if (fmpz_is_zero(cP))
            fmpz_one(cP);
        if (n_randint(state, 2))
            fmpz_neg(cP, cP);
        /* |coefQ| >= 4 |coefP| keeps S away from 0 (|sum| <= 2) */
        _randbig(cQ, state, 100);
        fmpz_abs(u, cP);
        fmpz_addmul_ui(cQ, u, 4);
        if (n_randint(state, 2))
            fmpz_neg(cQ, cQ);
        _randbig(cD, state, (iter % 2) ? 200 : 20);
        if (fmpz_is_zero(cD))
            fmpz_one(cD);
        if (n_randint(state, 2))
            fmpz_neg(cD, cD);

        test_series_init(&t, power, cP, cQ, cD, P, Plen, Q, Qlen, R, Rlen);
        fball_hypgeom_series(x, &t.s, n);
        test_series_clear(&t);

        M = FLINT_BITS * n + fmpz_bits(cP) + fmpz_bits(cD) + 100;
        _reference(b, power, cP, cQ, cD, P, Plen, Q, Qlen, R, Rlen, M,
            FLINT_BITS * n + fmpz_bits(cP) + fmpz_bits(cD) + 100);
        fball_get_arb(a, x);

        if (!arb_overlaps(a, b)
            || arb_rel_accuracy_bits(a) < FLINT_BITS * (n - 1) - 40)
        {
            flint_printf("FAIL: random series, iter = %wd, n = %wd, "
                "kind = %d, power = %d\n", iter, n, kind, power);
            flint_printf("P = "); _fmpz_vec_print(P, Plen);
            flint_printf("\nQ = "); _fmpz_vec_print(Q, Qlen);
            flint_printf("\nR = "); _fmpz_vec_print(R, Rlen);
            flint_printf("\ncP = "); fmpz_print(cP);
            flint_printf(", cQ = "); fmpz_print(cQ);
            flint_printf(", cD = "); fmpz_print(cD);
            flint_printf("\n");
            arb_printd(a, 50); flint_printf("\n");
            arb_printd(b, 50); flint_printf("\n");
            flint_abort();
        }

        _fmpz_vec_clear(P, Qlen);
        _fmpz_vec_clear(Q, Qlen);
        _fmpz_vec_clear(R, Qlen);
        fmpz_clear(cP); fmpz_clear(cQ); fmpz_clear(cD);
        fmpz_clear(c); fmpz_clear(u);
        fball_clear(x);
        arb_clear(a); arb_clear(b);
    }

    flint_set_num_threads(1);
    TEST_FUNCTION_END(state);
}
