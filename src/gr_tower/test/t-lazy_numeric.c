/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h>
#include <stdio.h>
#include "test_helpers.h"
#include "fmpq.h"
#include "acb.h"
#include "gr_special.h"
#include "gr_tower.h"
#include "gr_tower_lazy.h"

/*
    Random expressions with transcendental functions (exp, log, pi, sqrt,
    trigonometric functions, conjugates, real parts) of algebraic and
    transcendental numbers, evaluated in the lazy tower field and, in
    parallel, in ball arithmetic (independently, operation by operation):
    the enclosure of the tower value must overlap the ball; the zero test
    must never find zero a value whose ball excludes zero; identities
    which hold by construction (exp(log x) = x, sqrt(x)^2 = x, x = re x +
    i im x, and so on) must be decided.
*/

#define PREC 128

static char trace_num[10000];

static void
_trace_num(const char * s)
{
    if (strlen(trace_num) + strlen(s) < 9000)
        strcat(trace_num, s);
}

/* a random expression, evaluated in K and in balls; depth bounds the tree */
static int
_random_expr_num(gr_ptr x, acb_t z, slong depth, flint_rand_t state, gr_ctx_t K)
{
    int status = GR_SUCCESS;
    int op = n_randint(state, depth == 0 ? 4 : 18);
    char buf[64];

    flint_sprintf(buf, "%d ", op);
    _trace_num(buf);

    switch (op)
    {
        case 0:   /* small rational */
        {
            fmpq_t c;
            fmpq_init(c);
            fmpq_set_si(c, (slong) n_randint(state, 11) - 5, 1 + n_randint(state, 4));
            flint_sprintf(buf, "(%wd/%wd) ", (slong) fmpz_get_si(fmpq_numref(c)), (slong) fmpz_get_si(fmpq_denref(c)));
            _trace_num(buf);
            status |= gr_set_fmpq(x, c, K);
            acb_set_fmpq(z, c, PREC);
            fmpq_clear(c);
            break;
        }
        case 1:   /* i */
            status |= gr_i(x, K);
            acb_onei(z);
            break;
        case 2:   /* pi */
            status |= gr_pi(x, K);
            acb_const_pi(z, PREC);
            break;
        case 3:   /* sqrt of a small positive integer */
        {
            ulong a = 1 + n_randint(state, 12);
            flint_sprintf(buf, "(sqrt %wu) ", a);
            _trace_num(buf);
            status |= gr_set_ui(x, a, K);
            status |= gr_sqrt(x, x, K);
            acb_set_ui(z, a);
            acb_sqrt(z, z, PREC);
            break;
        }
        case 4: case 5:   /* sum, difference */
        {
            gr_ptr x2;
            acb_t z2;
            x2 = gr_heap_init(K);
            acb_init(z2);
            status |= _random_expr_num(x, z, depth - 1, state, K);
            status |= _random_expr_num(x2, z2, depth - 1, state, K);
            if (op == 4)
            {
                status |= gr_add(x, x, x2, K);
                acb_add(z, z, z2, PREC);
            }
            else
            {
                status |= gr_sub(x, x, x2, K);
                acb_sub(z, z, z2, PREC);
            }
            gr_heap_clear(x2, K);
            acb_clear(z2);
            break;
        }
        case 6: case 7:   /* product, quotient */
        {
            gr_ptr x2;
            acb_t z2;
            x2 = gr_heap_init(K);
            acb_init(z2);
            status |= _random_expr_num(x, z, depth - 1, state, K);
            status |= _random_expr_num(x2, z2, depth - 1, state, K);
            if (op == 6 || acb_contains_zero(z2))
            {
                status |= gr_mul(x, x, x2, K);
                acb_mul(z, z, z2, PREC);
            }
            else
            {
                status |= gr_div(x, x, x2, K);
                acb_div(z, z, z2, PREC);
            }
            gr_heap_clear(x2, K);
            acb_clear(z2);
            break;
        }
        case 8:   /* exp */
            status |= _random_expr_num(x, z, depth - 1, state, K);
            if (arb_is_negative(acb_realref(z)) || arf_cmpabs_2exp_si(arb_midref(acb_realref(z)), 6) < 0)
            {
                status |= gr_exp(x, x, K);
                acb_exp(z, z, PREC);
            }
            break;
        case 9:   /* log (away from zero) */
            status |= _random_expr_num(x, z, depth - 1, state, K);
            if (!acb_contains_zero(z))
            {
                status |= gr_log(x, x, K);
                acb_log(z, z, PREC);
            }
            break;
        case 10:  /* sqrt */
            status |= _random_expr_num(x, z, depth - 1, state, K);
            status |= gr_sqrt(x, x, K);
            acb_sqrt(z, z, PREC);
            break;
        case 11:  /* conjugate */
            status |= _random_expr_num(x, z, depth - 1, state, K);
            status |= gr_conj(x, x, K);
            acb_conj(z, z);
            break;
        case 12:  /* real part */
            status |= _random_expr_num(x, z, depth - 1, state, K);
            status |= gr_re(x, x, K);
            arb_zero(acb_imagref(z));
            break;
        case 13:  /* sin or cos */
            status |= _random_expr_num(x, z, depth - 1, state, K);
            if (arf_cmpabs_2exp_si(arb_midref(acb_imagref(z)), 6) < 0)
            {
                if (n_randint(state, 2))
                {
                    status |= gr_sin(x, x, K);
                    acb_sin(z, z, PREC);
                }
                else
                {
                    status |= gr_cos(x, x, K);
                    acb_cos(z, z, PREC);
                }
            }
            break;
        case 14:  /* absolute value */
            status |= _random_expr_num(x, z, depth - 1, state, K);
            status |= gr_abs(x, x, K);
            acb_abs(acb_realref(z), z, PREC);
            arb_zero(acb_imagref(z));
            break;
        case 15:  /* imaginary part */
            status |= _random_expr_num(x, z, depth - 1, state, K);
            status |= gr_im(x, x, K);
            arb_swap(acb_realref(z), acb_imagref(z));
            arb_zero(acb_imagref(z));
            break;
        case 16:  /* rational power (principal branch), away from zero */
        {
            fmpq_t e;
            acb_t w;
            fmpq_init(e);
            acb_init(w);
            /* (small exponents: powers of compound expressions grow fast) */
            fmpq_set_si(e, (slong) n_randint(state, 3) - 1, 2 + n_randint(state, 2));
            flint_sprintf(buf, "(^%wd/%wd) ", (slong) fmpz_get_si(fmpq_numref(e)), (slong) fmpz_get_si(fmpq_denref(e)));
            _trace_num(buf);
            status |= _random_expr_num(x, z, depth - 1, state, K);
            if (!acb_contains_zero(z) && !fmpq_is_zero(e))
            {
                status |= gr_pow_fmpq(x, x, e, K);
                acb_set_fmpq(w, e, PREC);
                acb_pow(z, z, w, PREC);
            }
            fmpq_clear(e);
            acb_clear(w);
            break;
        }
        case 17:  /* arctangent of a real number */
            status |= _random_expr_num(x, z, depth - 1, state, K);
            status |= gr_re(x, x, K);
            arb_zero(acb_imagref(z));
            status |= gr_atan(x, x, K);
            acb_atan(z, z, PREC);
            break;
    }

    return status;
}

TEST_FUNCTION_START(gr_tower_lazy_numeric, state)
{
    gr_ctx_t QQ, K;
    slong iter;

    gr_ctx_init_fmpq(QQ);
    gr_ctx_init_tower_lazy(K, QQ, 0);

    for (iter = 0; iter < 200 * flint_test_multiplier(); iter++)
    {
        gr_ptr x, y, d;
        acb_t z, w, v;
        slong depth = 1 + n_randint(state, 4);
        int status;

        /* (a fresh context every ten iterations; K is initialized before
           the loop so that it is valid whatever the test multiplier) */
        if (iter % 10 == 0 && iter > 0)
        {
            gr_ctx_clear(K);
            gr_ctx_init_tower_lazy(K, QQ, n_randint(state, 2) ? GR_TOWER_MERGE_EXPRESS : 0);
        }

        x = gr_heap_init(K);
        y = gr_heap_init(K);
        d = gr_heap_init(K);
        acb_init(z);
        acb_init(w);
        acb_init(v);

        trace_num[0] = 0;
        status = _random_expr_num(x, z, depth, state, K);

        if (status == GR_SUCCESS && acb_is_finite(z))
        {
            /* the enclosure overlaps the ball */
            if (gr_tower_lazy_get_acb(w, x, PREC, K) == GR_SUCCESS && !acb_overlaps(w, z))
            {
                flint_printf("FAIL: enclosures do not overlap\n%s\n", trace_num);
                gr_println(x, K);
                acb_printd(w, 20); flint_printf("\n");
                acb_printd(z, 20); flint_printf("\n");
                flint_abort();
            }

            /* a value whose ball excludes zero is not zero */
            if (!acb_contains_zero(z) && gr_is_zero(x, K) == T_TRUE)
            {
                flint_printf("FAIL: nonzero value found zero\n%s\n", trace_num);
                gr_println(x, K);
                acb_printd(z, 20); flint_printf("\n");
                flint_abort();
            }

            /* identities holding by construction */
            switch (n_randint(state, 7))
            {
                case 0:   /* exp(log(x)) = x */
                    if (!acb_contains_zero(z))
                    {
                        status = gr_log(y, x, K);
                        status |= gr_exp(y, y, K);
                        status |= gr_sub(d, y, x, K);
                        if (status == GR_SUCCESS && gr_is_zero(d, K) != T_TRUE)
                        {
                            flint_printf("FAIL: exp(log(x)) != x\n%s\n", trace_num);
                            gr_println(x, K); gr_println(d, K); { slong lv; gr_tower_print(gr_tower_lazy_get_tower(&lv, d, K)); }
                            flint_abort();
                        }
                    }
                    break;
                case 1:   /* sqrt(x)^2 = x */
                    status = gr_sqrt(y, x, K);
                    status |= gr_sqr(y, y, K);
                    status |= gr_sub(d, y, x, K);
                    if (status == GR_SUCCESS && gr_is_zero(d, K) != T_TRUE)
                    {
                        flint_printf("FAIL: sqrt(x)^2 != x\n%s\n", trace_num);
                        gr_println(x, K); gr_println(d, K); { slong lv; gr_tower_print(gr_tower_lazy_get_tower(&lv, d, K)); }
                        flint_abort();
                    }
                    break;
                case 2:   /* x = re(x) + i im(x) */
                    status = gr_re(y, x, K);
                    status |= gr_im(d, x, K);
                    {
                        gr_ptr ii;
                        ii = gr_heap_init(K);
                        status |= gr_i(ii, K);
                        status |= gr_mul(d, d, ii, K);
                        gr_heap_clear(ii, K);
                    }
                    status |= gr_add(y, y, d, K);
                    status |= gr_sub(d, y, x, K);
                    if (status == GR_SUCCESS && gr_is_zero(d, K) != T_TRUE)
                    {
                        flint_printf("FAIL: x != re(x) + i im(x)\n%s\n", trace_num);
                        gr_println(x, K); gr_println(d, K); { slong lv; gr_tower_print(gr_tower_lazy_get_tower(&lv, d, K)); }
                        flint_abort();
                    }
                    break;
                case 3:   /* conj(conj(x)) = x, x conj(x) = |x|^2 */
                    status = gr_conj(y, x, K);
                    status |= gr_conj(y, y, K);
                    status |= gr_sub(d, y, x, K);
                    if (status == GR_SUCCESS && gr_is_zero(d, K) != T_TRUE)
                    {
                        flint_printf("FAIL: conj(conj(x)) != x\n%s\n", trace_num);
                        gr_println(x, K); gr_println(d, K); { slong lv; gr_tower_print(gr_tower_lazy_get_tower(&lv, d, K)); }
                        flint_abort();
                    }
                    status = gr_conj(y, x, K);
                    status |= gr_mul(y, y, x, K);
                    status |= gr_abs(d, x, K);
                    status |= gr_sqr(d, d, K);
                    status |= gr_sub(d, d, y, K);
                    if (status == GR_SUCCESS && gr_is_zero(d, K) != T_TRUE)
                    {
                        flint_printf("FAIL: x conj(x) != |x|^2\n%s\n", trace_num);
                        gr_println(x, K); gr_println(d, K); { slong lv; gr_tower_print(gr_tower_lazy_get_tower(&lv, d, K)); }
                        flint_abort();
                    }
                    break;
                case 4:   /* sin(x)^2 + cos(x)^2 = 1 */
                    if (arf_cmpabs_2exp_si(arb_midref(acb_imagref(z)), 6) < 0)
                    {
                        status = gr_sin(y, x, K);
                        status |= gr_sqr(y, y, K);
                        status |= gr_cos(d, x, K);
                        status |= gr_sqr(d, d, K);
                        status |= gr_add(d, d, y, K);
                        status |= gr_sub_ui(d, d, 1, K);
                        if (status == GR_SUCCESS && gr_is_zero(d, K) != T_TRUE)
                        {
                            flint_printf("FAIL: sin^2 + cos^2 != 1\n%s\n", trace_num);
                            gr_println(x, K); gr_println(d, K); { slong lv; gr_tower_print(gr_tower_lazy_get_tower(&lv, d, K)); }
                            flint_abort();
                        }
                    }
                    break;
                case 6:   /* tan(x) = sin(x)/cos(x), atan(tan(x)) = x for small real x */
                    if (arb_is_zero(acb_imagref(z)) && arf_cmpabs_2exp_si(arb_midref(acb_realref(z)), 0) < 0 &&
                        mag_cmp_2exp_si(arb_radref(acb_realref(z)), -20) < 0)
                    {
                        status = gr_tan(y, x, K);
                        status |= gr_sin(d, x, K);
                        {
                            gr_ptr cc;
                            cc = gr_heap_init(K);
                            status |= gr_cos(cc, x, K);
                            status |= gr_div(d, d, cc, K);
                            gr_heap_clear(cc, K);
                        }
                        status |= gr_sub(d, d, y, K);
                        if (status == GR_SUCCESS && gr_is_zero(d, K) != T_TRUE)
                        {
                            flint_printf("FAIL: tan(x) != sin(x)/cos(x)\n%s\n", trace_num);
                            gr_println(x, K); gr_println(d, K); { slong lv; gr_tower_print(gr_tower_lazy_get_tower(&lv, d, K)); }
                            flint_abort();
                        }
                        status = gr_tan(y, x, K);
                        status |= gr_atan(y, y, K);
                        status |= gr_sub(d, y, x, K);
                        if (status == GR_SUCCESS && gr_is_zero(d, K) != T_TRUE)
                        {
                            flint_printf("FAIL: atan(tan(x)) != x\n%s\n", trace_num);
                            gr_println(x, K); gr_println(d, K); { slong lv; gr_tower_print(gr_tower_lazy_get_tower(&lv, d, K)); }
                            flint_abort();
                        }
                    }
                    break;
                case 5:   /* exp(x) exp(-x) = 1, log(exp(x)) = x for small x */
                    if (arf_cmpabs_2exp_si(arb_midref(acb_realref(z)), 6) < 0)
                    {
                        status = gr_exp(y, x, K);
                        status |= gr_neg(d, x, K);
                        status |= gr_exp(d, d, K);
                        status |= gr_mul(d, d, y, K);
                        status |= gr_sub_ui(d, d, 1, K);
                        if (status == GR_SUCCESS && gr_is_zero(d, K) != T_TRUE)
                        {
                            flint_printf("FAIL: exp(x) exp(-x) != 1\n%s\n", trace_num);
                            gr_println(x, K); gr_println(d, K); { slong lv; gr_tower_print(gr_tower_lazy_get_tower(&lv, d, K)); }
                            flint_abort();
                        }
                        if (arb_is_zero(acb_imagref(z)) ||
                            (arf_cmpabs_2exp_si(arb_midref(acb_imagref(z)), 1) < 0 && mag_cmp_2exp_si(arb_radref(acb_imagref(z)), -10) < 0))
                        {
                            /* Im(x) in (-pi, pi): log(exp(x)) = x */
                            status = gr_exp(y, x, K);
                            status |= gr_log(y, y, K);
                            status |= gr_sub(d, y, x, K);
                            if (status == GR_SUCCESS && gr_is_zero(d, K) != T_TRUE)
                            {
                                flint_printf("FAIL: log(exp(x)) != x\n%s\n", trace_num);
                                gr_println(x, K); gr_println(d, K); { slong lv; gr_tower_print(gr_tower_lazy_get_tower(&lv, d, K)); }
                                flint_abort();
                            }
                        }
                    }
                    break;
            }
        }

        gr_heap_clear(x, K);
        gr_heap_clear(y, K);
        gr_heap_clear(d, K);
        acb_clear(z);
        acb_clear(w);
        acb_clear(v);
    }

    /* rounding, csgn, cmpabs and get_si/get_ui on (possibly nonreal)
       elements whose real part is irrational, rational, or an integer or
       half-integer only after cancellation; checked against balls (the
       real part of x is exactly r = a + b (sqrt(2)^2 - 2) + c sqrt(2)) */
    {
        gr_ptr x, y, t, u;
        acb_t z;
        arb_t re, im, fl;
        GR_TMP_INIT4(x, y, t, u, K);
        acb_init(z);
        arb_init(re);
        arb_init(im);
        arb_init(fl);

        for (iter = 0; iter < 30 * flint_test_multiplier(); iter++)
        {
            fmpq_t a, c, e;
            slong b = (slong) n_randint(state, 5) - 2;
            int which, sg, imz, ok;
            slong n;
            ulong un;

            fmpq_init(a);
            fmpq_init(c);
            fmpq_init(e);
            fmpq_set_si(a, (slong) n_randint(state, 21) - 10, 1 + n_randint(state, 2));
            fmpq_set_si(c, n_randint(state, 3) == 0 ? (slong) n_randint(state, 5) - 2 : 0, 1 + n_randint(state, 3));
            fmpq_set_si(e, n_randint(state, 2) ? (slong) n_randint(state, 5) - 2 : 0, 1 + n_randint(state, 3));
            imz = fmpq_is_zero(e);

            /* x = a + b (sqrt(2)^2 - 2) + c sqrt(2) + e i sqrt(3) (+ b i (sqrt(3)^2 - 3)) */
            GR_MUST_SUCCEED(gr_set_ui(x, 2, K));
            GR_MUST_SUCCEED(gr_sqrt(x, x, K));
            GR_MUST_SUCCEED(gr_mul_fmpq(u, x, c, K));
            GR_MUST_SUCCEED(gr_sqr(x, x, K));
            GR_MUST_SUCCEED(gr_sub_ui(x, x, 2, K));
            GR_MUST_SUCCEED(gr_mul_si(x, x, b, K));
            GR_MUST_SUCCEED(gr_add_fmpq(x, x, a, K));
            GR_MUST_SUCCEED(gr_add(x, x, u, K));
            GR_MUST_SUCCEED(gr_set_ui(t, 3, K));
            GR_MUST_SUCCEED(gr_sqrt(t, t, K));
            GR_MUST_SUCCEED(gr_sqr(u, t, K));
            GR_MUST_SUCCEED(gr_sub_ui(u, u, 3, K));
            GR_MUST_SUCCEED(gr_mul_si(u, u, b, K));
            GR_MUST_SUCCEED(gr_mul_fmpq(t, t, e, K));
            GR_MUST_SUCCEED(gr_add(t, t, u, K));
            GR_MUST_SUCCEED(gr_i(u, K));
            GR_MUST_SUCCEED(gr_mul(t, t, u, K));
            GR_MUST_SUCCEED(gr_add(x, x, t, K));

            /* the exact real and imaginary parts as balls */
            arb_set_ui(re, 2);
            arb_sqrt(re, re, 4 * PREC);
            arb_mul_fmpz(re, re, fmpq_numref(c), 4 * PREC);
            arb_div_fmpz(re, re, fmpq_denref(c), 4 * PREC);
            {
                arb_t q;
                arb_init(q);
                arb_set_fmpq(q, a, 4 * PREC);
                arb_add(re, re, q, 4 * PREC);
                arb_clear(q);
            }
            arb_set_ui(im, 3);
            arb_sqrt(im, im, 4 * PREC);
            arb_mul_fmpz(im, im, fmpq_numref(e), 4 * PREC);
            arb_div_fmpz(im, im, fmpq_denref(e), 4 * PREC);

            which = n_randint(state, 4);
            if (which == 0)
                ok = gr_floor(y, x, K);
            else if (which == 1)
                ok = gr_ceil(y, x, K);
            else if (which == 2)
                ok = gr_trunc(y, x, K);
            else
                ok = gr_nint(y, x, K);
            if (ok != GR_SUCCESS)
            {
                flint_printf("FAIL: rounding (%d) status %d\n", which, ok);
                gr_println(x, K);
                flint_abort();
            }
            /* the result r0 must satisfy the defining inequalities */
            {
                fmpz_t r0;
                arb_t d;
                int bad = 0;
                fmpz_init(r0);
                arb_init(d);
                if (gr_get_fmpz(r0, y, K) != GR_SUCCESS)
                    bad = 1;
                else
                {
                    arb_sub_fmpz(d, re, r0, 4 * PREC);      /* d = re - r0 */
                    if (which == 0)         /* 0 <= d < 1 */
                    {
                        bad = arb_is_negative(d);
                        arb_sub_ui(d, d, 1, 4 * PREC);
                        bad = bad || arb_is_nonnegative(d);
                    }
                    else if (which == 1)    /* -1 < d <= 0 */
                    {
                        bad = arb_is_positive(d);
                        arb_add_ui(d, d, 1, 4 * PREC);
                        bad = bad || arb_is_nonpositive(d);
                    }
                    else if (which == 2)    /* |d| < 1, sign(d) = sign(re) or 0 */
                    {
                        bad = (arb_is_positive(re) && arb_is_negative(d)) || (arb_is_negative(re) && arb_is_positive(d));
                        arb_abs(d, d);
                        arb_sub_ui(d, d, 1, 4 * PREC);
                        bad = bad || arb_is_nonnegative(d);
                    }
                    else                    /* |d| <= 1/2 */
                    {
                        arb_abs(d, d);
                        arb_mul_2exp_si(d, d, 1);
                        arb_sub_ui(d, d, 1, 4 * PREC);
                        bad = arb_is_positive(d);
                    }
                }
                if (bad)
                {
                    flint_printf("FAIL: rounding (%d) value\n", which);
                    gr_println(x, K);
                    gr_println(y, K);
                    flint_abort();
                }
                fmpz_clear(r0);
                arb_clear(d);
            }

            /* csgn */
            GR_MUST_SUCCEED(gr_csgn(y, x, K));
            GR_MUST_SUCCEED(gr_get_si(&n, y, K));
            sg = arb_is_positive(re) ? 1 : arb_is_negative(re) ? -1 :
                 arb_is_positive(im) ? 1 : arb_is_negative(im) ? -1 : 0;
            if ((arb_is_exact(re) || !arb_contains_zero(re)) && (arb_is_exact(im) || !arb_contains_zero(im) || !arb_contains_zero(re)) && n != sg)
            {
                flint_printf("FAIL: csgn = %wd, expected %d\n", n, sg);
                gr_println(x, K);
                flint_abort();
            }

            /* get_si, get_ui: integers exactly when the real part is an
               integer and the imaginary part zero */
            ok = gr_get_si(&n, x, K);
            if (imz && fmpq_is_zero(c) && fmpz_is_one(fmpq_denref(a)))
            {
                if (ok != GR_SUCCESS || n != fmpz_get_si(fmpq_numref(a)))
                {
                    flint_printf("FAIL: get_si\n");
                    gr_println(x, K);
                    flint_abort();
                }
                ok = gr_get_ui(&un, x, K);
                if ((fmpz_sgn(fmpq_numref(a)) >= 0) != (ok == GR_SUCCESS))
                {
                    flint_printf("FAIL: get_ui\n");
                    gr_println(x, K);
                    flint_abort();
                }
            }
            else if (ok != GR_DOMAIN)
            {
                flint_printf("FAIL: get_si of a non-integer (%d)\n", ok);
                gr_println(x, K);
                flint_abort();
            }

            /* cmpabs with |1 + i| = sqrt(2) */
            GR_MUST_SUCCEED(gr_set_str(y, "1+i", K));
            {
                int cres;
                arb_t m;
                arb_init(m);
                GR_MUST_SUCCEED(gr_cmpabs(&cres, x, y, K));
                arb_sqr(m, re, 4 * PREC);
                arb_addmul(m, im, im, 4 * PREC);
                arb_sub_ui(m, m, 2, 4 * PREC);
                if ((arb_is_positive(m) && cres != 1) || (arb_is_negative(m) && cres != -1) ||
                    (arb_is_zero(m) && cres != 0) || (cres == 0 && !arb_contains_zero(m)))
                {
                    flint_printf("FAIL: cmpabs = %d\n", cres);
                    gr_println(x, K);
                    flint_abort();
                }
                arb_clear(m);
            }

            fmpq_clear(a);
            fmpq_clear(c);
            fmpq_clear(e);
        }

        GR_TMP_CLEAR4(x, y, t, u, K);
        acb_clear(z);
        arb_clear(re);
        arb_clear(im);
        arb_clear(fl);
    }

    gr_ctx_clear(K);
    gr_ctx_clear(QQ);

    TEST_FUNCTION_END(state);
}
