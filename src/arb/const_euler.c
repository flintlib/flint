/*
    Copyright (C) 2012, 2013, 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "thread_support.h"
#include "fmpz_poly.h"
#include "arb.h"
#include "hypgeom.h"
#include "arb_hypgeom.h"

typedef struct
{
    arb_t P;
    arb_t Q;
    arb_t T;
    arb_t C;
    arb_t D;
    arb_t V;
    slong a;
    slong b;
} euler_bsplit_1_struct;

typedef euler_bsplit_1_struct euler_bsplit_1_t[1];

static void euler_bsplit_1_init(euler_bsplit_1_t s, void * args)
{
    arb_init(s->P);
    arb_init(s->Q);
    arb_init(s->T);
    arb_init(s->C);
    arb_init(s->D);
    arb_init(s->V);
}

static void euler_bsplit_1_clear(euler_bsplit_1_t s, void * args)
{
    arb_clear(s->P);
    arb_clear(s->Q);
    arb_clear(s->T);
    arb_clear(s->C);
    arb_clear(s->D);
    arb_clear(s->V);
}

typedef struct
{
    slong N;
    slong prec;
    slong a;
    slong b;
}
bsplit_args_t;

static void
euler_bsplit_1_merge(euler_bsplit_1_t s, euler_bsplit_1_t L, euler_bsplit_1_t R, bsplit_args_t * args)
{
    arb_t t, u, v;

    slong wp = args->prec;

    slong b = R->b;

    int cont = (b != args->b);

    arb_init(t);
    arb_init(u);
    arb_init(v);

    if (cont)
        arb_mul(s->P, L->P, R->P, wp);

    arb_mul(s->Q, L->Q, R->Q, wp);
    arb_mul(s->D, L->D, R->D, wp);

    /* T = LP RT + RQ LT*/
    arb_mul(t, L->P, R->T, wp);
    arb_mul(v, R->Q, L->T, wp);
    arb_add(s->T, t, v, wp);

    /* C = LC RD + RC LD */
    if (cont)
    {
        arb_mul(s->C, L->C, R->D, wp);
        arb_addmul(s->C, R->C, L->D, wp);
    }

    /* V = RD (RQ LV + LC LP RT) + LD LP RV */
    arb_mul(u, L->P, R->V, wp);
    arb_mul(u, u, L->D, wp);
    arb_mul(v, R->Q, L->V, wp);
    arb_addmul(v, t, L->C, wp);
    arb_mul(v, v, R->D, wp);
    arb_add(s->V, u, v, wp);

    arb_clear(t);
    arb_clear(u);
    arb_clear(v);

    s->a = L->a;
    s->b = R->b;
}

static void
euler_bsplit_1_basecase(euler_bsplit_1_t s, slong n1, slong n2, bsplit_args_t * args)
{
    if (n2 - n1 == 1)
    {
        slong wp = args->prec;

        arb_set_si(s->P, args->N); /* p = N^2 */
        arb_mul(s->P, s->P, s->P, wp);
        arb_set_si(s->Q, n1 + 1); /* q = (k + 1)^2 */
        arb_mul(s->Q, s->Q, s->Q, wp);
        arb_set_si(s->C, 1);
        arb_set_si(s->D, n1 + 1);
        arb_set(s->T, s->P);
        arb_set(s->V, s->P);

        s->a = n1;
        s->b = n2;
    }
    else
    {
        euler_bsplit_1_t L, R;
        slong m = n1 + (n2 - n1) / 2;

        euler_bsplit_1_init(L, args);
        euler_bsplit_1_init(R, args);
        euler_bsplit_1_basecase(L, n1, m, args);
        euler_bsplit_1_basecase(R, m, n2, args);
        euler_bsplit_1_merge(s, L, R, args);
        euler_bsplit_1_clear(L, args);
        euler_bsplit_1_clear(R, args);
    }
}

static void
euler_bsplit_1(euler_bsplit_1_t s, slong n1, slong n2, slong N, slong wp, int cont)
{
    bsplit_args_t args;

    args.N = N;
    args.prec = wp;
    args.a = n1;
    args.b = n2;

    flint_parallel_binary_splitting(s,
        (bsplit_basecase_func_t) euler_bsplit_1_basecase,
        (bsplit_merge_func_t) euler_bsplit_1_merge,
        sizeof(euler_bsplit_1_struct),
        (bsplit_init_func_t) euler_bsplit_1_init,
        (bsplit_clear_func_t) euler_bsplit_1_clear,
        &args, n1, n2, 4, -1, 0);
}

typedef struct
{
    arb_t P;
    arb_t Q;
    arb_t T;
    slong a;
    slong b;
} euler_bsplit_2_struct;

typedef euler_bsplit_2_struct euler_bsplit_2_t[1];

static void euler_bsplit_2_init(euler_bsplit_2_t s, void * args)
{
    arb_init(s->P);
    arb_init(s->Q);
    arb_init(s->T);
}

static void euler_bsplit_2_clear(euler_bsplit_2_t s, void * args)
{
    arb_clear(s->P);
    arb_clear(s->Q);
    arb_clear(s->T);
}

static void
euler_bsplit_2_merge(euler_bsplit_2_t s, euler_bsplit_2_t L, euler_bsplit_2_t R, bsplit_args_t * args)
{
    arb_ptr P = s->P;
    arb_ptr Q = s->Q;
    arb_ptr T = s->T;

    arb_ptr P2 = R->P;
    arb_ptr Q2 = R->Q;
    arb_ptr T2 = R->T;

    slong wp = args->prec;

    slong b = R->b;

    int cont = (b != args->b);

    arb_mul(T, T, Q2, wp);
    arb_mul(T2, T2, P, wp);
    arb_add(T, T, T2, wp);

    if (cont)
        arb_mul(P, P, P2, wp);

    arb_mul(Q, Q, Q2, wp);

    s->a = L->a;
    s->b = R->b;
}

static void
euler_bsplit_2_basecase(euler_bsplit_2_t s, slong n1, slong n2, bsplit_args_t * args)
{
    if (n2 - n1 == 1)
    {
        slong wp = args->prec;
        slong N = args->N;

        arb_ptr P = s->P;
        arb_ptr Q = s->Q;
        arb_ptr T = s->T;

        if (n2 - n1 != 1)
            flint_throw(FLINT_ERROR, "(%s)\n", __func__);

        if (n1 == 0)
        {
            arb_set_si(P, 1);
            arb_set_si(Q, 4 * N);
            arb_set_si(T, 1);
        }
        else
        {
            arb_si_pow_ui(P, 1 - 2*n1, 3, wp);
            arb_neg(P, P);

            arb_set_si(Q, 32 * n1);
            arb_mul_ui(Q, Q, N, wp);
            arb_mul_ui(Q, Q, N, wp);
        }

        arb_set(T, P);

        s->a = n1;
        s->b = n2;
    }
    else
    {
        euler_bsplit_2_t R;
        slong m = n1 + (n2 - n1) / 2;

        euler_bsplit_2_init(R, args);
        euler_bsplit_2_basecase(s, n1, m, args);
        euler_bsplit_2_basecase(R, m, n2, args);
        euler_bsplit_2_merge(s, s, R, args);
        euler_bsplit_2_clear(R, args);
    }
}

static void
euler_bsplit_2(arb_t P, arb_t Q, arb_t T, slong n1, slong n2,
                        slong N, slong wp, int cont)
{
    euler_bsplit_2_t s;
    bsplit_args_t args;

    args.N = N;
    args.prec = wp;
    args.a = n1;
    args.b = n2;

    *s->P = *P;
    *s->Q = *Q;
    *s->T = *T;

    flint_parallel_binary_splitting(s,
        (bsplit_basecase_func_t) euler_bsplit_2_basecase,
        (bsplit_merge_func_t) euler_bsplit_2_merge,
        sizeof(euler_bsplit_2_struct),
        (bsplit_init_func_t) euler_bsplit_2_init,
        (bsplit_clear_func_t) euler_bsplit_2_clear,
        &args, n1, n2, 4, -1, FLINT_PARALLEL_BSPLIT_LEFT_INPLACE);

    *P = *s->P;
    *Q = *s->Q;
    *T = *s->T;
}

static void
atanh_bsplit(arb_t s, ulong c, slong a, slong prec)
{
    arb_t t;
    hypgeom_t series;
    hypgeom_init(series);
    arb_init(t);

    fmpz_poly_set_ui(series->A, 1);
    fmpz_poly_set_coeff_ui(series->B, 0, 1);
    fmpz_poly_set_coeff_ui(series->B, 1, 2);
    fmpz_poly_set_ui(series->P, 1);
    fmpz_poly_set_ui(series->Q, c * c);

    arb_hypgeom_infsum(s, t, series, prec, prec);
    arb_mul_si(s, s, a, prec);
    arb_mul_ui(t, t, c, prec);
    arb_div(s, s, t, prec);

    arb_clear(t);
    hypgeom_clear(series);
}

static ulong
next_smooth(ulong n)
{
    ulong t, k;

    for (k = n; ; k++)
    {
        t = k;
        while (t % 2 == 0) t /= 2;
        while (t % 3 == 0) t /= 3;
        while (t % 5 == 0) t /= 5;
        if (t == 1)
            return k;
    }
}

static void
arb_log_ui_smooth(arb_t s, ulong n, slong prec)
{
    ulong m, i, j, k;
    arb_t t;

    m = n;
    i = j = k = 0;
    while (m % 2 == 0) { m /= 2; i++; }
    while (m % 3 == 0) { m /= 3; j++; }
    while (m % 5 == 0) { m /= 5; k++; }

    if (m != 1)
        flint_throw(FLINT_ERROR, "(%s)\n", __func__);

    arb_init(t);

    prec += FLINT_CLOG2(prec);

    atanh_bsplit(s, 31, 14*i + 22*j + 32*k, prec);
    atanh_bsplit(t, 49, 10*i + 16*j + 24*k, prec);
    arb_add(s, s, t, prec);
    atanh_bsplit(t, 161, 6*i + 10*j + 14*k, prec);
    arb_add(s, s, t, prec);

    arb_clear(t);
}

void
arb_const_euler_eval(arb_t res, slong prec)
{
    euler_bsplit_1_t sum;
    arb_t t, u, v, P2, T2, Q2;
    slong bits, wp, wp2, n, N, M;

    bits = prec + 10;
    n = 0.086643397569993163677 * bits + 1;  /* log(2) / 8 */

    /* round n to have many trailing zeros, speeding up arithmetic,
       and make it smooth to allow computing the logarithm cheaply */
    if (n > 256)
    {
        int b = FLINT_BIT_COUNT(n);
        n = next_smooth((n >> (b-4)) + 1) << (b-4);
    }
    else
    {
        n = next_smooth(n);
    }

    /* As shown in the paper, it is sufficient to take
       N >= alpha n + 1 where alpha = 3/W(3/e) = 4.970625759544... */
    {
        fmpz_t a;
        fmpz_init(a);
        fmpz_set_ui(a, n);
        fmpz_mul_ui(a, a, 4970626);
        fmpz_cdiv_q_ui(a, a, 1000000);
        fmpz_add_ui(a, a, 1);
        N = fmpz_get_ui(a);
        fmpz_clear(a);
    }

    M = 2 * n;

    wp  = bits   + 2 * FLINT_BIT_COUNT(n);
    wp2 = bits/2 + 2 * FLINT_BIT_COUNT(n);

    euler_bsplit_1_init(sum, NULL);
    arb_init(P2);
    arb_init(T2);
    arb_init(Q2);
    arb_init(t);
    arb_init(u);
    arb_init(v);

    /* Compute S0 = V / (Q D), I0 = 1 + T / Q */
    euler_bsplit_1(sum, 0, N, n, wp, 0);

    /* I0 = T / Q */
    arb_add(sum->T, sum->T, sum->Q, wp);

    /* Compute S0 / I0 = V / (D T) */
    arb_mul(t, sum->T, sum->D, wp);
    arb_div(res, sum->V, t, wp);

    /* Compute K0 (actually I_0(2n) K_0(2n)) = T2 / Q2 */
    euler_bsplit_2(P2, Q2, T2, 0, M, n, wp2, 0);

    /* Compute K0 / I^2 = Q^2 * T2 / (Q2 * T^2) */
    arb_set_round(t, sum->Q, wp2);
    arb_mul(t, t, t, wp2);
    arb_mul(t, t, T2, wp2);
    arb_set_round(u, sum->T, wp2);
    arb_mul(u, u, u, wp2);
    arb_mul(u, u, Q2, wp2);
    arb_div(t, t, u, wp2);

    arb_sub(res, res, t, wp);

    /* subtract log(n) */
    arb_log_ui_smooth(u, n, wp);
    arb_sub(res, res, u, wp);

    /* add error bound 24 exp(-8n) */
    {
        mag_t b;
        mag_init(b);

        /* exp(-8) < 737690121 * 2^-41 */
        mag_set_ui_2exp_si(b, 737690121, -41);

        mag_pow_ui(b, b, n);
        mag_mul_ui(b, b, 24);

        mag_add(arb_radref(res), arb_radref(res), b);
        mag_clear(b);
    }

    arb_clear(P2);
    arb_clear(T2);
    arb_clear(Q2);
    arb_clear(t);
    arb_clear(u);
    arb_clear(v);

    euler_bsplit_1_clear(sum, NULL);
}

ARB_DEF_CACHED_CONSTANT(arb_const_euler_brent_mcmillan, arb_const_euler_eval)

FLINT_DLL extern const mp_limb_t arb_hypgeom_gamma_tab_limbs[];

void
arb_const_euler(arb_t res, slong prec)
{
    if (prec < ARB_HYPGEOM_GAMMA_TAB_PREC - 16)
    {
        slong exp;
        mp_size_t n;

        n = ARB_HYPGEOM_GAMMA_TAB_PREC / FLINT_BITS;

        /* just reading the table is known to give the correct rounding */
        _arf_set_round_mpn(arb_midref(res), &exp, arb_hypgeom_gamma_tab_limbs + n, n, 0, prec, ARF_RND_NEAR);
        _fmpz_set_si_small(ARF_EXPREF(arb_midref(res)), exp);

        /* 1/2 ulp error */
        _fmpz_set_si_small(MAG_EXPREF(arb_radref(res)), exp - prec);
        MAG_MAN(arb_radref(res)) = MAG_ONE_HALF;
    }
    else
    {
        arb_const_euler_brent_mcmillan(res, prec);
    }
}

