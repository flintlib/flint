/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"

/* With parameter n, the error is bounded by 3/(3+sqrt(8))^n */
#define ERROR_A 1.5849625007211561815 /* log2(3) */
#define ERROR_B 2.5431066063272239453 /* log2(3+sqrt(8)) */

typedef struct
{
    arb_t A;
    arb_t B;
    arb_t C;
    arb_t D;
    arb_t Q1;
    arb_t Q2;
    arb_t Q3;
}
zeta_bsplit_state;

typedef zeta_bsplit_state zeta_bsplit_t[1];

static inline void
zeta_bsplit_init(zeta_bsplit_t S)
{
    arb_init(S->A);
    arb_init(S->B);
    arb_init(S->C);
    arb_init(S->D);
    arb_init(S->Q1);
    arb_init(S->Q2);
    arb_init(S->Q3);
}

static inline void
zeta_bsplit_clear(zeta_bsplit_t S)
{
    arb_clear(S->A);
    arb_clear(S->B);
    arb_clear(S->C);
    arb_clear(S->D);
    arb_clear(S->Q1);
    arb_clear(S->Q2);
    arb_clear(S->Q3);
}


static inline void
zeta_coeff_k(zeta_bsplit_t S, slong k, slong n, slong s)
{
    arb_set_si(S->D, 2 * (n + k));
    arb_mul_si(S->D, S->D, n - k, ARF_PREC_EXACT);
    arb_set_si(S->Q1, k + 1);
    arb_mul_si(S->Q1, S->Q1, 2*k + 1, ARF_PREC_EXACT);

    if (k == 0)
    {
        arb_zero(S->A);
        arb_one(S->Q2);
    }
    else
    {
        arb_set_si(S->A, k % 2 ? 1 : -1);
        arb_mul(S->A, S->A, S->Q1, ARF_PREC_EXACT);
        arb_ui_pow_ui(S->Q2, k, s, ARF_PREC_EXACT);
    }

    arb_mul(S->Q3, S->Q1, S->Q2, ARF_PREC_EXACT);
    arb_zero(S->B);
    arb_set(S->C, S->Q1);
}

static void
zeta_bsplit(zeta_bsplit_t L, slong a, slong b,
    slong n, slong s, int cont, slong bits)
{
    if (a + 1 == b)
    {
        zeta_coeff_k(L, a, n, s);
    }
    else
    {
        zeta_bsplit_t R;

        slong m = (a + b) / 2;

        zeta_bsplit(L, m, b, n, s, 1, bits);

        zeta_bsplit_init(R);
        zeta_bsplit(R, a, m, n, s, 1, bits);

        arb_mul(L->B, L->B, R->D, bits);
        arb_addmul(L->B, L->A, R->C, bits);

        arb_mul(L->B, L->B, R->Q2, bits);
        arb_addmul(L->B, R->B, L->Q3, bits);

        arb_mul(L->A, L->A, R->Q3, bits);
        arb_addmul(L->A, R->A, L->Q3, bits);

        arb_mul(L->C, L->C, R->D, bits);
        arb_addmul(L->C, R->C, L->Q1, bits);

        if (cont)
        {
            arb_mul(L->D, L->D, R->D, bits);
            arb_mul(L->Q2, L->Q2, R->Q2, bits);
        }

        arb_mul(L->Q1, L->Q1, R->Q1, bits);
        arb_mul(L->Q3, L->Q3, R->Q3, bits);

        zeta_bsplit_clear(R);
    }
}

/* The error for eta(s) is bounded by 3/(3+sqrt(8))^n */
void
mag_borwein_error(mag_t err, slong n)
{
    /* upper bound for 1/(3+sqrt(8)) */
    mag_set_ui_2exp_si(err, 736899889, -32);

    mag_pow_ui(err, err, n);
    mag_mul_ui(err, err, 3);
}

void
arb_zeta_ui_borwein_bsplit(arb_t x, ulong s, slong prec)
{
    zeta_bsplit_t sum;
    mag_t err;
    slong wp, n;

    /* zeta(0) = -1/2 */
    if (s == 0)
    {
        arb_set_si(x, -1);
        arb_mul_2exp_si(x, x, -1);
        return;
    }

    if (s == 1)
    {
        flint_throw(FLINT_ERROR, "zeta_ui_borwein_bsplit: zeta(1)");
    }

    n = prec / ERROR_B + 2;
    wp = prec + 30;

    zeta_bsplit_init(sum);
    zeta_bsplit(sum, 0, n + 1, n, s, 0, wp);

    /*  A/Q3 - B/Q3 / (C/Q1) = (A*C - B*Q1) / (Q3*C)    */
    arb_mul(sum->A, sum->A, sum->C, wp);
    arb_mul(sum->B, sum->B, sum->Q1, wp);
    arb_sub(sum->A, sum->A, sum->B, wp);
    arb_mul(sum->Q3, sum->Q3, sum->C, wp);
    arb_div(sum->C, sum->A, sum->Q3, wp);

    mag_init(err);
    mag_borwein_error(err, n);
    mag_add(arb_radref(sum->C), arb_radref(sum->C), err);
    mag_clear(err);

    /* convert from eta(s) to zeta(s) */
    arb_div_2expm1_ui(x, sum->C, s - 1, wp);
    arb_mul_2exp_si(x, x, s - 1);

    zeta_bsplit_clear(sum);
}

