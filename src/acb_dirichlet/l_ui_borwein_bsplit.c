/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    L(s, chi) for integer s >= 1 by Borwein's algorithm (the Chebyshev
    acceleration of alternating series of Cohen, Rodriguez Villegas and
    Zagier) applied to alternating residue chains, with binary splitting.

    Let chi be primitive modulo q.

    If 4 | q, the 2-component of chi is primitive, which implies
    chi(n + q/2) = -chi(n), so with Q = q/2,

        L(s, chi) = sum_{a <= Q, gcd(a,q)=1} chi(a) T_a,
        T_a = sum_{j >= 0} (-1)^j (Q j + a)^(-s).

    If q is odd, then with Q = q, (-1)^(n-1) = (-1)^j (-1)^(a-1) for
    n = qj + a, so

        (1 - chi(2) 2^(1-s)) L(s, chi) = sum_n (-1)^(n-1) chi(n) n^(-s)
                                     = sum_{a <= q} (-1)^(a-1) chi(a) T_a.

    (There are no primitive characters with q = 2 mod 4.) The terms
    (Qj + a)^(-s) are the moments of a positive measure on [0,1] with total
    mass a^(-s), so Borwein's bound applies to each chain: the error is at
    most 3 a^(-s) / (3+sqrt(8))^n with n terms. The chains share the
    Borwein weights, so they are computed together by one binary splitting
    over vectors (or a few, in parallel).
*/

#include "acb.h"
#include "acb_dirichlet.h"
#include "thread_support.h"
#include "arb/impl.h"

/* log2(3+sqrt(8)) */
#define ERROR_B 2.5431066063272239453

typedef struct
{
    slong m;                    /* number of chains */
    arb_ptr A, B, Q2, Q3;       /* per chain */
    arb_t C, D, Q1;             /* shared */
}
lchain_state_struct;

static void
lchain_state_init(lchain_state_struct * S, slong m)
{
    S->m = m;
    S->A = _arb_vec_init(4 * m);
    S->B = S->A + m;
    S->Q2 = S->A + 2 * m;
    S->Q3 = S->A + 3 * m;
    arb_init(S->C);
    arb_init(S->D);
    arb_init(S->Q1);
}

static void
lchain_state_clear(lchain_state_struct * S)
{
    _arb_vec_clear(S->A, 4 * S->m);
    arb_clear(S->C);
    arb_clear(S->D);
    arb_clear(S->Q1);
}

static void
lchain_leaf(lchain_state_struct * S, slong k, slong n, ulong s, ulong Q, const ulong * res)
{
    slong i;

    arb_set_si(S->D, 2 * (n + k));
    arb_mul_si(S->D, S->D, n - k, ARF_PREC_EXACT);
    arb_set_si(S->Q1, k + 1);
    arb_mul_si(S->Q1, S->Q1, 2 * k + 1, ARF_PREC_EXACT);
    arb_set(S->C, S->Q1);

    for (i = 0; i < S->m; i++)
    {
        if (k == 0)
        {
            arb_zero(S->A + i);
            arb_one(S->Q2 + i);
        }
        else
        {
            /* the term (Q (k-1) + a)^(-s) with sign (-1)^(k-1) */
            if (k % 2)
                arb_set(S->A + i, S->Q1);
            else
                arb_neg(S->A + i, S->Q1);

            arb_ui_pow_ui(S->Q2 + i, Q * (k - 1) + res[i], s, ARF_PREC_EXACT);
        }

        arb_mul(S->Q3 + i, S->Q1, S->Q2 + i, ARF_PREC_EXACT);
        arb_zero(S->B + i);
    }
}

static void
lchain_bsplit(lchain_state_struct * L, slong a, slong b, slong n, ulong s,
    ulong Q, const ulong * res, int cont, slong bits)
{
    if (a + 1 == b)
    {
        lchain_leaf(L, a, n, s, Q, res);
    }
    else
    {
        lchain_state_struct R;
        slong m = (a + b) / 2, i;

        lchain_bsplit(L, m, b, n, s, Q, res, 1, bits);

        lchain_state_init(&R, L->m);
        lchain_bsplit(&R, a, m, n, s, Q, res, 1, bits);

        for (i = 0; i < L->m; i++)
        {
            arb_mul(L->B + i, L->B + i, R.D, bits);
            arb_addmul(L->B + i, L->A + i, R.C, bits);
            arb_mul(L->B + i, L->B + i, R.Q2 + i, bits);
            arb_addmul(L->B + i, R.B + i, L->Q3 + i, bits);

            arb_mul(L->A + i, L->A + i, R.Q3 + i, bits);
            arb_addmul(L->A + i, R.A + i, L->Q3 + i, bits);

            if (cont)
                arb_mul(L->Q2 + i, L->Q2 + i, R.Q2 + i, bits);

            arb_mul(L->Q3 + i, L->Q3 + i, R.Q3 + i, bits);
        }

        arb_mul(L->C, L->C, R.D, bits);
        arb_addmul(L->C, R.C, L->Q1, bits);

        if (cont)
            arb_mul(L->D, L->D, R.D, bits);

        arb_mul(L->Q1, L->Q1, R.Q1, bits);

        lchain_state_clear(&R);
    }
}

/* T_a for the chains a = res[0], ..., res[m-1] */
static void
lchain_eval(arb_ptr T, ulong s, ulong Q, const ulong * res, slong m, slong n, slong wp)
{
    lchain_state_struct S;
    arb_t t;
    slong i;

    arb_init(t);
    lchain_state_init(&S, m);
    lchain_bsplit(&S, 0, n + 1, n, s, Q, res, 0, wp);

    for (i = 0; i < m; i++)
    {
        /* T_a = (A C - B Q1) / (Q3 C) */
        arb_mul(S.A + i, S.A + i, S.C, wp);
        arb_mul(S.B + i, S.B + i, S.Q1, wp);
        arb_sub(S.A + i, S.A + i, S.B + i, wp);
        arb_mul(t, S.Q3 + i, S.C, wp);
        arb_div(T + i, S.A + i, t, wp);
    }

    lchain_state_clear(&S);
    arb_clear(t);
}

typedef struct
{
    arb_ptr T;
    ulong s;
    ulong Q;
    const ulong * res;
    slong m;
    slong blocks;
    slong n;
    slong wp;
}
lchain_work_struct;

static void
lchain_worker(slong i, void * arg)
{
    lchain_work_struct * W = (lchain_work_struct *) arg;
    slong lo, hi;

    lo = (W->m * i) / W->blocks;
    hi = (W->m * (i + 1)) / W->blocks;

    if (hi > lo)
        lchain_eval(W->T + lo, W->s, W->Q, W->res + lo, hi - lo, W->n, W->wp);
}

void
acb_dirichlet_l_ui_borwein_bsplit(acb_t res, ulong s, const dirichlet_group_t G,
    const dirichlet_char_t chi, slong prec)
{
    ulong q, Q, a, e;
    slong m, n, wp, i, blocks;
    ulong * resid;
    ulong * expo;
    arb_ptr T;
    acb_t c, z;
    mag_t err, t;
    int real;
    lchain_work_struct W;
    acb_dirichlet_roots_t roots;

    q = (G == NULL) ? 1 : G->q;
    real = (G == NULL) || dirichlet_char_is_real(G, chi);

    if (s == 0 || (G != NULL && !dirichlet_char_is_primitive(G, chi)) || q % 4 == 2)
    {
        acb_indeterminate(res);
        return;
    }

    /* 1 - chi(2) 2^(1-s) = 0 */
    if (s == 1 && q % 2 == 1 && (q == 1 || dirichlet_chi(G, chi, 2) == 0))
    {
        acb_indeterminate(res);
        return;
    }

    Q = (q % 4 == 0) ? q / 2 : q;

    wp = prec + 10 + 2 * FLINT_BIT_COUNT(q);
    if (s == 1 && G != NULL)
        wp += FLINT_BIT_COUNT(G->expo);  /* 1 - chi(2) may be small */

    n = wp / ERROR_B + 2;

    resid = flint_malloc(sizeof(ulong) * (Q + 1));
    expo = flint_malloc(sizeof(ulong) * (Q + 1));

    if (G != NULL)
        dirichlet_chi_vec(expo, G, chi, Q + 1);

    m = 0;
    for (a = 1; a <= Q; a++)
    {
        if (G == NULL)
            e = 0;
        else
            e = expo[a];

        if (e == DIRICHLET_CHI_NULL)
            continue;

        resid[m] = a;
        expo[m] = e;
        m++;
    }

    T = _arb_vec_init(m);

    /* compute the chains, in parallel blocks if worthwhile */
    /* the chains are independent; threads pay off already at a few
       thousand bits of total work (e.g. from about 700 bits for q = 7
       and at all precisions for q >= 100) */
    blocks = 1;
    if (m >= 2 && prec * m >= 4000)
        blocks = FLINT_MIN(m, flint_get_num_threads());

    W.T = T;
    W.s = s;
    W.Q = Q;
    W.res = resid;
    W.m = m;
    W.blocks = blocks;
    W.n = n;
    W.wp = wp;

    flint_parallel_do(lchain_worker, &W, blocks, -1, FLINT_PARALLEL_STRIDED);

    /* combine: sum of sign(a) chi(a) T_a */
    acb_init(c);
    acb_init(z);
    mag_init(err);
    mag_init(t);

    if (real)
    {
        arb_zero(acb_imagref(res));
        arb_zero(acb_realref(res));
    }
    else
    {
        acb_zero(res);
        acb_dirichlet_roots_init(roots, G->expo, m, wp);
    }

    for (i = 0; i < m; i++)
    {
        int neg = (q % 2 == 1) && (resid[i] % 2 == 0);

        if (real)
        {
            /* chi(a) = +/- 1 */
            if (expo[i] != 0)
                neg = !neg;

            if (neg)
                arb_sub(acb_realref(res), acb_realref(res), T + i, wp);
            else
                arb_add(acb_realref(res), acb_realref(res), T + i, wp);
        }
        else
        {
            acb_dirichlet_root(c, roots, expo[i], wp);
            if (neg)
                acb_neg(c, c);
            acb_addmul_arb(res, c, T + i, wp);
        }

        /* truncation error of the chain: 3 a^(-s) (3+sqrt(8))^(-n) */
        mag_set_ui_lower(t, resid[i]);
        mag_pow_ui_lower(t, t, s);
        mag_inv(t, t);
        mag_add(err, err, t);
    }

    mag_borwein_error(t, n);
    mag_mul(err, err, t);

    if (real)
        arb_add_error_mag(acb_realref(res), err);
    else
        acb_add_error_mag(res, err);

    if (q % 2 == 1)
    {
        /* divide by 1 - chi(2) 2^(1-s) */
        if (q == 1)
            acb_one(c);
        else
            acb_dirichlet_chi(c, G, chi, 2, wp);

        acb_mul_2exp_si(c, c, 1 - (slong) s);
        acb_sub_ui(c, c, 1, wp);
        acb_neg(c, c);
        acb_div(res, res, c, wp);
    }

    acb_set_round(res, res, prec);

    if (!real)
        acb_dirichlet_roots_clear(roots);

    _arb_vec_clear(T, m);
    flint_free(resid);
    flint_free(expo);
    acb_clear(c);
    acb_clear(z);
    mag_clear(err);
    mag_clear(t);
}
