/*
    Copyright (C) 2012, 2022 Fredrik Johansson

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

static inline void
fmpz_poly_evaluate_si(fmpz_t y, const fmpz_poly_t poly, slong x)
{
    fmpz_set_si(y, x);
    fmpz_poly_evaluate_fmpz(y, poly, y);
}

static void
bsplit_recursive_fmpz(fmpz_t P, fmpz_t Q, fmpz_t B, fmpz_t T,
    const hypgeom_t hyp, slong a, slong b, int cont)
{
    if (b - a == 1)
    {
        if (a == 0)
        {
            fmpz_one(P);
            fmpz_one(Q);
        }
        else
        {
            fmpz_poly_evaluate_si(P, hyp->P, a);
            fmpz_poly_evaluate_si(Q, hyp->Q, a);
        }

        fmpz_poly_evaluate_si(B, hyp->B, a);
        fmpz_poly_evaluate_si(T, hyp->A, a);
        fmpz_mul(T, T, P);
    }
    else
    {
        slong m;
        fmpz_t P2, Q2, B2, T2;

        m = (a + b) / 2;

        fmpz_init(P2);
        fmpz_init(Q2);
        fmpz_init(B2);
        fmpz_init(T2);

        bsplit_recursive_fmpz(P, Q, B, T, hyp, a, m, 1);
        bsplit_recursive_fmpz(P2, Q2, B2, T2, hyp, m, b, 1);

        if (fmpz_is_one(B) && fmpz_is_one(B2))
        {
            fmpz_mul(T, T, Q2);
            fmpz_addmul(T, P, T2);
        }
        else
        {
            fmpz_mul(T, T, B2);
            fmpz_mul(T, T, Q2);
            fmpz_mul(T2, T2, B);
            fmpz_addmul(T, P, T2);
        }

        fmpz_mul(B, B, B2);
        fmpz_mul(Q, Q, Q2);
        if (cont)
            fmpz_mul(P, P, P2);

        fmpz_clear(P2);
        fmpz_clear(Q2);
        fmpz_clear(B2);
        fmpz_clear(T2);
    }
}

typedef struct
{
    arb_struct P;
    arb_struct Q;
    arb_struct B;
    arb_struct T;
    slong a;
    slong b;
}
bsplit_res_t;

typedef struct
{
    const hypgeom_struct * hyp;
    slong prec;
    slong a;
    slong b;
}
bsplit_args_t;

static void
bsplit_init(bsplit_res_t * x, void * args)
{
    arb_init(&x->P);
    arb_init(&x->Q);
    arb_init(&x->B);
    arb_init(&x->T);
}

static void
bsplit_clear(bsplit_res_t * x, void * args)
{
    arb_clear(&x->P);
    arb_clear(&x->Q);
    arb_clear(&x->B);
    arb_clear(&x->T);
}

static void
bsplit_basecase(bsplit_res_t * res, slong a, slong b, bsplit_args_t * args)
{
    fmpz_t PP, QQ, BB, TT;
    int cont;

    fmpz_init(PP);
    fmpz_init(QQ);
    fmpz_init(BB);
    fmpz_init(TT);

    cont = (b != args->b);

    bsplit_recursive_fmpz(PP, QQ, BB, TT, args->hyp, a, b, cont);

    arb_set_fmpz(&res->P, PP);
    arb_set_fmpz(&res->Q, QQ);
    arb_set_fmpz(&res->B, BB);
    arb_set_fmpz(&res->T, TT);

    res->a = a;
    res->b = b;

    fmpz_clear(PP);
    fmpz_clear(QQ);
    fmpz_clear(BB);
    fmpz_clear(TT);
}

/* res = left */
static void
bsplit_merge(bsplit_res_t * res, bsplit_res_t * left, bsplit_res_t * right, bsplit_args_t * args)
{
    arb_ptr P = &res->P;
    arb_ptr Q = &res->Q;
    arb_ptr B = &res->B;
    arb_ptr T = &res->T;

    arb_ptr P2 = &right->P;
    arb_ptr Q2 = &right->Q;
    arb_ptr B2 = &right->B;
    arb_ptr T2 = &right->T;

    slong prec = args->prec;

    slong b = right->b;

    int cont = b != args->b;

    if (res != left)
        flint_throw(FLINT_ERROR, "(%s)\n", __func__);

    if (arb_is_one(B) && arb_is_one(B2))
    {
        arb_mul(T, T, Q2, prec);
        arb_addmul(T, P, T2, prec);
    }
    else
    {
        arb_mul(T, T, B2, prec);
        arb_mul(T, T, Q2, prec);
        arb_mul(T2, T2, B, prec);
        arb_addmul(T, P, T2, prec);
    }

    arb_mul(B, B, B2, prec);
    arb_mul(Q, Q, Q2, prec);
    if (cont)
        arb_mul(P, P, P2, prec);

    res->b = right->b;
}

#define WANT_PARALLEL 1

static void
bsplit_recursive_arb(arb_t P, arb_t Q, arb_t B, arb_t T,
    const hypgeom_t hyp, slong a, slong b, int cont, slong prec)
{
    if (WANT_PARALLEL)
    {
        bsplit_res_t res;
        bsplit_args_t args;

        res.P = *P;
        res.Q = *Q;
        res.B = *B;
        res.T = *T;

        args.hyp = hyp;
        args.prec = prec;
        args.a = a;
        args.b = b;

        flint_parallel_binary_splitting(&res,
            (bsplit_basecase_func_t) bsplit_basecase,
            (bsplit_merge_func_t) bsplit_merge,
            sizeof(bsplit_res_t),
            (bsplit_init_func_t) bsplit_init,
            (bsplit_clear_func_t) bsplit_clear,
            &args, a, b, 4, -1, FLINT_PARALLEL_BSPLIT_LEFT_INPLACE);

        *P = res.P;
        *Q = res.Q;
        *B = res.B;
        *T = res.T;
    }
    else
    {
        if (b - a < 4)
        {
            fmpz_t PP, QQ, BB, TT;

            fmpz_init(PP);
            fmpz_init(QQ);
            fmpz_init(BB);
            fmpz_init(TT);

            bsplit_recursive_fmpz(PP, QQ, BB, TT, hyp, a, b, cont);

            arb_set_fmpz(P, PP);
            arb_set_fmpz(Q, QQ);
            arb_set_fmpz(B, BB);
            arb_set_fmpz(T, TT);

            fmpz_clear(PP);
            fmpz_clear(QQ);
            fmpz_clear(BB);
            fmpz_clear(TT);
        }
        else
        {
            slong m;
            arb_t P2, Q2, B2, T2;

            m = (a + b) / 2;

            arb_init(P2);
            arb_init(Q2);
            arb_init(B2);
            arb_init(T2);

            bsplit_recursive_arb(P, Q, B, T, hyp, a, m, 1, prec);
            bsplit_recursive_arb(P2, Q2, B2, T2, hyp, m, b, 1, prec);

            if (arb_is_one(B) && arb_is_one(B2))
            {
                arb_mul(T, T, Q2, prec);
                arb_addmul(T, P, T2, prec);
            }
            else
            {
                arb_mul(T, T, B2, prec);
                arb_mul(T, T, Q2, prec);
                arb_mul(T2, T2, B, prec);
                arb_addmul(T, P, T2, prec);
            }

            arb_mul(B, B, B2, prec);
            arb_mul(Q, Q, Q2, prec);
            if (cont)
                arb_mul(P, P, P2, prec);

            arb_clear(P2);
            arb_clear(Q2);
            arb_clear(B2);
            arb_clear(T2);
        }
    }
}

void
arb_hypgeom_sum(arb_t P, arb_t Q, const hypgeom_t hyp, slong n, slong prec)
{
    if (n < 1)
    {
        arb_zero(P);
        arb_one(Q);
    }
    else
    {
        arb_t B, T;
        arb_init(B);
        arb_init(T);
        bsplit_recursive_arb(P, Q, B, T, hyp, 0, n, 0, prec);
        if (!arb_is_one(B))
            arb_mul(Q, Q, B, prec);
        arb_swap(P, T);
        arb_clear(B);
        arb_clear(T);
    }
}

void
arb_hypgeom_infsum(arb_t P, arb_t Q, hypgeom_t hyp, slong target_prec, slong prec)
{
    mag_t err, z;
    slong n;

    mag_init(err);
    mag_init(z);

    mag_set_fmpz(z, hyp->P->coeffs + hyp->P->length - 1);
    mag_div_fmpz(z, z, hyp->Q->coeffs + hyp->Q->length - 1);

    if (!hyp->have_precomputed)
    {
        hypgeom_precompute(hyp);
        hyp->have_precomputed = 1;
    }

    n = hypgeom_bound(err, hyp->r, hyp->boundC, hyp->boundD,
        hyp->boundK, hyp->MK, z, target_prec);

    arb_hypgeom_sum(P, Q, hyp, n, prec);

    if (arf_sgn(arb_midref(Q)) < 0)
    {
        arb_neg(P, P);
        arb_neg(Q, Q);
    }

    /* We have p/q = s + err i.e. (p + q*err)/q = s */
    {
        mag_t u;
        mag_init(u);
        arb_get_mag(u, Q);
        mag_mul(u, u, err);
        mag_add(arb_radref(P), arb_radref(P), u);
        mag_clear(u);
    }

    mag_clear(z);
    mag_clear(err);
}

