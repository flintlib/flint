/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "acb_poly.h"
#include "acb_hypgeom.h"

static void
evaluate(acb_poly_t A, acb_srcptr a, slong p, const acb_t z, slong n, slong prec)
{
    acb_poly_fit_length(A, p + 1);

    if (p == 1)
    {
        acb_add_ui(A->coeffs, a, n, prec);
        if (z != NULL)
            acb_mul(A->coeffs, A->coeffs, z, prec);
    }
    else if (p == 2)
    {
        acb_add(A->coeffs, a + 0, a + 1, prec);
        acb_add_ui(A->coeffs + 1, A->coeffs, 2 * n, prec);
        acb_add_ui(A->coeffs, A->coeffs, n, prec);
        acb_mul_ui(A->coeffs, A->coeffs, n, prec);
        acb_addmul(A->coeffs, a + 0, a + 1, prec);
        if (z != NULL)
        {
            acb_mul(A->coeffs, A->coeffs, z, prec);
            acb_mul(A->coeffs + 1, A->coeffs + 1, z, prec);
        }
    }
    else if (p == 3)
    {
        acb_t t, u;
        acb_init(t);
        acb_init(u);

        acb_add(t, a + 0, a + 1, prec);
        acb_add(t, t, a + 2, prec);

        acb_mul(u, a + 0, a + 1, prec);
        acb_mul(A->coeffs, u, a + 2, prec);

        acb_addmul(u, a + 0, a + 2, prec);
        acb_addmul(u, a + 1, a + 2, prec);

        /*
        (a0 + n)(a1 + n)(a2 + n) = a0 a1 a2 + (a0 a1 + a0 a2 + a1 a2) n + (a0 + a1 + a2) n^2 + n^3
        (a0 a1 + a0 a2 + a1 a2) + 2 (a0 + a1 + a2) n + 3 n^2
        (a0 + a1 + a2) + 3n
        1
        */

        acb_addmul_ui(A->coeffs, u, n, prec);
        acb_addmul_ui(A->coeffs, t, n * n, prec);
        acb_add_ui(A->coeffs, A->coeffs, n * n * n, prec);

        acb_set(A->coeffs + 1, u);
        acb_addmul_ui(A->coeffs + 1, t, 2 * n, prec);
        acb_add_ui(A->coeffs + 1, A->coeffs + 1, 3 * n * n, prec);

        acb_add_ui(A->coeffs + 2, t, 3 * n, prec);

        if (z != NULL)
        {
            acb_mul(A->coeffs + 0, A->coeffs + 0, z, prec);
            acb_mul(A->coeffs + 1, A->coeffs + 1, z, prec);
            acb_mul(A->coeffs + 2, A->coeffs + 2, z, prec);
        }

        acb_clear(t);
        acb_clear(u);
    }
    else if (p != 0)
    {
        flint_throw(FLINT_ERROR, "(%s)\n", __func__);
    }

    if (z != NULL)
        acb_set(A->coeffs + p, z);
    else
        acb_one(A->coeffs + p);

    _acb_poly_set_length(A, p + 1);
    _acb_poly_normalise(A);
}

/* todo: do this using underscore methods */
static void
bsplit(acb_poly_t A, acb_poly_t B, acb_poly_t C,
    acb_srcptr a, slong p, acb_srcptr b, slong q,
    const acb_t z, slong an, slong bn, slong prec)
{
    if (bn - an == 1)
    {
        evaluate(A, a, p, z, an, prec);
        evaluate(B, b, q, NULL, an, prec);
        acb_poly_set(C, B);
    }
    else if (bn - an == 2)  /* inlined */
    {
        acb_poly_t A2, B2;

        acb_poly_init(A2);
        acb_poly_init(B2);

        evaluate(A, a, p, z, an, prec);
        evaluate(A2, a, p, z, an + 1, prec);
        evaluate(B, b, q, NULL, an, prec);
        evaluate(B2, b, q, NULL, an + 1, prec);

        acb_poly_mul(C, B, B2, prec);
        acb_poly_set(B, C);
        acb_poly_mul(C, A, B2, prec);
        acb_poly_add(C, C, B, prec);
        acb_poly_mul(A2, A, A2, prec);
        acb_poly_swap(A, A2);

        acb_poly_clear(A2);
        acb_poly_clear(B2);
    }
    else
    {
        slong m = an + (bn - an) / 2;

        acb_poly_t A2, B2, C2, T;

        acb_poly_init(A2);
        acb_poly_init(B2);
        acb_poly_init(C2);
        acb_poly_init(T);

        bsplit(A, B, C, a, p, b, q, z, an, m, prec);
        bsplit(A2, B2, C2, a, p, b, q, z, m, bn, prec);

        acb_poly_mul(T, B2, C, prec);
        acb_poly_mul(C, A, C2, prec);
        acb_poly_add(C, C, T, prec);
        acb_poly_mul(C2, B, B2, prec);
        acb_poly_swap(B, C2);
        acb_poly_mul(B2, A, A2, prec);
        acb_poly_swap(A, B2);

        acb_poly_clear(A2);
        acb_poly_clear(B2);
        acb_poly_clear(C2);
        acb_poly_clear(T);
    }
}

void
acb_hypgeom_pfq_sum_fme(acb_t s, acb_t t,
    acb_srcptr a, slong p, acb_srcptr b, slong q,
    const acb_t z, slong n, slong prec)
{
    acb_poly_t A, B, C;
    acb_ptr ks, As, Bs, Cs;
    acb_t u, v;
    acb_ptr * tree;
    slong i, k, m, w;

    /* we compute to n-1 instead of n to avoid dividing by 0 in the
       denominator when computing a hypergeometric polynomial
       that terminates right before a pole */
    if (n > 4)
    {
        m = n_sqrt(n - 1) / 4;  /* tuning parameter */
        w = (n - 1) / FLINT_MAX(m, 1);
    }
    else
    {
        m = w = 0;
    }

    if (m < 1 || w < 1 || p > 3 || q > 3)
    {
        acb_hypgeom_pfq_sum_forward(s, t, a, p, b, q, z, n, prec);
        return;
    }

    acb_poly_init(A);
    acb_poly_init(B);
    acb_poly_init(C);

    acb_init(u);
    acb_init(v);

    ks = _acb_vec_init(w);
    As = _acb_vec_init(w);
    Bs = _acb_vec_init(w);
    Cs = _acb_vec_init(w);

    bsplit(A, B, C, a, p, b, q, z, 0, m, prec);

    for (i = 0; i < w; i++)
        acb_set_ui(ks + i, i * m);

    tree = _acb_poly_tree_alloc(w);
    _acb_poly_tree_build(tree, ks, w, prec);
    _acb_poly_evaluate_vec_fast_precomp(As, A->coeffs, A->length, tree, w, prec);
    _acb_poly_evaluate_vec_fast_precomp(Bs, B->coeffs, B->length, tree, w, prec);
    _acb_poly_evaluate_vec_fast_precomp(Cs, C->coeffs, C->length, tree, w, prec);
    _acb_poly_tree_free(tree, w);

    /* todo: use binary splitting here for improved numerical stability */
    for (i = 1; i < w; i++)
    {
        acb_mul(Cs, Cs, Bs + i, prec);
        acb_addmul(Cs, As, Cs + i, prec);
        acb_mul(As, As, As + i, prec);
        acb_mul(Bs, Bs, Bs + i, prec);
    }

    acb_div(s, Cs, Bs, prec);
    acb_div(t, As, Bs, prec);

    for (k = w * m; k < n && !acb_is_zero(t); k++)
    {
        acb_add(s, s, t, prec);

        if (p > 0)
        {
            acb_add_ui(u, a, k, prec);

            for (i = 1; i < p; i++)
            {
                acb_add_ui(v, a + i, k, prec);
                acb_mul(u, u, v, prec);
            }

            acb_mul(t, t, u, prec);
        }

        if (q > 0)
        {
            acb_add_ui(u, b, k, prec);

            for (i = 1; i < q; i++)
            {
                acb_add_ui(v, b + i, k, prec);
                acb_mul(u, u, v, prec);
            }

            acb_div(t, t, u, prec);
        }

        acb_mul(t, t, z, prec);
    }

    acb_clear(u);
    acb_clear(v);

    _acb_vec_clear(ks, w);
    _acb_vec_clear(As, w);
    _acb_vec_clear(Bs, w);
    _acb_vec_clear(Cs, w);

    acb_poly_clear(A);
    acb_poly_clear(B);
    acb_poly_clear(C);
}

