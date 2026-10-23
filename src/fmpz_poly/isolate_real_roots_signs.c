/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpq.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "fmpz_poly/impl.h"
#include "arb.h"
#include "arb_poly.h"

/*
    Real root isolation by locating sign changes, certified by Descartes'
    rule of signs.

    Let V be the number of sign variations in the coefficients of p. Then p
    has at most V positive roots. If we find V disjoint intervals with a
    sign change of p at the endpoints (or exact roots), each of them must
    contain exactly one root and there can be no other positive roots.

    If p has only real roots, Descartes' rule of signs is exact, so this
    always succeeds given enough subdivision. This is useful for real-rooted
    polynomials with roots spread over many orders of magnitude (for
    example, Eulerian polynomials), where VCA is very slow because the
    polynomials transformed to small subintervals have huge coefficients,
    while evaluating p at a point in ball arithmetic is cheap.

    We start with the points 2^L, ..., 2^K (L, K being lower and upper
    bounds for the positive roots) and repeatedly bisect all intervals until
    the count is reached, or give up when the next round of bisections would
    make the number of evaluations exceed the given budget.
*/

typedef struct
{
    fmpz_t c;       /* the point c 2^k */
    slong k;
    int sl;         /* sign of p just to the left of the point */
    int sr;         /* sign of p just to the right of the point */
}
sign_point_struct;

typedef struct
{
    const fmpz * p;
    const fmpz * dp;
    slong len;
    arb_ptr pa;     /* p rounded to prec */
    arb_ptr dpa;
    slong prec;
    slong bits;
    slong evals;
}
sign_eval_struct;

/* Sign of (poly, len) at c 2^k, with poly rounded to prec in polya. */
static int
_sign_at(const fmpz * poly, arb_srcptr polya, slong len, slong bits,
    slong prec0, const fmpz_t c, slong k)
{
    arb_t x, y;
    arb_ptr t;
    slong prec, i;
    int s = 2;

    arb_init(x);
    arb_init(y);
    arb_set_fmpz(x, c);
    arb_mul_2exp_si(x, x, k);

    _arb_poly_evaluate_horner(y, polya, len, x, prec0);
    if (arb_is_positive(y))
        s = 1;
    else if (arb_is_negative(y))
        s = -1;

    for (prec = 2 * prec0; s == 2 && prec <= 4 * (bits + prec0); prec *= 2)
    {
        t = _arb_vec_init(len);
        for (i = 0; i < len; i++)
            arb_set_round_fmpz(t + i, poly + i, prec);
        _arb_poly_evaluate_horner(y, t, len, x, prec);
        _arb_vec_clear(t, len);

        if (arb_is_positive(y))
            s = 1;
        else if (arb_is_negative(y))
            s = -1;
    }

    if (s == 2)
    {
        fmpq_t q, v;
        fmpq_init(q);
        fmpq_init(v);
        fmpz_set(fmpq_numref(q), c);
        if (k >= 0)
            fmpq_mul_2exp(q, q, k);
        else
            fmpq_div_2exp(q, q, -k);
        _fmpz_poly_evaluate_fmpq(fmpq_numref(v), fmpq_denref(v), poly, len,
            fmpq_numref(q), fmpq_denref(q));
        s = fmpz_sgn(fmpq_numref(v));
        fmpq_clear(q);
        fmpq_clear(v);
    }

    arb_clear(x);
    arb_clear(y);
    return s;
}

static void
_sign_point_set_signs(sign_point_struct * P, sign_eval_struct * E)
{
    int s;

    E->evals++;
    s = _sign_at(E->p, E->pa, E->len, E->bits, E->prec, P->c, P->k);

    if (s != 0)
    {
        P->sl = P->sr = s;
    }
    else
    {
        /* exact root; p is squarefree so p' is nonzero there */
        E->evals++;
        s = _sign_at(E->dp, E->dpa, E->len - 1, E->bits, E->prec, P->c, P->k);
        P->sl = -s;
        P->sr = s;
    }
}

/* Isolates the positive roots of (pol, len) where pol[0] != 0, given the
   number V of sign variations; the output is in increasing order.
   Returns 0 on failure (the output is then undefined). */
static int
_isolate_positive_signs(fmpq * exact_roots, slong * n_exact, fmpz * c_array,
    slong * k_array, slong * n_interval, const fmpz * pol, slong len,
    slong V, slong budget)
{
    sign_eval_struct E[1];
    sign_point_struct * P, * Q;
    fmpz * tmp;
    slong i, j, np, nq, K, L, S, X, kk;
    fmpz_t a, b;
    int ok = 0;

    *n_exact = *n_interval = 0;

    if (V == 0)
        return 1;

    /* 2^L < positive roots < 2^K */
    tmp = _fmpz_vec_init(len);
    _fmpz_vec_set(tmp, pol, len);
    K = _fmpz_poly_scale_positive_roots_0_1(tmp, len);
    for (i = 0; i < len; i++)
        fmpz_set(tmp + i, pol + len - 1 - i);
    L = -_fmpz_poly_scale_positive_roots_0_1(tmp, len);

    /* the bounds show that there are no positive roots */
    if (K <= L)
    {
        _fmpz_vec_clear(tmp, len);
        return 1;
    }

    E->p = pol;
    E->len = len;
    _fmpz_poly_derivative(tmp, pol, len);
    E->dp = tmp;
    E->bits = FLINT_ABS(_fmpz_vec_max_bits(pol, len));
    E->prec = 64 + 2 * FLINT_BIT_COUNT(len);
    E->pa = _arb_vec_init(len);
    E->dpa = _arb_vec_init(len - 1);
    for (i = 0; i < len; i++)
        arb_set_round_fmpz(E->pa + i, pol + i, E->prec);
    for (i = 0; i < len - 1; i++)
        arb_set_round_fmpz(E->dpa + i, tmp + i, E->prec);
    E->evals = 0;

    fmpz_init(a);
    fmpz_init(b);

    np = K - L + 1;
    P = flint_malloc(np * sizeof(sign_point_struct));
    for (i = 0; i < np; i++)
    {
        fmpz_init_set_ui(P[i].c, 1);
        P[i].k = L + i;
        _sign_point_set_signs(P + i, E);
    }

    while (1)
    {
        /* count exact roots and sign changes */
        S = X = 0;
        for (i = 0; i < np; i++)
        {
            if (P[i].sl != P[i].sr)
                X++;
            if (i + 1 < np && P[i].sr != P[i + 1].sl)
                S++;
        }

        if (S + X == V)
        {
            ok = 1;
            break;
        }

        /* each round doubles the number of evaluations, so check the
           budget before starting a round */
        if (S + X > V || E->evals + (np - 1) > budget)
            break;

        /* bisect all intervals */
        nq = 2 * np - 1;
        Q = flint_malloc(nq * sizeof(sign_point_struct));

        for (i = 0; i < np; i++)
        {
            fmpz_init(Q[2 * i].c);
            fmpz_swap(Q[2 * i].c, P[i].c);
            Q[2 * i].k = P[i].k;
            Q[2 * i].sl = P[i].sl;
            Q[2 * i].sr = P[i].sr;

            if (i + 1 < np)
            {
                /* the midpoint of c1 2^k1 and c2 2^k2 */
                sign_point_struct * M = Q + 2 * i + 1;
                kk = FLINT_MIN(P[i].k, P[i + 1].k);
                fmpz_mul_2exp(a, Q[2 * i].c, Q[2 * i].k - kk);
                fmpz_mul_2exp(b, P[i + 1].c, P[i + 1].k - kk);
                fmpz_init(M->c);
                fmpz_add(M->c, a, b);
                M->k = kk - 1;
                j = fmpz_val2(M->c);
                fmpz_fdiv_q_2exp(M->c, M->c, j);
                M->k += j;
                _sign_point_set_signs(M, E);
            }
        }

        for (i = 0; i < np; i++)
            fmpz_clear(P[i].c);
        flint_free(P);
        P = Q;
        np = nq;
    }

    if (ok)
    {
        for (i = 0; i < np; i++)
        {
            if (P[i].sl != P[i].sr)
            {
                if (exact_roots != NULL)
                {
                    fmpz_set(fmpq_numref(exact_roots + *n_exact), P[i].c);
                    fmpz_one(fmpq_denref(exact_roots + *n_exact));
                    if (P[i].k >= 0)
                        fmpq_mul_2exp(exact_roots + *n_exact, exact_roots + *n_exact, P[i].k);
                    else
                        fmpq_div_2exp(exact_roots + *n_exact, exact_roots + *n_exact, -P[i].k);
                }
                (*n_exact)++;
            }

            if (i + 1 < np && P[i].sr != P[i + 1].sl)
            {
                /* the interval between the points is (c 2^k, (c + 1) 2^k)
                   where 2^k is its width */
                kk = FLINT_MIN(P[i].k, P[i + 1].k);
                fmpz_mul_2exp(a, P[i].c, P[i].k - kk);
                fmpz_mul_2exp(b, P[i + 1].c, P[i + 1].k - kk);
                fmpz_sub(b, b, a);
                j = fmpz_val2(b);
                FLINT_ASSERT(fmpz_bits(b) == j + 1);

                if (c_array != NULL && k_array != NULL)
                {
                    fmpz_fdiv_q_2exp(c_array + *n_interval, a, j);
                    k_array[*n_interval] = kk + j;
                }
                (*n_interval)++;
            }
        }
    }

    for (i = 0; i < np; i++)
        fmpz_clear(P[i].c);
    flint_free(P);
    fmpz_clear(a);
    fmpz_clear(b);
    _arb_vec_clear(E->pa, len);
    _arb_vec_clear(E->dpa, len - 1);
    _fmpz_vec_clear(tmp, len);

    return ok;
}

int
_fmpz_poly_isolate_real_roots_signs(fmpq * exact_roots, slong * n_exact,
    fmpz * c_array, slong * k_array, slong * n_interval,
    const fmpz * pol, slong len, int positive_only, slong budget)
{
    fmpz * q;
    slong i, n_zeros, Vp, Vn, ne, ni, tmp;
    int s, sp, sn, ok;

    for (n_zeros = 0; n_zeros < len && fmpz_is_zero(pol + n_zeros); n_zeros++)
        ;
    pol += n_zeros;
    len -= n_zeros;

    if (len == 0)
        flint_throw(FLINT_ERROR, "(%s): zero polynomial\n", __func__);

    /* sign variations of p(x) and p(-x) */
    Vp = Vn = 0;
    sp = sn = fmpz_sgn(pol);
    for (i = 1; i < len; i++)
    {
        s = fmpz_sgn(pol + i);
        if (s != 0)
        {
            if (s != sp)
            {
                Vp++;
                sp = s;
            }
            if (i % 2 == 1)
                s = -s;
            if (s != sn)
            {
                Vn++;
                sn = s;
            }
        }
    }

    *n_exact = *n_interval = 0;
    ok = 1;

    if (!positive_only && Vn != 0)
    {
        q = _fmpz_vec_init(len);
        for (i = 0; i < len; i++)
        {
            if (i % 2 == 1)
                fmpz_neg(q + i, pol + i);
            else
                fmpz_set(q + i, pol + i);
        }

        ok = _isolate_positive_signs(exact_roots, &ne, c_array, k_array,
            &ni, q, len, Vn, budget);
        _fmpz_vec_clear(q, len);

        if (!ok)
            return 0;

        /* negate and reverse: (c 2^k, (c+1) 2^k) -> (-(c+1) 2^k, -c 2^k) */
        if (c_array != NULL && k_array != NULL)
        {
            for (i = 0; i < ni; i++)
            {
                fmpz_add_ui(c_array + i, c_array + i, 1);
                fmpz_neg(c_array + i, c_array + i);
            }
            for (i = 0; i < ni / 2; i++)
            {
                fmpz_swap(c_array + i, c_array + ni - 1 - i);
                tmp = k_array[i];
                k_array[i] = k_array[ni - 1 - i];
                k_array[ni - 1 - i] = tmp;
            }
        }

        if (exact_roots != NULL)
        {
            for (i = 0; i < ne; i++)
                fmpq_neg(exact_roots + i, exact_roots + i);
            for (i = 0; i < ne / 2; i++)
                fmpq_swap(exact_roots + i, exact_roots + ne - 1 - i);
        }

        *n_exact = ne;
        *n_interval = ni;
    }

    if (!positive_only)
    {
        if (exact_roots != NULL)
            for (i = 0; i < n_zeros; i++)
                fmpq_zero(exact_roots + *n_exact + i);
        *n_exact += n_zeros;
    }

    if (Vp != 0)
    {
        ok = _isolate_positive_signs(
            exact_roots == NULL ? NULL : exact_roots + *n_exact, &ne,
            c_array == NULL ? NULL : c_array + *n_interval,
            k_array == NULL ? NULL : k_array + *n_interval, &ni,
            pol, len, Vp, budget);

        *n_exact += ne;
        *n_interval += ni;
    }

    return ok;
}
