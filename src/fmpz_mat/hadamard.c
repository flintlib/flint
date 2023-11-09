/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "nmod_poly.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mat.h"
#include "fq_nmod.h"

static ulong
_fq_nmod_rank(const fq_nmod_t x, const fq_nmod_ctx_t ctx)
{
    slong i;
    ulong t = 0;

    for (i = x->length - 1; i >= 0; i--)
        t = t * ctx->mod.n + x->coeffs[i];

    return t;
}

static void
_fq_nmod_unrank(fq_nmod_t x, ulong r, const fq_nmod_ctx_t ctx)
{
    slong i;

    nmod_poly_zero(x);
    nmod_poly_fit_length(x, fq_nmod_ctx_degree(ctx));

    for (i = 0; r != 0; i++)
    {
        x->coeffs[i] = r % ctx->mod.n;
        x->length = i + 1;
        r = r / ctx->mod.n;
    }
}

static int
n_is_prime_power(mp_limb_t * p, mp_limb_t n)
{
    n_factor_t fac;

    if (n < 2)
        return 0;

    n_factor_init(&fac);
    n_factor(&fac, n, 1);

    if (fac.num == 1)
    {
        if (p != NULL)
            *p = fac.p[0];
        return fac.exp[0];
    }

    return 0;
}

/* Jacobsthal matrix of order q = p^d */
/* Could speed up greatly for d = 1. */
void
fmpz_mat_jacobsthal(fmpz_mat_t Q)
{
    int * quadratic;
    fmpz_t pp;
    ulong r, c, q, p, d;
    fq_nmod_ctx_t ctx;
    fq_nmod_t x, y, x2;

    q = fmpz_mat_nrows(Q);

    if (!(d = n_is_prime_power(&p, q)) || q % 2 == 0)
        flint_throw(FLINT_ERROR, "Not an odd prime power in %s\n", __func__);

    fmpz_init_set_ui(pp, p);
    fq_nmod_ctx_init(ctx, pp, d, "x");
    fq_nmod_init(x, ctx);
    fq_nmod_init(y, ctx);
    fq_nmod_init(x2, ctx);
    quadratic = flint_malloc(q * sizeof(int));

    for (r = 1; r < q; r++)
        quadratic[r] = -1;

    for (r = 1; r < q; r++)
    {
        _fq_nmod_unrank(x, r, ctx);
        fq_nmod_sqr(x2, x, ctx);
        quadratic[_fq_nmod_rank(x2, ctx)] = 1;
    }

    quadratic[0] = 0;

    for (r = 0; r < q; r++)
    {
        _fq_nmod_unrank(x, r, ctx);

        for (c = r; c < q; c++)
        {
            _fq_nmod_unrank(y, c, ctx);
            fq_nmod_sub(x2, x, y, ctx);
            fmpz_set_si(fmpz_mat_entry(Q, r, c),
                quadratic[_fq_nmod_rank(x2, ctx)]);

            if (q % 4 == 1)
                fmpz_set(fmpz_mat_entry(Q, c, r), fmpz_mat_entry(Q, r, c));
            else
                fmpz_neg(fmpz_mat_entry(Q, c, r), fmpz_mat_entry(Q, r, c));
        }
    }

    fq_nmod_clear(x, ctx);
    fq_nmod_clear(y, ctx);
    fq_nmod_clear(x2, ctx);
    fq_nmod_ctx_clear(ctx);
    flint_free(quadratic);
    fmpz_clear(pp);
}

/* 0 -- not possible */
/* 1 -- n = 2^v * (p^e + 1) */
/* 2 -- n = 2^v * 2*(p^e + 1) */
/* 3 -- n = 2^v */
static int
paley_construction(mp_limb_t * q, mp_limb_t n)
{
    int i, v;

    v = flint_ctz(n);

    if (UWORD(1) << v == n)
        return 3;

    if (n % 4 != 0)
        return 0;

    for (i = v - 1; i >= 0; i--)
    {
        *q = (n >> i) - 1;

        if (n_is_prime_power(NULL, *q) != 0)
        {
            if (*q % 4 == 3)
                return 1;
            else
                return 2;
        }
    }

    return 0;
}

static void
fmpz_mat_set2x2(fmpz_mat_t A, slong i, slong j,
                                  slong a, slong b, slong c, slong d)
{
    fmpz_set_si(fmpz_mat_entry(A, i, j), a);
    fmpz_set_si(fmpz_mat_entry(A, i, j + 1), b);
    fmpz_set_si(fmpz_mat_entry(A, i + 1, j), c);
    fmpz_set_si(fmpz_mat_entry(A, i + 1, j + 1), d);
}

int
fmpz_mat_hadamard(fmpz_mat_t A)
{
    slong n, m, i, j;
    mp_limb_t q;
    int kind;

    n = fmpz_mat_nrows(A);

    if (n != fmpz_mat_ncols(A))
        return 0;

    if (n == 0)
        return 1;

    kind = paley_construction(&q, n);

    if (kind == 0)
        return 0;

    if (kind == 3)
    {
        fmpz_one(fmpz_mat_entry(A, 0, 0));
        m = 1;
    }
    else
    {
        fmpz_mat_t Q;

        fmpz_mat_init(Q, q, q);
        fmpz_mat_jacobsthal(Q);

        if (kind == 1)
        {
            fmpz_zero(fmpz_mat_entry(A, 0, 0));

            for (i = 0; i < q; i++)
            {
                fmpz_set_si(fmpz_mat_entry(A, 0, i+1), 1);
                fmpz_set_si(fmpz_mat_entry(A, i+1, 0), -1);
            }

            for (i = 0; i < q; i++)
                for (j = 0; j < q; j++)
                    fmpz_set(fmpz_mat_entry(A, i+1, j+1),
                        fmpz_mat_entry(Q, i, j));

            for (i = 0; i < q + 1; i++)
                fmpz_add_ui(fmpz_mat_entry(A, i, i),
                    fmpz_mat_entry(A, i, i), 1);

            m = q + 1;
        }
        else
        {
            for (i = 0; i < q + 1; i++)
            {
                for (j = 0; j < q + 1; j++)
                {
                    if (i == j)
                        fmpz_mat_set2x2(A, 2 * i, 2 * j, 1, -1, -1, -1);
                    else if (i == 0 || j == 0
                            || fmpz_is_one(fmpz_mat_entry(Q, i - 1, j - 1)))
                        fmpz_mat_set2x2(A, 2 * i, 2 * j, 1, 1, 1, -1);
                    else
                        fmpz_mat_set2x2(A, 2 * i, 2 * j, -1, -1, -1, 1);
                }
            }

            m = 2 * (q + 1);
        }

        fmpz_mat_clear(Q);
    }

    for ( ; m < n; m *= 2)
    {
        for (i = 0; i < m; i++)
        {
            _fmpz_vec_set(A->rows[i] + m, A->rows[i], m);
            _fmpz_vec_set(A->rows[i + m], A->rows[i], m);
            _fmpz_vec_neg(A->rows[i + m] + m, A->rows[i], m);
        }
    }

    return 1;
}

