/*
    Copyright (C) 2007 David Howden
    Copyright (C) 2007, 2008, 2009, 2010 William Hart
    Copyright (C) 2008 Richard Howell-Peak
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2012 Lina Kulakova
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "templates.h"

#include "perm.h"

static void
TEMPLATE(T, to_mat_col) (TEMPLATE(T, mat_t) mat, slong col,
                         TEMPLATE(T, poly_t) poly,
                         const TEMPLATE(T, ctx_t) ctx)
{
    slong i;

    for (i = 0; i < poly->length; i++)
        TEMPLATE(T, set) (TEMPLATE(T, mat_entry) (mat, i, col),
                          poly->coeffs + i, ctx);

    for (; i < mat->r; i++)
        TEMPLATE(T, zero) (TEMPLATE(T, mat_entry) (mat, i, col), ctx);

}

static void
TEMPLATE(T, mat_col_to_shifted) (TEMPLATE(T, poly_t) poly,
                                 TEMPLATE(T, mat_t) mat,
                                 slong col, slong * shift,
                                 const TEMPLATE(T, ctx_t) ctx)
{
    slong i, j, rows = mat->r;

    TEMPLATE(T, poly_fit_length) (poly, rows, ctx);

    for (i = 0, j = 0; j < rows; j++)
    {
        if (shift[j])
            TEMPLATE(T, zero) (poly->coeffs + j, ctx);
        else
        {
            TEMPLATE(T, set) (poly->coeffs + j,
                              TEMPLATE(T, mat_entry) (mat, i, col), ctx);
            i++;
        }
    }

    _TEMPLATE(T, poly_set_length) (poly, rows, ctx);
    _TEMPLATE(T, poly_normalise) (poly, ctx);
}

static void
__TEMPLATE(T, poly_factor_berlekamp) (TEMPLATE(T, poly_factor_t) factors,
                                      flint_rand_t state,
                                      const TEMPLATE(T, poly_t) f,
                                      const TEMPLATE(T, ctx_t) ctx)
{
    const slong n = TEMPLATE(T, poly_degree) (f, ctx);

    TEMPLATE(T, poly_factor_t) fac1, fac2;
    TEMPLATE(T, poly_t) x, x_q;
    TEMPLATE(T, poly_t) x_qi, x_qi2;
    TEMPLATE(T, poly_t) Q, r;

    TEMPLATE(T, mat_t) matrix;
    TEMPLATE(T, t) mul, coeff, neg_one;
    fmpz_t p, q, s, pow;
    slong i, nullity, col, row;
    slong *shift;

    TEMPLATE(T, poly_t) * basis;

    if (f->length <= 2)
    {
        TEMPLATE(T, poly_factor_insert) (factors, f, 1, ctx);
        return;
    }

    TEMPLATE(T, init) (coeff, ctx);
    TEMPLATE(T, init) (neg_one, ctx);
    TEMPLATE(T, init) (mul, ctx);

    fmpz_init_set(p, TEMPLATE(T, ctx_prime) (ctx));
    fmpz_init(q);
    TEMPLATE(T, ctx_order) (q, ctx);

    TEMPLATE(T, one) (neg_one, ctx);
    TEMPLATE(T, neg) (neg_one, neg_one, ctx);


    /* s = q - 1 */
    fmpz_init_set(s, q);
    fmpz_sub_ui(s, s, 1);

    /* pow = (q-1)/2 */
    fmpz_init(pow);
    if (fmpz_cmp_ui(p, 3) > 0)
    {
        fmpz_set(pow, s);
        fmpz_divexact_ui(pow, pow, 2);
    }

    /* Step 1, compute x^q mod f in F_p[X]/<f> */
    TEMPLATE(T, poly_init) (x, ctx);
    TEMPLATE(T, poly_init) (x_q, ctx);

    TEMPLATE(T, poly_gen) (x, ctx);
    TEMPLATE(T, poly_powmod_fmpz_binexp) (x_q, x, q, f, ctx);
    TEMPLATE(T, poly_clear) (x, ctx);

    /* Step 2, compute the matrix for the Berlekamp Map */
    TEMPLATE(T, mat_init) (matrix, n, n, ctx);
    TEMPLATE(T, poly_init) (x_qi, ctx);
    TEMPLATE(T, poly_init) (x_qi2, ctx);
    TEMPLATE(T, poly_one) (x_qi, ctx);

    for (i = 0; i < n; i++)
    {
        /* Q - I */
        TEMPLATE(T, poly_set) (x_qi2, x_qi, ctx);
        TEMPLATE(T, poly_get_coeff) (coeff, x_qi2, i, ctx);
        TEMPLATE(T, sub_one) (coeff, coeff, ctx);
        TEMPLATE(T, poly_set_coeff) (x_qi2, i, coeff, ctx);
        TEMPLATE(T, to_mat_col) (matrix, i, x_qi2, ctx);
        TEMPLATE(T, poly_mulmod) (x_qi, x_qi, x_q, f, ctx);
    }

    TEMPLATE(T, poly_clear) (x_q, ctx);
    TEMPLATE(T, poly_clear) (x_qi, ctx);
    TEMPLATE(T, poly_clear) (x_qi2, ctx);

    /* Row reduce Q - I */
    nullity = n - TEMPLATE(T, mat_rref) (matrix, ctx);

    /* Find a basis for the nullspace */
    basis = flint_malloc(nullity * sizeof(TEMPLATE(T, poly_t)));
    shift = (slong *) flint_calloc(n, sizeof(slong));

    col = 1;                    /* first column is always zero */
    row = 0;
    shift[0] = 1;

    for (i = 1; i < nullity; i++)
    {
        TEMPLATE(T, poly_init) (basis[i], ctx);
        while (!TEMPLATE(T, is_zero)
               (TEMPLATE(T, mat_entry) (matrix, row, col), ctx))
        {
            row++;
            col++;
        }
        TEMPLATE(T, mat_col_to_shifted) (basis[i], matrix, col, shift, ctx);
        TEMPLATE(T, poly_set_coeff) (basis[i], col, neg_one, ctx);
        shift[col] = 1;
        col++;
    }

    flint_free(shift);
    TEMPLATE(T, mat_clear) (matrix, ctx);

    /* we are done */
    if (nullity == 1)
    {
        TEMPLATE(T, poly_factor_insert) (factors, f, 1, ctx);
    }
    else
    {
        /* Generate random linear combinations */
        TEMPLATE(T, poly_t) factor, b, power, g;
        TEMPLATE(T, poly_init) (factor, ctx);
        TEMPLATE(T, poly_init) (b, ctx);
        TEMPLATE(T, poly_init) (power, ctx);
        TEMPLATE(T, poly_init) (g, ctx);

        while (1)
        {
            do
            {
                TEMPLATE(T, poly_zero) (factor, ctx);
                for (i = 1; i < nullity; i++)
                {
                    TEMPLATE(T, randtest) (mul, state, ctx);
                    TEMPLATE(T, TEMPLATE(poly_scalar_mul, T)) (b, basis[i],
                                                               mul, ctx);
                    TEMPLATE(T, poly_add) (factor, factor, b, ctx);
                }

                TEMPLATE(T, randtest) (coeff, state, ctx);
                TEMPLATE(T, poly_set_coeff) (factor, 0, coeff, ctx);
                if (!TEMPLATE(T, poly_is_zero) (factor, ctx))
                    TEMPLATE(T, poly_make_monic) (factor, factor, ctx);
            }
            while (TEMPLATE(T, poly_is_zero) (factor, ctx) ||
                   (factor->length < 2
                    && TEMPLATE(T, is_one) (factor->coeffs, ctx)));

            TEMPLATE(T, poly_gcd) (g, f, factor, ctx);

            if (TEMPLATE(T, poly_length) (g, ctx) != 1)
                break;

            if (fmpz_cmp_ui(p, 3) > 0)
                TEMPLATE(T, poly_powmod_fmpz_binexp) (power, factor, pow, f,
                                                      ctx);
            else
                TEMPLATE(T, poly_set) (power, factor, ctx);

            TEMPLATE(T, sub_one) (power->coeffs, power->coeffs, ctx);

            _TEMPLATE(T, poly_normalise) (power, ctx);
            TEMPLATE(T, poly_gcd) (g, power, f, ctx);

            if (TEMPLATE(T, poly_length) (g, ctx) != 1)
                break;
        }

        TEMPLATE(T, poly_clear) (power, ctx);
        TEMPLATE(T, poly_clear) (factor, ctx);
        TEMPLATE(T, poly_clear) (b, ctx);

        if (!TEMPLATE(T, poly_is_zero) (g, ctx))
            TEMPLATE(T, poly_make_monic) (g, g, ctx);

        TEMPLATE(T, poly_factor_init) (fac1, ctx);
        TEMPLATE(T, poly_factor_init) (fac2, ctx);
        __TEMPLATE(T, poly_factor_berlekamp) (fac1, state, g, ctx);
        TEMPLATE(T, poly_init) (Q, ctx);
        TEMPLATE(T, poly_init) (r, ctx);
        TEMPLATE(T, poly_divrem) (Q, r, f, g, ctx);
        TEMPLATE(T, poly_clear) (r, ctx);

        if (!TEMPLATE(T, poly_is_zero) (Q, ctx))
            TEMPLATE(T, poly_make_monic) (Q, Q, ctx);

        __TEMPLATE(T, poly_factor_berlekamp) (fac2, state, Q, ctx);
        TEMPLATE(T, poly_factor_concat) (factors, fac1, ctx);
        TEMPLATE(T, poly_factor_concat) (factors, fac2, ctx);
        TEMPLATE(T, poly_factor_clear) (fac1, ctx);
        TEMPLATE(T, poly_factor_clear) (fac2, ctx);
        TEMPLATE(T, poly_clear) (Q, ctx);
        TEMPLATE(T, poly_clear) (g, ctx);
    }

    for (i = 1; i < nullity; i++)
        TEMPLATE(T, poly_clear) (basis[i], ctx);
    flint_free(basis);

    TEMPLATE(T, clear) (coeff, ctx);
    TEMPLATE(T, clear) (neg_one, ctx);
    TEMPLATE(T, clear) (mul, ctx);
    fmpz_clear(pow);
    fmpz_clear(p);
    fmpz_clear(q);
    fmpz_clear(s);
}

void
TEMPLATE(T, poly_factor_berlekamp) (TEMPLATE(T, poly_factor_t) factors,
                                    const TEMPLATE(T, poly_t) f,
                                    const TEMPLATE(T, ctx_t) ctx)
{
    slong i;
    flint_rand_t r;
    TEMPLATE(T, poly_t) v;
    TEMPLATE(T, poly_factor_t) sq_free;

    TEMPLATE(T, poly_init) (v, ctx);

    TEMPLATE(T, poly_make_monic) (v, f, ctx);

    /* compute squarefree factorisation */
    TEMPLATE(T, poly_factor_init) (sq_free, ctx);
    TEMPLATE(T, poly_factor_squarefree) (sq_free, v, ctx);

    /* run Berlekamp algorithm for all squarefree factors */
    flint_randinit(r);
    for (i = 0; i < sq_free->num; i++)
    {
        __TEMPLATE(T, poly_factor_berlekamp) (factors, r, sq_free->poly + i,
                                              ctx);
    }
    flint_randclear(r);

    /* compute multiplicities of factors in f */
    for (i = 0; i < factors->num; i++)
        factors->exp[i] = TEMPLATE(T, poly_remove) (v, factors->poly + i, ctx);

    TEMPLATE(T, poly_clear) (v, ctx);
    TEMPLATE(T, poly_factor_clear) (sq_free, ctx);
}


#endif
