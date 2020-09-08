/*
    Copyright (C) 2007 David Howden
    Copyright (C) 2007, 2008, 2009, 2010 William Hart
    Copyright (C) 2008 Richard Howell-Peak
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2012 Lina Kulakova

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include "fmpz_mod_poly.h"
#include "fmpz_mat.h"
#include "ulong_extras.h"
#include "profiler.h"
#include "perm.h"

static void
fmpz_mod_poly_to_fmpz_mat_col(fmpz_mat_t mat, slong col, fmpz_mod_poly_t poly)
{
    slong i;

    for (i = 0; i < poly->length; i++)
        fmpz_set(fmpz_mat_entry(mat, i, col), poly->coeffs + i);

    for (; i < mat->r; i++)
        fmpz_zero(fmpz_mat_entry(mat, i, col));

}

static void
fmpz_mat_col_to_fmpz_mod_poly_shifted(fmpz_mod_poly_t poly, fmpz_mat_t mat,
                             slong col, slong *shift, const fmpz_mod_ctx_t ctx)
{
    slong i, j, rows = mat->r;

    fmpz_mod_poly_fit_length(poly, rows, ctx);

    for (i = 0, j = 0; j < rows; j++)
    {
        if (shift[j])
            fmpz_zero(poly->coeffs + j);
        else
        {
            fmpz_set(poly->coeffs + j, fmpz_mat_entry(mat, i, col));
            i++;
        }
    }

    _fmpz_mod_poly_set_length(poly, rows);
    _fmpz_mod_poly_normalise(poly);
}

static void
__fmpz_mod_poly_factor_berlekamp(fmpz_mod_poly_factor_t factors,
         flint_rand_t state, const fmpz_mod_poly_t f, const fmpz_mod_ctx_t ctx)
{
    const slong n = fmpz_mod_poly_degree(f, ctx);
    const fmpz * p = fmpz_mod_ctx_modulus(ctx);
    fmpz_mod_poly_factor_t fac1, fac2;
    fmpz_mod_poly_t x, x_p;
    fmpz_mod_poly_t x_pi, x_pi2;
    fmpz_mod_poly_t Q, r;
    fmpz_mat_t matrix;
    fmpz_t coeff, q, mul, pow;
    slong i, nullity, col, row;
    slong *shift, *perm;
    fmpz_mod_poly_t *basis;

    if (f->length <= 2)
    {
        fmpz_mod_poly_factor_insert(factors, f, 1, ctx);
        return;
    }

    fmpz_init(coeff);
    fmpz_init(mul);

    /* q = p - 1 */
    fmpz_init_set(q, p);
    fmpz_sub_ui(q, q, 1);
    fmpz_mod(q, q, p);

    /* pow = (p-1)/2 */
    fmpz_init(pow);
    if (fmpz_cmp_ui(p, 3) > 0)
    {
        fmpz_set(pow, q);
        fmpz_divexact_ui(pow, pow, 2);
    }

    /* Step 1, compute x^p mod f in F_p[X]/<f> */
    fmpz_mod_poly_init(x, ctx);
    fmpz_mod_poly_init(x_p, ctx);

    fmpz_mod_poly_set_coeff_ui(x, 1, 1, ctx);
    fmpz_mod_poly_powmod_fmpz_binexp(x_p, x, p, f, ctx);
    fmpz_mod_poly_clear(x, ctx);

    /* Step 2, compute the matrix for the Berlekamp Map */
    fmpz_mat_init(matrix, n, n);
    fmpz_mod_poly_init(x_pi, ctx);
    fmpz_mod_poly_init(x_pi2, ctx);
    fmpz_mod_poly_set_coeff_ui(x_pi, 0, 1, ctx);

    for (i = 0; i < n; i++)
    {
        /* Q - I */
        fmpz_mod_poly_set(x_pi2, x_pi, ctx);
        fmpz_mod_poly_get_coeff_fmpz(coeff, x_pi2, i, ctx);
        if (!fmpz_is_zero(coeff))
        {
            fmpz_sub_ui(coeff, coeff, 1);
            fmpz_mod(coeff, coeff, p);
            fmpz_mod_poly_set_coeff_fmpz(x_pi2, i, coeff, ctx);
        }
        else
        {
            fmpz_mod_poly_set_coeff_fmpz(x_pi2, i, q, ctx);
        }
        fmpz_mod_poly_to_fmpz_mat_col(matrix, i, x_pi2);
        fmpz_mod_poly_mulmod(x_pi, x_pi, x_p, f, ctx);
    }

    fmpz_mod_poly_clear(x_p, ctx);
    fmpz_mod_poly_clear(x_pi, ctx);
    fmpz_mod_poly_clear(x_pi2, ctx);

    /* Row reduce Q - I */
    perm = _perm_init(n);
    nullity = n - fmpz_mat_rref_mod(perm, matrix, p);
    _perm_clear(perm);

    /* Find a basis for the nullspace */
    basis =
        (fmpz_mod_poly_t *) flint_malloc(nullity * sizeof(fmpz_mod_poly_t));
    shift = (slong *) flint_calloc(n, sizeof(slong));

    col = 1;                    /* first column is always zero */
    row = 0;
    shift[0] = 1;

    for (i = 1; i < nullity; i++)
    {
        fmpz_mod_poly_init(basis[i], ctx);
        while (!fmpz_is_zero(fmpz_mat_entry(matrix, row, col)))
        {
            row++;
            col++;
        }
        fmpz_mat_col_to_fmpz_mod_poly_shifted(basis[i], matrix, col, shift, ctx);
        fmpz_mod_poly_set_coeff_fmpz(basis[i], col, q, ctx);
        shift[col] = 1;
        col++;
    }

    flint_free(shift);
    fmpz_mat_clear(matrix);

    /* we are done */
    if (nullity == 1)
    {
        fmpz_mod_poly_factor_insert(factors, f, 1, ctx);
    }
    else
    {
        /* Generate random linear combinations */
        fmpz_mod_poly_t factor, b, power, g;
        fmpz_mod_poly_init(factor, ctx);
        fmpz_mod_poly_init(b, ctx);
        fmpz_mod_poly_init(power, ctx);
        fmpz_mod_poly_init(g, ctx);

        while (1)
        {
            do
            {
                fmpz_mod_poly_zero(factor, ctx);
                for (i = 1; i < nullity; i++)
                {
                    fmpz_randm(mul, state, p);
                    fmpz_mod_poly_scalar_mul_fmpz(b, basis[i], mul, ctx);
                    fmpz_mod_poly_add(factor, factor, b, ctx);
                }

                fmpz_randm(coeff, state, p);
                fmpz_mod_poly_set_coeff_fmpz(factor, 0, coeff, ctx);
                if (!fmpz_mod_poly_is_zero(factor, ctx))
                    fmpz_mod_poly_make_monic(factor, factor, ctx);
            }
            while (fmpz_mod_poly_is_zero(factor, ctx) ||
                   (factor->length < 2 && fmpz_is_one(factor->coeffs)));

            fmpz_mod_poly_gcd(g, f, factor, ctx);

            if (fmpz_mod_poly_length(g, ctx) != 1)
                break;

            if (fmpz_cmp_ui(p, 3) > 0)
                fmpz_mod_poly_powmod_fmpz_binexp(power, factor, pow, f, ctx);
            else
                fmpz_mod_poly_set(power, factor, ctx);

            fmpz_add(power->coeffs, power->coeffs, q);
            fmpz_mod(power->coeffs, power->coeffs, p);
            _fmpz_mod_poly_normalise(power);
            fmpz_mod_poly_gcd(g, power, f, ctx);

            if (fmpz_mod_poly_length(g, ctx) != 1)
                break;
        }

        fmpz_mod_poly_clear(power, ctx);
        fmpz_mod_poly_clear(factor, ctx);
        fmpz_mod_poly_clear(b, ctx);

        if (!fmpz_mod_poly_is_zero(g, ctx))
            fmpz_mod_poly_make_monic(g, g, ctx);

        fmpz_mod_poly_factor_init(fac1, ctx);
        fmpz_mod_poly_factor_init(fac2, ctx);
        __fmpz_mod_poly_factor_berlekamp(fac1, state, g, ctx);
        fmpz_mod_poly_init(Q, ctx);
        fmpz_mod_poly_init(r, ctx);
        fmpz_mod_poly_divrem(Q, r, f, g, ctx);
        fmpz_mod_poly_clear(r, ctx);

        if (!fmpz_mod_poly_is_zero(Q, ctx))
            fmpz_mod_poly_make_monic(Q, Q, ctx);

        __fmpz_mod_poly_factor_berlekamp(fac2, state, Q, ctx);
        fmpz_mod_poly_factor_concat(factors, fac1, ctx);
        fmpz_mod_poly_factor_concat(factors, fac2, ctx);
        fmpz_mod_poly_factor_clear(fac1, ctx);
        fmpz_mod_poly_factor_clear(fac2, ctx);
        fmpz_mod_poly_clear(Q, ctx);
        fmpz_mod_poly_clear(g, ctx);
    }

    for (i = 1; i < nullity; i++)
        fmpz_mod_poly_clear(basis[i], ctx);
    flint_free(basis);

    fmpz_clear(coeff);
    fmpz_clear(q);
    fmpz_clear(mul);
    fmpz_clear(pow);
}

void
fmpz_mod_poly_factor_berlekamp(fmpz_mod_poly_factor_t factors,
                             const fmpz_mod_poly_t f, const fmpz_mod_ctx_t ctx)
{
    slong i;
    flint_rand_t r;
    fmpz_mod_poly_t v;
    fmpz_mod_poly_factor_t sq_free;

    fmpz_mod_poly_init(v, ctx);

    fmpz_mod_poly_make_monic(v, f, ctx);

    /* compute squarefree factorisation */
    fmpz_mod_poly_factor_init(sq_free, ctx);
    fmpz_mod_poly_factor_squarefree(sq_free, v, ctx);

    /* run Berlekamp algorithm for all squarefree factors */
    flint_randinit(r);
    for (i = 0; i < sq_free->num; i++)
    {
        __fmpz_mod_poly_factor_berlekamp(factors, r, sq_free->poly + i, ctx);
    }
    flint_randclear(r);

    /* compute multiplicities of factors in f */
    for (i = 0; i < factors->num; i++)
        factors->exp[i] = fmpz_mod_poly_remove(v, factors->poly + i, ctx);

    fmpz_mod_poly_clear(v, ctx);
    fmpz_mod_poly_factor_clear(sq_free, ctx);
}
