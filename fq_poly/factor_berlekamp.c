/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2013 Mike Hansen

******************************************************************************/

#include "fq_poly.h"
#include "fq_mat.h"
#include "perm.h"

static void
fq_to_mat_col(fq_mat_t mat, slong col, fq_poly_t poly, const fq_ctx_t ctx)
{
    slong i;

    for (i = 0; i < poly->length; i++)
        fq_set(fq_mat_entry(mat, i, col), poly->coeffs + i, ctx);

    for (; i < mat->r; i++)
        fq_zero(fq_mat_entry(mat, i, col), ctx);

}

static void
fq_mat_col_to_shifted(fq_poly_t poly, fq_mat_t mat, slong col, slong * shift,
                      const fq_ctx_t ctx)
{
    slong i, j, rows = mat->r;

    fq_poly_fit_length(poly, rows, ctx);

    for (i = 0, j = 0; j < rows; j++)
    {
        if (shift[j])
            fq_zero(poly->coeffs + j, ctx);
        else
        {
            fq_set(poly->coeffs + j, fq_mat_entry(mat, i, col), ctx);
            i++;
        }
    }

    _fq_poly_set_length(poly, rows, ctx);
    _fq_poly_normalise(poly, ctx);
}

static void
__fq_poly_factor_berlekamp(fq_poly_factor_t factors, flint_rand_t state,
                           const fq_poly_t f, const fq_ctx_t ctx)
{
    const slong n = fq_poly_degree(f, ctx);

    fq_poly_factor_t fac1, fac2;
    fq_poly_t x, x_q;
    fq_poly_t x_qi, x_qi2;
    fq_poly_t Q, r;

    fq_mat_t matrix;
    fq_t mul, coeff, neg_one;
    fmpz_t p, q, s, pow;
    slong i, nullity, col, row;
    slong *shift;

    fq_poly_t *basis;

    if (f->length <= 2)
    {
        fq_poly_factor_insert(factors, f, 1, ctx);
        return;
    }

    fq_init(coeff, ctx);
    fq_init(neg_one, ctx);
    fq_init(mul, ctx);

    fmpz_init_set(p, fq_ctx_prime(ctx));
    fmpz_init(q);
    fq_ctx_order(q, ctx);

    fq_one(neg_one, ctx);
    fq_neg(neg_one, neg_one, ctx);


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
    fq_poly_init(x, ctx);
    fq_poly_init(x_q, ctx);

    fq_poly_gen(x, ctx);
    fq_poly_powmod_fmpz_binexp(x_q, x, q, f, ctx);
    fq_poly_clear(x, ctx);

    /* Step 2, compute the matrix for the Berlekamp Map */
    fq_mat_init(matrix, n, n, ctx);
    fq_poly_init(x_qi, ctx);
    fq_poly_init(x_qi2, ctx);
    fq_poly_one(x_qi, ctx);

    for (i = 0; i < n; i++)
    {
        /* Q - I */
        fq_poly_set(x_qi2, x_qi, ctx);
        fq_poly_get_coeff(coeff, x_qi2, i, ctx);
        fq_sub_one(coeff, coeff, ctx);
        fq_poly_set_coeff(x_qi2, i, coeff, ctx);
        fq_to_mat_col(matrix, i, x_qi2, ctx);
        fq_poly_mulmod(x_qi, x_qi, x_q, f, ctx);
    }

    fq_poly_clear(x_q, ctx);
    fq_poly_clear(x_qi, ctx);
    fq_poly_clear(x_qi2, ctx);

    /* Row reduce Q - I */
    nullity = n - fq_mat_rref(matrix, ctx);

    /* Find a basis for the nullspace */
    basis = flint_malloc(nullity * sizeof(fq_poly_t));
    shift = (slong *) flint_calloc(n, sizeof(slong));

    col = 1;                    /* first column is always zero */
    row = 0;
    shift[0] = 1;

    for (i = 1; i < nullity; i++)
    {
        fq_poly_init(basis[i], ctx);
        while (!fq_is_zero(fq_mat_entry(matrix, row, col), ctx))
        {
            row++;
            col++;
        }
        fq_mat_col_to_shifted(basis[i], matrix, col, shift, ctx);
        fq_poly_set_coeff(basis[i], col, neg_one, ctx);
        shift[col] = 1;
        col++;
    }

    flint_free(shift);
    fq_mat_clear(matrix, ctx);

    /* we are done */
    if (nullity == 1)
    {
        fq_poly_factor_insert(factors, f, 1, ctx);
    }
    else
    {
        /* Generate random linear combinations */
        fq_poly_t factor, b, power, g;
        fq_poly_init(factor, ctx);
        fq_poly_init(b, ctx);
        fq_poly_init(power, ctx);
        fq_poly_init(g, ctx);

        while (1)
        {
            do
            {
                fq_poly_zero(factor, ctx);
                for (i = 1; i < nullity; i++)
                {
                    fq_randtest(mul, state, ctx);
                    fq_poly_scalar_mul_fq(b, basis[i], mul, ctx);
                    fq_poly_add(factor, factor, b, ctx);
                }

                fq_randtest(coeff, state, ctx);
                fq_poly_set_coeff(factor, 0, coeff, ctx);
                if (!fq_poly_is_zero(factor, ctx))
                    fq_poly_make_monic(factor, factor, ctx);
            }
            while (fq_poly_is_zero(factor, ctx) ||
                   (factor->length < 2 && fq_is_one(factor->coeffs, ctx)));

            fq_poly_gcd(g, f, factor, ctx);

            if (fq_poly_length(g, ctx) != 1)
                break;

            if (fmpz_cmp_ui(p, 3) > 0)
                fq_poly_powmod_fmpz_binexp(power, factor, pow, f, ctx);
            else
                fq_poly_set(power, factor, ctx);

            fq_sub_one(power->coeffs, power->coeffs, ctx);

            _fq_poly_normalise(power, ctx);
            fq_poly_gcd(g, power, f, ctx);

            if (fq_poly_length(g, ctx) != 1)
                break;
        }

        fq_poly_clear(power, ctx);
        fq_poly_clear(factor, ctx);
        fq_poly_clear(b, ctx);

        if (!fq_poly_is_zero(g, ctx))
            fq_poly_make_monic(g, g, ctx);

        fq_poly_factor_init(fac1, ctx);
        fq_poly_factor_init(fac2, ctx);
        __fq_poly_factor_berlekamp(fac1, state, g, ctx);
        fq_poly_init(Q, ctx);
        fq_poly_init(r, ctx);
        fq_poly_divrem(Q, r, f, g, ctx);
        fq_poly_clear(r, ctx);

        if (!fq_poly_is_zero(Q, ctx))
            fq_poly_make_monic(Q, Q, ctx);

        __fq_poly_factor_berlekamp(fac2, state, Q, ctx);
        fq_poly_factor_concat(factors, fac1, ctx);
        fq_poly_factor_concat(factors, fac2, ctx);
        fq_poly_factor_clear(fac1, ctx);
        fq_poly_factor_clear(fac2, ctx);
        fq_poly_clear(Q, ctx);
        fq_poly_clear(g, ctx);
    }

    for (i = 1; i < nullity; i++)
        fq_poly_clear(basis[i], ctx);
    flint_free(basis);

    fq_clear(coeff, ctx);
    fq_clear(neg_one, ctx);
    fq_clear(mul, ctx);
    fmpz_clear(pow);
    fmpz_clear(p);
    fmpz_clear(s);
}

void
fq_poly_factor_berlekamp(fq_poly_factor_t factors, const fq_poly_t f,
                         const fq_ctx_t ctx)
{
    slong i;
    flint_rand_t r;
    fq_poly_t v;
    fq_poly_factor_t sq_free;

    fq_poly_init(v, ctx);

    fq_poly_make_monic(v, f, ctx);

    /* compute squarefree factorisation */
    fq_poly_factor_init(sq_free, ctx);
    fq_poly_factor_squarefree(sq_free, v, ctx);

    /* run Berlekamp algorithm for all squarefree factors */
    flint_randinit(r);
    for (i = 0; i < sq_free->num; i++)
    {
        __fq_poly_factor_berlekamp(factors, r, sq_free->poly + i, ctx);
    }
    flint_randclear(r);

    /* compute multiplicities of factors in f */
    for (i = 0; i < factors->num; i++)
        factors->exp[i] = fq_poly_remove(v, factors->poly + i, ctx);

    fq_poly_clear(v, ctx);
    fq_poly_factor_clear(sq_free, ctx);
}
