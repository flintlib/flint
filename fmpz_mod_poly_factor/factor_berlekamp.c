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

    Copyright (C) 2007 David Howden
    Copyright (C) 2007, 2008, 2009, 2010 William Hart
    Copyright (C) 2008 Richard Howell-Peak
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2012 Lina Kulakova

******************************************************************************/

#include <stdlib.h>
#include "fmpz_mod_poly_factor.h"
#include "fmpz_mat.h"
#include "ulong_extras.h"
#include "profiler.h"
#include "perm.h"

static void
fmpz_mod_poly_to_fmpz_mat_col(fmpz_mat_t mat, len_t col, fmpz_mod_poly_t poly)
{
    len_t i;

    for (i = 0; i < poly->length; i++)
        fmpz_set(fmpz_mat_entry(mat, i, col), poly->coeffs + i);

    for (; i < mat->r; i++)
        fmpz_zero(fmpz_mat_entry(mat, i, col));

}

static void
fmpz_mat_col_to_fmpz_mod_poly_shifted(fmpz_mod_poly_t poly, fmpz_mat_t mat,
                                      len_t col, len_t *shift)
{
    len_t i, j, rows = mat->r;

    fmpz_mod_poly_fit_length(poly, rows);

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
                                 flint_rand_t state, const fmpz_mod_poly_t f)
{
    const len_t n = fmpz_mod_poly_degree(f);

    fmpz_mod_poly_factor_t fac1, fac2;
    fmpz_mod_poly_t x, x_p;
    fmpz_mod_poly_t x_pi, x_pi2;
    fmpz_mod_poly_t Q, r;
    fmpz_mat_t matrix;
    fmpz_t coeff, p, q, mul, pow;
    len_t i, nullity, col, row;
    len_t *shift, *perm;
    fmpz_mod_poly_t *basis;

    if (f->length <= 2)
    {
        fmpz_mod_poly_factor_insert(factors, f, 1);
        return;
    }

    fmpz_init(coeff);
    fmpz_init(mul);

    fmpz_init_set(p, &f->p);

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
    fmpz_mod_poly_init(x, p);
    fmpz_mod_poly_init(x_p, p);

    fmpz_mod_poly_set_coeff_ui(x, 1, 1);
    fmpz_mod_poly_powmod_fmpz_binexp(x_p, x, p, f);
    fmpz_mod_poly_clear(x);

    /* Step 2, compute the matrix for the Berlekamp Map */
    fmpz_mat_init(matrix, n, n);
    fmpz_mod_poly_init(x_pi, p);
    fmpz_mod_poly_init(x_pi2, p);
    fmpz_mod_poly_set_coeff_ui(x_pi, 0, 1);

    for (i = 0; i < n; i++)
    {
        /* Q - I */
        fmpz_mod_poly_set(x_pi2, x_pi);
        fmpz_mod_poly_get_coeff_fmpz(coeff, x_pi2, i);
        if (!fmpz_is_zero(coeff))
        {
            fmpz_sub_ui(coeff, coeff, 1);
            fmpz_mod(coeff, coeff, p);
            fmpz_mod_poly_set_coeff_fmpz(x_pi2, i, coeff);
        }
        else
        {
            fmpz_mod_poly_set_coeff_fmpz(x_pi2, i, q);
        }
        fmpz_mod_poly_to_fmpz_mat_col(matrix, i, x_pi2);
        fmpz_mod_poly_mulmod(x_pi, x_pi, x_p, f);
    }

    fmpz_mod_poly_clear(x_p);
    fmpz_mod_poly_clear(x_pi);
    fmpz_mod_poly_clear(x_pi2);

    /* Row reduce Q - I */
    perm = _perm_init(n);
    nullity = n - fmpz_mat_rref_mod(perm, matrix, p);
    _perm_clear(perm);

    /* Find a basis for the nullspace */
    basis =
        (fmpz_mod_poly_t *) flint_malloc(nullity * sizeof(fmpz_mod_poly_t));
    shift = (len_t *) flint_calloc(n, sizeof(len_t));

    col = 1;                    /* first column is always zero */
    row = 0;
    shift[0] = 1;

    for (i = 1; i < nullity; i++)
    {
        fmpz_mod_poly_init(basis[i], p);
        while (!fmpz_is_zero(fmpz_mat_entry(matrix, row, col)))
        {
            row++;
            col++;
        }
        fmpz_mat_col_to_fmpz_mod_poly_shifted(basis[i], matrix, col, shift);
        fmpz_mod_poly_set_coeff_fmpz(basis[i], col, q);
        shift[col] = 1;
        col++;
    }

    flint_free(shift);
    fmpz_mat_clear(matrix);

    /* we are done */
    if (nullity == 1)
    {
        fmpz_mod_poly_factor_insert(factors, f, 1);
    }
    else
    {
        /* Generate random linear combinations */
        fmpz_mod_poly_t factor, b, power, g;
        fmpz_mod_poly_init(factor, p);
        fmpz_mod_poly_init(b, p);
        fmpz_mod_poly_init(power, p);
        fmpz_mod_poly_init(g, p);

        while (1)
        {
            do
            {
                fmpz_mod_poly_zero(factor);
                for (i = 1; i < nullity; i++)
                {
                    fmpz_randm(mul, state, p);
                    fmpz_mod_poly_scalar_mul_fmpz(b, basis[i], mul);
                    fmpz_mod_poly_add(factor, factor, b);
                }

                fmpz_randm(coeff, state, p);
                fmpz_mod_poly_set_coeff_fmpz(factor, 0, coeff);
                if (!fmpz_mod_poly_is_zero(factor))
                    fmpz_mod_poly_make_monic(factor, factor);
            }
            while (fmpz_mod_poly_is_zero(factor) ||
                   (factor->length < 2 && fmpz_is_one(factor->coeffs)));

            fmpz_mod_poly_gcd(g, f, factor);

            if (fmpz_mod_poly_length(g) != 1)
                break;

            if (fmpz_cmp_ui(p, 3) > 0)
                fmpz_mod_poly_powmod_fmpz_binexp(power, factor, pow, f);
            else
                fmpz_mod_poly_set(power, factor);

            fmpz_add(power->coeffs, power->coeffs, q);
            fmpz_mod(power->coeffs, power->coeffs, p);
            _fmpz_mod_poly_normalise(power);
            fmpz_mod_poly_gcd(g, power, f);

            if (fmpz_mod_poly_length(g) != 1)
                break;
        }

        fmpz_mod_poly_clear(power);
        fmpz_mod_poly_clear(factor);
        fmpz_mod_poly_clear(b);

        if (!fmpz_mod_poly_is_zero(g))
            fmpz_mod_poly_make_monic(g, g);

        fmpz_mod_poly_factor_init(fac1);
        fmpz_mod_poly_factor_init(fac2);
        __fmpz_mod_poly_factor_berlekamp(fac1, state, g);
        fmpz_mod_poly_init(Q, p);
        fmpz_mod_poly_init(r, p);
        fmpz_mod_poly_divrem(Q, r, f, g);
        fmpz_mod_poly_clear(r);

        if (!fmpz_mod_poly_is_zero(Q))
            fmpz_mod_poly_make_monic(Q, Q);

        __fmpz_mod_poly_factor_berlekamp(fac2, state, Q);
        fmpz_mod_poly_factor_concat(factors, fac1);
        fmpz_mod_poly_factor_concat(factors, fac2);
        fmpz_mod_poly_factor_clear(fac1);
        fmpz_mod_poly_factor_clear(fac2);
        fmpz_mod_poly_clear(Q);
        fmpz_mod_poly_clear(g);
    }

    for (i = 1; i < nullity; i++)
        fmpz_mod_poly_clear(basis[i]);
    flint_free(basis);

    fmpz_clear(coeff);
    fmpz_clear(p);
    fmpz_clear(q);
    fmpz_clear(mul);
    fmpz_clear(pow);
}

void
fmpz_mod_poly_factor_berlekamp(fmpz_mod_poly_factor_t factors,
                               const fmpz_mod_poly_t f)
{
    len_t i;
    flint_rand_t r;
    fmpz_mod_poly_t v;
    fmpz_mod_poly_factor_t sq_free;

    fmpz_mod_poly_init(v, &f->p);

    fmpz_mod_poly_make_monic(v, f);

    /* compute squarefree factorisation */
    fmpz_mod_poly_factor_init(sq_free);
    fmpz_mod_poly_factor_squarefree(sq_free, v);

    /* run Berlekamp algorithm for all squarefree factors */
    flint_randinit(r);
    for (i = 0; i < sq_free->num; i++)
    {
        __fmpz_mod_poly_factor_berlekamp(factors, r, sq_free->poly + i);
    }
    flint_randclear(r);

    /* compute multiplicities of factors in f */
    for (i = 0; i < factors->num; i++)
        factors->exp[i] = fmpz_mod_poly_remove(v, factors->poly + i);

    fmpz_mod_poly_clear(v);
    fmpz_mod_poly_factor_clear(sq_free);
}
