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

******************************************************************************/

#include <stdlib.h>
#include "nmod_poly.h"
#include "nmod_mat.h"
#include "ulong_extras.h"
#include "profiler.h"
#include "perm.h"

static void
nmod_poly_to_nmod_mat_col(nmod_mat_t mat, len_t col, nmod_poly_t poly)
{
    len_t i;

    for (i = 0; i < poly->length; i++)
        nmod_mat_entry(mat, i, col) = poly->coeffs[i];

    for ( ; i < mat->r; i++)
        nmod_mat_entry(mat, i, col) = 0UL;
}

static void
nmod_mat_col_to_nmod_poly_shifted(nmod_poly_t poly, nmod_mat_t mat,
                                            len_t col, len_t * shift)
{
    len_t i, j, rows = mat->r;

    nmod_poly_fit_length(poly, rows);

    for (i = 0, j = 0; j < rows; j++)
    {
        if (shift[j])
            poly->coeffs[j] = 0;
        else
        {
            poly->coeffs[j] = nmod_mat_entry(mat, i, col);
            i++;
        }
    }

    poly->length = rows;
    _nmod_poly_normalise(poly);
}

static void
__nmod_poly_factor_berlekamp(nmod_poly_factor_t factors,
    flint_rand_t state, const nmod_poly_t f)
{
    const mp_limb_t p = nmod_poly_modulus(f);
    const len_t n      = nmod_poly_degree(f);

    nmod_poly_factor_t fac1, fac2;
    nmod_poly_t x, x_p;
    nmod_poly_t x_pi, x_pi2;
    nmod_poly_t Q;
    nmod_mat_t matrix;
    mp_limb_t coeff;
    len_t i, nullity, col, row, *shift;
    nmod_poly_t *basis;

    if (f->length <= 2)
    {
        nmod_poly_factor_insert(factors, f, 1);
        return;
    }

    /* Step 1, we compute x^p mod f in F_p[X]/<f> */
    nmod_poly_init(x, p);
    nmod_poly_init(x_p, p);

    nmod_poly_set_coeff_ui(x, 1, 1);
    nmod_poly_powmod_ui_binexp(x_p, x, p, f);
    nmod_poly_clear(x);

    /* Step 2, compute the matrix for the Berlekamp Map */
    nmod_mat_init(matrix, n, n, p);
    nmod_poly_init(x_pi, p);
    nmod_poly_init(x_pi2, p);
    nmod_poly_set_coeff_ui(x_pi, 0, 1);

    for (i = 0; i < n; i++)
    {
        /* Q - I */
        nmod_poly_set(x_pi2, x_pi);
        coeff = nmod_poly_get_coeff_ui(x_pi2, i);
        if (coeff)
            nmod_poly_set_coeff_ui(x_pi2, i, coeff - 1);
        else
            nmod_poly_set_coeff_ui(x_pi2, i, p - 1);
        nmod_poly_to_nmod_mat_col(matrix, i, x_pi2);
        nmod_poly_mulmod(x_pi, x_pi, x_p, f);
    }

    nmod_poly_clear(x_p);
    nmod_poly_clear(x_pi);
    nmod_poly_clear(x_pi2);

    /* Row reduce Q - I */
    nullity = n - nmod_mat_rref(matrix);

    /* Find a basis for the nullspace */
    basis = (nmod_poly_t *) flint_malloc(nullity * sizeof(nmod_poly_t));
    shift = (len_t *) flint_calloc(n, sizeof(len_t));

    col = 1;  /* first column is always zero */
    row = 0;
    shift[0] = 1;

    for (i = 1; i < nullity; i++)
    {
        nmod_poly_init(basis[i], p);
        while (nmod_mat_entry(matrix, row, col))
        {
            row++;
            col++;
        }
        nmod_mat_col_to_nmod_poly_shifted(basis[i], matrix, col, shift);
        nmod_poly_set_coeff_ui(basis[i], col, p - 1);
        shift[col] = 1;
        col++;
    }

    flint_free(shift);
    nmod_mat_clear(matrix);

    /* we are done */
    if (nullity == 1)
    {
        nmod_poly_factor_insert(factors, f, 1);
        flint_free(basis);
    }
    else
    {
        /* Generate random linear combinations */
        nmod_poly_t factor, b, power, g;
        nmod_poly_init(factor, p);
        nmod_poly_init(b, p);
        nmod_poly_init(power, p);
        nmod_poly_init(g, p);

        while (1)
        {
            do
            {
                nmod_poly_zero(factor);
                for (i = 1; i < nullity; i++)
                {
                    nmod_poly_scalar_mul_nmod(b, basis[i], n_randint(state, p));
                    nmod_poly_add(factor, factor, b);
                }

                nmod_poly_set_coeff_ui(factor, 0, n_randint(state, p));
                if (!nmod_poly_is_zero(factor))
                    nmod_poly_make_monic(factor, factor);
            }
            while (nmod_poly_is_one(factor) || nmod_poly_is_zero(factor));

            nmod_poly_gcd(g, f, factor);

            if (nmod_poly_length(g) != 1) break;

            if (p > 3)
                nmod_poly_powmod_ui_binexp(power, factor, p >> 1, f);
            else
                nmod_poly_set(power, factor);

            power->coeffs[0] = n_addmod(power->coeffs[0], p - 1, p);
            _nmod_poly_normalise(power);
            nmod_poly_gcd(g, power, f);

            if (nmod_poly_length(g) != 1)
                break;
        }

        for (i = 1; i < nullity; i++)
            nmod_poly_clear(basis[i]);

        flint_free(basis);
        nmod_poly_clear(power);
        nmod_poly_clear(factor);
        nmod_poly_clear(b);

        if (!nmod_poly_is_zero(g))
            nmod_poly_make_monic(g, g);

        nmod_poly_factor_init(fac1);
        nmod_poly_factor_init(fac2);
        __nmod_poly_factor_berlekamp(fac1, state, g);
        nmod_poly_init(Q, p);
        nmod_poly_div(Q, f, g);

        if (!nmod_poly_is_zero(Q))
            nmod_poly_make_monic(Q, Q);

        __nmod_poly_factor_berlekamp(fac2, state, Q);
        nmod_poly_factor_concat(factors, fac1);
        nmod_poly_factor_concat(factors, fac2);
        nmod_poly_factor_clear(fac1);
        nmod_poly_factor_clear(fac2);
        nmod_poly_clear(Q);
        nmod_poly_clear(g);
    }
}

void
nmod_poly_factor_berlekamp(nmod_poly_factor_t factors, const nmod_poly_t f)
{
    flint_rand_t r;
    flint_randinit(r);
    __nmod_poly_factor_berlekamp(factors, r, f);
    flint_randclear(r);
}
