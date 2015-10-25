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

    Copyright (C) 2013 Sebastian Pancratz

******************************************************************************/

#include <math.h>

#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "fmpz_mat.h"
#include "nmod_mat.h"
#include "nmod_poly.h"

#define M_LOG2E  1.44269504088896340736  /* log2(e) */

static __inline__ long double _log2(const long double x)
{
    return log(x) * M_LOG2E;
}

static void _fmpz_mat_charpoly_small_2x2(fmpz *rop, fmpz ** const x)
{
    fmpz_one   (rop + 2);
    fmpz_add   (rop + 1, &x[0][0], &x[1][1]);
    fmpz_neg   (rop + 1, rop + 1);
    fmpz_mul   (rop + 0, &x[0][0], &x[1][1]);
    fmpz_submul(rop + 0, &x[0][1], &x[1][0]);
}

static void _fmpz_mat_charpoly_small_3x3(fmpz *rop, fmpz ** const x)
{
    fmpz a[2];
    fmpz_init(a + 0);
    fmpz_init(a + 1);

    fmpz_mul(   a + 0,   &x[1][0], &x[2][1]);
    fmpz_submul(a + 0,   &x[1][1], &x[2][0]);
    fmpz_mul(   rop + 0, a + 0,    &x[0][2]);
    fmpz_neg(   rop + 0, rop + 0);
    fmpz_mul(   rop + 1, &x[2][0], &x[0][2]);
    fmpz_neg(   rop + 1, rop + 1);

    fmpz_mul(   a + 0,   &x[1][2], &x[2][0]);
    fmpz_submul(a + 0,   &x[1][0], &x[2][2]);
    fmpz_submul(rop + 0, a + 0,    &x[0][1]);
    fmpz_submul(rop + 1, &x[1][0], &x[0][1]);

    fmpz_mul(   a + 0,   &x[1][1], &x[2][2]);
    fmpz_add(   a + 1,   &x[1][1], &x[2][2]);
    fmpz_neg(   a + 1,   a + 1);
    fmpz_submul(a + 0,   &x[1][2], &x[2][1]);

    fmpz_submul(rop + 0, a + 0,    &x[0][0]);
    fmpz_submul(rop + 1, a + 1,    &x[0][0]);
    fmpz_add(   rop + 1, rop + 1,  a + 0);
    fmpz_sub(   rop + 2, a + 1,    &x[0][0]);
    fmpz_one(   rop + 3);

    fmpz_clear(a + 0);
    fmpz_clear(a + 1);
}

void fmpz_mat_charpoly_small(fmpz_poly_t rop, const fmpz_mat_t op)
{
    fmpz_poly_fit_length(rop, op->r + 1);

    if (op->r == 0)
    {
        fmpz_one(rop->coeffs + 0);
    }
    else if (op->r == 1)
    {
        fmpz_one(rop->coeffs + 1);
        fmpz_neg(rop->coeffs + 0, &(op->rows[0][0]));
    }
    else if (op->r == 2)
    {
        _fmpz_mat_charpoly_small_2x2(rop->coeffs, op->rows);
    }
    else  /* op->r == 3 */
    {
        _fmpz_mat_charpoly_small_3x3(rop->coeffs, op->rows);
    }
    _fmpz_poly_set_length(rop, op->r + 1);
}

void fmpz_mat_charpoly_modular(fmpz_poly_t rop, const fmpz_mat_t op)
{
    const long n = op->r;

    if (n < 4)
    {
        fmpz_mat_charpoly_small(rop, op);
    }
    else
    {
        /*
            If $A$ is an $n \times n$ matrix with $n \geq 4$ and 
            coefficients bounded in absolute value by $B > 1$ then 
            the coefficients of the characteristic polynomial have 
            less than $\ceil{n/2 (\log_2(n) + \log_2(B^2) + 1.6669)}$ 
            bits.
            See Lemma 4.1 in Dumas, Pernet, and Wan, "Efficient computation 
            of the characteristic polynomial", 2008.
         */
        long bound;

        long pbits  = FLINT_BITS - 1;
        mp_limb_t p = (1UL << pbits);

        fmpz_t m;

        /* Determine the bound in bits */
        {
            long i, j;
            fmpz *ptr;
            double t;

            ptr = fmpz_mat_entry(op, 0, 0);
            for (i = 0; i < n; i++)
                for (j = 0; j < n; j++)
                    if (fmpz_cmpabs(ptr, fmpz_mat_entry(op, i, j)) < 0)
                        ptr = fmpz_mat_entry(op, i, j);

            if (fmpz_bits(ptr) == 0)  /* Zero matrix */
            {
                fmpz_poly_zero(rop);
                fmpz_poly_set_coeff_ui(rop, n, 1);
                return;
            }

            t = (fmpz_bits(ptr) <= FLINT_D_BITS) ? 
                _log2(FLINT_ABS(fmpz_get_d(ptr))) : fmpz_bits(ptr);

            bound = ceil( (n / 2.0) * (_log2(n) + 2.0 * t + 1.6669) );
        }

        fmpz_init_set_ui(m, 1);
        fmpz_poly_zero(rop);

        for ( ; fmpz_bits(m) < bound; )
        {
            nmod_mat_t mat;
            nmod_poly_t poly;

            p = n_nextprime(p, 0);

            nmod_mat_init(mat, n, n, p);
            nmod_poly_init(poly, p);

            fmpz_mat_get_nmod_mat(mat, op);
            nmod_mat_charpoly(poly, mat);

            fmpz_poly_CRT_ui(rop, rop, m, poly, 1);

            fmpz_mul_ui(m, m, p);

            nmod_mat_clear(mat);
            nmod_poly_clear(poly);
        }
        fmpz_clear(m);
    }
}

