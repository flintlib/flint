/*
    Copyright (C) 2013 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>

#include "ulong_extras.h"
#include "nmod_mat.h"
#include "nmod_poly.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "fmpz_mat.h"

#define CHARPOLY_M_LOG2E  1.44269504088896340736  /* log2(e) */

static inline long double _log2(const long double x)
{
    return log(x) * CHARPOLY_M_LOG2E;
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

void _fmpz_mat_charpoly_small(fmpz * rop, const fmpz_mat_t op)
{
    if (op->r == 0)
    {
        fmpz_one(rop + 0);
    }
    else if (op->r == 1)
    {
        fmpz_one(rop + 1);
        fmpz_neg(rop + 0, &(op->rows[0][0]));
    }
    else if (op->r == 2)
    {
        _fmpz_mat_charpoly_small_2x2(rop, op->rows);
    }
    else  /* op->r == 3 */
    {
        _fmpz_mat_charpoly_small_3x3(rop, op->rows);
    }
}

void _fmpz_mat_charpoly_modular(fmpz * rop, const fmpz_mat_t op)
{
    const slong n = op->r;

    if (n < 4)
    {
        _fmpz_mat_charpoly_small(rop, op);
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
        slong bound;

        slong pbits  = FLINT_BITS - 1;
        mp_limb_t p = (UWORD(1) << pbits);

        fmpz_t m;

        /* Determine the bound in bits */
        {
            slong i, j;
            fmpz *ptr;
            double t;

            ptr = fmpz_mat_entry(op, 0, 0);
            for (i = 0; i < n; i++)
                for (j = 0; j < n; j++)
                    if (fmpz_cmpabs(ptr, fmpz_mat_entry(op, i, j)) < 0)
                        ptr = fmpz_mat_entry(op, i, j);

            if (fmpz_bits(ptr) == 0)  /* Zero matrix */
            {
                for (i = 0; i < n; i++)
                   fmpz_zero(rop + i);
                fmpz_set_ui(rop + n, 1);
                return;
            }

            t = (fmpz_bits(ptr) <= FLINT_D_BITS) ?
                _log2(FLINT_ABS(fmpz_get_d(ptr))) : fmpz_bits(ptr);

            bound = ceil( (n / 2.0) * (_log2(n) + 2.0 * t + 1.6669) );
        }

        fmpz_init_set_ui(m, 1);

        for ( ; fmpz_bits(m) < bound; )
        {
            nmod_mat_t mat;
            nmod_poly_t poly;

            p = n_nextprime(p, 0);

            nmod_mat_init(mat, n, n, p);
            nmod_poly_init(poly, p);

            fmpz_mat_get_nmod_mat(mat, op);
            nmod_mat_charpoly(poly, mat);

            _fmpz_poly_CRT_ui(rop, rop, n + 1, m, poly->coeffs, n + 1, poly->mod.n, poly->mod.ninv, 1);

            fmpz_mul_ui(m, m, p);

            nmod_mat_clear(mat);
            nmod_poly_clear(poly);
        }

        fmpz_clear(m);
    }
}

void fmpz_mat_charpoly_modular(fmpz_poly_t cp, const fmpz_mat_t mat)
{
     fmpz_poly_fit_length(cp, mat->r + 1);
    _fmpz_poly_set_length(cp, mat->r + 1);

    _fmpz_mat_charpoly_modular(cp->coeffs, mat);
}
