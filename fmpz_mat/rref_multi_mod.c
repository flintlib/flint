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

    Copyright (C) 2014 Alex J. Best

******************************************************************************/

#include "fmpz_mat.h"
#include "fmpq_mat.h"
#include "perm.h"

/* returns a negative value if the list of pivots p1 is less optimal than the
   list p1, 0 if they are equal and a positive value otherwise */
int compare_pivots(const slong * p1, slong len1, const slong * p2, slong len2)
{
    slong i;

    if (len1 != len2)
        return (len1 - len2);
    for (i = 0; i < len1; i++)
    {
        if (p1[i] == p2[i]) continue;
        return (p2[i] - p1[i]);
    }
    return 0;
}

slong
fmpz_mat_rref_multi_mod(fmpz_mat_t R, fmpz_t den, const fmpz_mat_t A)
{
    fmpz_t bound, b, c, prod, one;
    fmpq_t tmp;
    fmpz_mat_t B;
    fmpq_mat_t C;
    nmod_mat_t Amod;
    mp_limb_t p;
    slong i, j, m, n, rank, new_rank;
    slong * pivots, * new_pivots;

    m = fmpz_mat_nrows(A);
    n = fmpz_mat_ncols(A);

    fmpz_init(bound);
    fmpz_init(b);
    fmpz_init(c);
    fmpz_init(prod);
    fmpz_init(one);
    fmpq_init(tmp);
    nmod_mat_init(Amod, m, n, 2);
    fmpq_mat_init(C, m, n);
    fmpz_mat_init(B, m, n);
    pivots = (slong *) flint_malloc(n * sizeof(slong));
    new_pivots = (slong *) flint_malloc(n * sizeof(slong));

    fmpz_one(one);
    /*TODO use R for B */

    /* compute bound */
    fmpz_one(c); /* TODO guess the height */
    for (i = 0; i < m; i++)
        for (j = 0; j < n; j++)
            if (fmpz_cmpabs(b, fmpz_mat_entry(A, i, j)) < 0)
                fmpz_set(b, fmpz_mat_entry(A, i, j));
    fmpz_mul(bound, b, c);
    fmpz_mul_si(bound, bound, n);
    fmpz_abs(bound, bound);
    fmpz_add_ui(bound, bound, UWORD(1));

    p = UWORD(1) << NMOD_MAT_OPTIMAL_MODULUS_BITS;
    rank = -1;
    while (1)
    {
        fmpz_one(prod);

        while (fmpz_cmp(prod, bound) <= 0)
        {
            p = n_nextprime(p, 0);

            _nmod_mat_set_mod(Amod, p);
            fmpz_mat_get_nmod_mat(Amod, A);
            new_rank = _nmod_mat_rref(Amod, new_pivots);
            /* new rank and column profile is better */
            if (compare_pivots(pivots, rank, new_pivots, new_rank) < 0)
            {
                rank = new_rank;
                for (i = 0; i < rank; i++)
                    pivots[i] = new_pivots[i];

                fmpz_mat_set_nmod_mat(B, Amod);
                fmpz_set_ui(prod, p);
            }
            else
            {
                fmpz_mat_CRT_ui(B, B, prod, Amod, 1);
                fmpz_mul_ui(prod, prod, p);
            }
        }

        if (!fmpq_mat_set_fmpz_mat_mod_fmpz(C, B, prod))
        {
            fmpz_mul2_uiui(bound, bound, p, p);
            fmpz_mul2_uiui(bound, bound, p, p);
            continue;
        }

        fmpz_one(den);
        for (i = 0; i < m; i++)
        {
            for (j = 0; j < n; j++)
                fmpz_lcm(den, den, fmpq_denref(fmpq_mat_entry(C, i, j)));
        }
        for (i = 0; i < m; i++)
        {
            for (j = 0; j < n; j++)
            {
                _fmpq_mul(fmpq_numref(tmp), fmpq_denref(tmp),
                        fmpq_mat_entry_num(C, i, j),
                        fmpq_mat_entry_den(C, i, j), den, one);
                if (fmpz_cmpabs(c, fmpq_numref(tmp)) < 0)
                    fmpz_abs(c, fmpq_numref(tmp));
                if (fmpz_cmpabs(c, fmpq_denref(tmp)) < 0)
                    fmpz_abs(c, fmpq_denref(tmp));
            }
        }

        fmpz_mul(bound, b, c);
        fmpz_mul_si(bound, bound, n);

        if (fmpz_cmp(prod, bound) > 0) break;
        fmpz_mul2_uiui(bound, bound, p, p);
    }

    /* cleanup */
    flint_free(pivots);
    flint_free(new_pivots);
    fmpz_mat_clear(B);
    nmod_mat_clear(Amod);
    fmpz_clear(b);
    fmpz_clear(c);
    fmpz_clear(prod);
    fmpz_clear(bound);
    fmpz_clear(one);
    fmpq_clear(tmp);

    /* multiply up */
    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
        {
            fmpz_mul(fmpz_mat_entry(R, i, j), den,
                    fmpq_numref(fmpq_mat_entry(C, i, j)));
            fmpz_divexact(fmpz_mat_entry(R, i, j), fmpz_mat_entry(R, i, j),
                    fmpq_denref(fmpq_mat_entry(C, i, j)));
        }
    }

    fmpq_mat_clear(C);

    return rank;
}

