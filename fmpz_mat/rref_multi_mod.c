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
    fmpz_t bound, b, c, prod, one, num, d, u, t;
    fmpq_mat_t C;
    nmod_mat_t Amod;
    mp_limb_t p;
    slong i, j, m, n, rank, new_rank;
    slong * pivots, * new_pivots, * P;

    m = fmpz_mat_nrows(A);
    n = fmpz_mat_ncols(A);

    fmpz_init(bound);
    fmpz_init(b);
    fmpz_init(c);
    fmpz_init(prod);
    fmpz_init(one);

    nmod_mat_init(Amod, m, n, 2);
    fmpq_mat_init(C, m, n);
    pivots = (slong *) flint_malloc(n * sizeof(slong));
    new_pivots = (slong *) flint_malloc(n * sizeof(slong));

    fmpz_one(one);

    /* compute bound */
    for (i = 0; i < m; i++)
        for (j = 0; j < n; j++)
            if (fmpz_cmpabs(b, fmpz_mat_entry(A, i, j)) < 0)
                fmpz_abs(b, fmpz_mat_entry(A, i, j));
    fmpz_set(c, b);
    fmpz_mul_2exp(c, c, m);
    fmpz_mul(bound, b, c);
    fmpz_mul_si(bound, bound, n);
    fmpz_add_ui(bound, bound, UWORD(1));

    p = UWORD(1) << NMOD_MAT_OPTIMAL_MODULUS_BITS;
    rank = -1;
    while (1)
    {
        int success = 1;
        fmpz_one(prod);

        while (fmpz_cmp(prod, bound) <= 0)
        {
            p = n_nextprime(p, 0);

            _nmod_mat_set_mod(Amod, p);
            fmpz_mat_get_nmod_mat(Amod, A);

            P = _perm_init(m);
            new_rank = _nmod_mat_rref(Amod, new_pivots, P);
            _perm_clear(P);

            /* new rank and column profile is better */
            if (compare_pivots(pivots, rank, new_pivots, new_rank) < 0)
            {
                rank = new_rank;
                for (i = 0; i < n; i++)
                    pivots[i] = new_pivots[i];

                fmpz_mat_set_nmod_mat(R, Amod);
                fmpz_set_ui(prod, p);
            }
            else
            {
                fmpz_mat_CRT_ui(R, R, prod, Amod, 1);
                fmpz_mul_ui(prod, prod, p);
            }
        }

        /* fill in the nonpivot columns */
        fmpz_init(num);
        fmpz_init(t);
        fmpz_init(d);
        fmpz_init(u);

        fmpz_one(d);

        for (i = 0; i < m && success; i++)
        {
            for (j = 0; j < n - rank; j++)
            {
                fmpz_mul(t, d, fmpz_mat_entry(R, i, pivots[rank + j]));
                fmpz_fdiv_qr(u, t, t, prod);

                success = _fmpq_reconstruct_fmpz(num, den, t, prod);

                fmpz_mul(den, den, d);
                fmpz_set(d, den);

                if (!success)
                    break;

                fmpz_set(fmpq_mat_entry_num(C, i, pivots[rank + j]), num);
                fmpz_set(fmpq_mat_entry_den(C, i, pivots[rank + j]), den);
                fmpq_canonicalise(fmpq_mat_entry(C, i, pivots[rank + j]));
            }
        }

        fmpz_clear(num);
        fmpz_clear(d);
        fmpz_clear(u);
        fmpz_clear(t);

        if (!success)
        {
            fmpz_mul2_uiui(bound, bound, p, p);
            fmpz_mul2_uiui(bound, bound, p, p);
            fmpz_mul(bound, bound, b);
            continue;
        }

        /* fill in the pivot columns */
        for (i = 0; i < m; i++)
        {
            for (j = 0; j < rank; j++)
            {
                if (i == j)
                    fmpq_one(fmpq_mat_entry(C, i, pivots[j]));
                else
                    fmpq_zero(fmpq_mat_entry(C, i, pivots[j]));
            }
        }

        fmpz_one(den);
        for (i = 0; i < m; i++)
            for (j = 0; j < n; j++)
                fmpz_lcm(den, den, fmpq_mat_entry_den(C, i, j));
        /* multiply up */
        for (i = 0; i < m; i++)
        {
            for (j = 0; j < n; j++)
            {
                fmpz_mul(fmpz_mat_entry(R, i, j), den,
                        fmpq_mat_entry_num(C, i, j));
                fmpz_divexact(fmpz_mat_entry(R, i, j), fmpz_mat_entry(R, i, j),
                        fmpq_mat_entry_den(C, i, j));
                if (fmpz_cmpabs(c, fmpq_mat_entry_num(C, i, j)) < 0)
                    fmpz_abs(c, fmpq_mat_entry_num(C, i, j));
                if (fmpz_cmpabs(c, fmpq_mat_entry_den(C, i, j)) < 0)
                    fmpz_abs(c, fmpq_mat_entry_den(C, i, j));
            }
        }

        fmpz_mul(bound, b, c);
        fmpz_mul_si(bound, bound, n);

        if (fmpz_cmp(prod, bound) > 0) break;
        fmpz_mul2_uiui(bound, bound, p, p);
        fmpz_mul2_uiui(bound, bound, p, p);
    }

    /* cleanup */
    flint_free(pivots);
    flint_free(new_pivots);
    nmod_mat_clear(Amod);
    fmpq_mat_clear(C);

    fmpz_clear(b);
    fmpz_clear(c);
    fmpz_clear(prod);
    fmpz_clear(bound);
    fmpz_clear(one);

    return rank;
}

