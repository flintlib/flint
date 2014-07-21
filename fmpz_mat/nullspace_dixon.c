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

int
_fmpz_mat_nullspace_element(fmpz_mat_t x, const fmpz_mat_t mat,
        ulong randbits)
{
    int success;
    /* TODO work out/tune for best initial prec */
    slong i, j, k, n, m, rank, nullity, prec = 4;
    slong * pivnonpiv;
    slong * pivots;
    slong * nonpivots;
    mp_limb_t p;
    nmod_mat_t matmod, bmod, T, Tb;
    fmpz_mat_t y, maty, b;
    fmpq_mat_t x_q;
    fmpz_t den, M, t;
    flint_rand_t state;

    success = 1;
    m = mat->r;
    n = mat->c;
    flint_randinit(state);
    fmpz_init(M);
    fmpz_init(den);
    fmpz_init(t);
    nmod_mat_init(matmod, m, n + m, 2);
    nmod_mat_init(bmod, m, 1, 2);
    nmod_mat_init(Tb, m, 1, 2);
    fmpq_mat_init(x_q, n, 1);
    fmpz_mat_init(y, n, 1);
    fmpz_mat_init(b, m, 1);
    fmpz_mat_init(maty, m, 1);
    pivnonpiv = flint_malloc(sizeof(slong) * FLINT_MAX(m, n));

    while (1)
    {
        p = n_randprime(state, randbits, 1);
        _nmod_mat_set_mod(matmod, p);
        _nmod_mat_set_mod(Tb, p);
        _nmod_mat_set_mod(bmod, p);
        /* set matmod to be mat mod p with an identity matrix on the right */
        for (i = 0; i < m; i++)
        {
            for (j = 0; j < n; j++)
                matmod->rows[i][j] = fmpz_fdiv_ui(fmpz_mat_entry(mat, i, j), p);
            for (j = n; j < n + m; j++)
            {
                if (j == n + i)
                    matmod->rows[i][j] = UWORD(1);
                else
                    matmod->rows[i][j] = UWORD(0);
            }
        }

        nmod_mat_rref(matmod);

        rank = 0;
        for (i = j = 0; i < m; i++, j++)
        {
            for (; j < n && nmod_mat_entry(matmod, i, j) == UWORD(0); j++);
            if (j == n)
                break;
            rank = i + 1;
        }
        nullity = n - rank;

        if (nullity == 0)
        {
            success = 0;
            break;
        }

        pivots = pivnonpiv;            /* length = rank */
        nonpivots = pivnonpiv + rank;  /* length = nullity */

        /* find pivots and nonpivots of rref mod p */
        for (i = j = k = 0; i < rank; i++)
        {
            while (nmod_mat_entry(matmod, i, j) == UWORD(0))
            {
                nonpivots[k] = j;
                k++;
                j++;
            }
            pivots[i] = j;
            j++;
        }
        while (k < nullity)
        {
            nonpivots[k] = j;
            k++;
            j++;
        }

        nmod_mat_window_init(T, matmod, 0, n, m, n + m);
        fmpz_one(M);
        fmpz_mat_zero(x);
        fmpz_mat_zero(b);

        for (k = 0; k < prec; k++)
        {
            fmpz_mat_get_nmod_mat(bmod, b);
            /* solve Ay = b (mod p) */
            nmod_mat_mul(Tb, T, bmod);
            /* for (i = 0; i < nullity; i++) */
            for (j = 0; j < rank; j++)
            {
                mp_limb_t c = nmod_mat_entry(matmod, j, nonpivots[0]);
                fmpz_set_ui(fmpz_mat_entry(y, pivots[j], 0),
                        nmod_sub(nmod_mat_entry(Tb, j, 0),
                            c, matmod->mod));
            }

            nmod_mat_entry(y, nonpivots[0], 0) = UWORD(1);

            for (j = 1; j < nullity; j++)
                nmod_mat_entry(y, nonpivots[j], 0) = UWORD(0);

            /* update b */
            fmpz_mat_mul(maty, mat, y);
            fmpz_mat_sub(b, b, maty);
            fmpz_mat_scalar_divexact_ui(b, b, p);

            fmpz_mat_scalar_addmul_fmpz(x, y, M);
            fmpz_mul_ui(M, M, p);
        }

        for (i = 0; i < n; i++)
            fmpz_mod(fmpz_mat_entry(x, i, 0), fmpz_mat_entry(x, i, 0), M);

        nmod_mat_window_clear(T);
        /* lift into place */
        if (fmpq_mat_set_fmpz_mat_mod_fmpz(x_q, x, M))
        {
            fmpz_one(den);
            for (i = 0; i < n; i++)
                fmpz_lcm(den, den, fmpq_mat_entry_den(x_q, i, 0));
            if (fmpz_is_one(den))
            {
                for (i = 0; i < n; i++)
                    fmpz_set(fmpz_mat_entry(x, i, 0),
                            fmpq_mat_entry_num(x_q, i, 0));
            }
            else
            {
                for (i = 0; i < n; i++)
                {
                    fmpz_divexact(t, den, fmpq_mat_entry_den(x_q, i, 0));
                    fmpz_mul(fmpz_mat_entry(x, i, 0),
                            fmpq_mat_entry_num(x_q, i, 0), t);
                }
            }
            fmpz_mat_mul(b, mat, x);
            if (fmpz_mat_is_zero(b))
                break;
        }
        prec++;
        if (randbits < NMOD_MAT_OPTIMAL_MODULUS_BITS)
            randbits++;
    }

    flint_free(pivnonpiv);
    flint_randclear(state);
    fmpz_clear(M);
    nmod_mat_clear(matmod);
    nmod_mat_clear(Tb);
    nmod_mat_clear(bmod);
    fmpz_mat_clear(b);
    fmpz_mat_clear(y);
    fmpz_mat_clear(maty);
    fmpq_mat_clear(x_q);
    fmpz_clear(t);
    fmpz_clear(den);

    return success;
}

slong
fmpz_mat_nullspace_dixon(fmpz_mat_t res, const fmpz_mat_t mat)
{
    slong i, j, n, m, nullity;
    ulong randbits = 6; /*NMOD_MAT_OPTIMAL_MODULUS_BITS;*/
    fmpz_mat_t tmp, x;

    m = mat->r;
    n = mat->c;

    fmpz_mat_init(tmp, m + n, n);
    fmpz_mat_init(x, n, 1);
    for (j = 0; j < n; j++)
    {
        for (i = 0; i < m; i++)
            fmpz_set(fmpz_mat_entry(tmp, i, j), fmpz_mat_entry(mat, i, j));
        for (i = m; i < m + n; i++)
            fmpz_zero(fmpz_mat_entry(tmp, i, j));
    }
    nullity = 0;
    while (_fmpz_mat_nullspace_element(x, tmp, randbits))
    {
        /* add to output and ensure this vector won't be found again */
        for (j = 0; j < n; j++)
        {
            fmpz_set(fmpz_mat_entry(tmp, m + nullity, j),
                    fmpz_mat_entry(x, j, 0));
            fmpz_set(fmpz_mat_entry(res, j, nullity),
                    fmpz_mat_entry(x, j, 0));
        }
        nullity++;
    }
    fmpz_mat_clear(x);
    fmpz_mat_clear(tmp);

    return nullity;
}
