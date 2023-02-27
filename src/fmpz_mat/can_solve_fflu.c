/*
    Copyright (C) 2010 Fredrik Johansson
    Copyright (C) 2021 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mat.h"
#include "perm.h"

void _fmpz_mat_window_with_perm_init(fmpz_mat_t Ap, slong * perm,
		                               const fmpz_mat_t A, slong start)
{
    slong i, n = A->r;

    Ap->entries = NULL;
    if (n > start)
        Ap->rows = (fmpz **) flint_malloc((n - start)*sizeof(fmpz *));
    else
        Ap->rows = NULL;

    for (i = 0; i < n - start; i++)
        Ap->rows[i] = A->rows[perm[start + i]];

    Ap->r = n - start;
    Ap->c = A->c;
}

void _fmpz_mat_window_with_perm_clear(fmpz_mat_t Ap)
{
    if (Ap->r != 0)
        flint_free(Ap->rows); 
}

int
fmpz_mat_can_solve_fflu(fmpz_mat_t X, fmpz_t den,
                    const fmpz_mat_t A, const fmpz_mat_t B)
{
    fmpz_mat_t LU, Ap, Bp, P1, P2;
    slong n, rank, *perm;
    int result;

    if (A->r == 0)
    {
        fmpz_mat_zero(X);
	fmpz_one(den);
        return 1;
    }

    if (A->c == 0)
    {
        fmpz_mat_zero(X);
        result = fmpz_mat_is_zero(B);
	    fmpz_set_ui(den, result);
	    return result;
    }

    n = fmpz_mat_nrows(A);

    perm = _perm_init(n);
    
    fmpz_mat_init_set(LU, A);
    rank = fmpz_mat_fflu(LU, den, perm, LU, 0);
    
    result = !fmpz_is_zero(den) && fmpz_mat_solve_fflu_precomp(X, perm, LU, B);

    if (result)
    {
        if (_perm_parity(perm, n))
        {
           fmpz_neg(den, den);

           fmpz_mat_neg(X, X);
        }

        if (rank < n)
	    {
           _fmpz_mat_window_with_perm_init(Ap, perm, A, rank);
           _fmpz_mat_window_with_perm_init(Bp, perm, B, rank);

	       fmpz_mat_init(P1, fmpz_mat_nrows(Ap), fmpz_mat_ncols(X));
           fmpz_mat_init(P2, fmpz_mat_nrows(Bp), fmpz_mat_ncols(Bp));

	       fmpz_mat_mul(P1, Ap, X);
           fmpz_mat_scalar_mul_fmpz(P2, Bp, den);

	       result = fmpz_mat_equal(P1, P2);

           fmpz_mat_clear(P1);
           fmpz_mat_clear(P2);

	       _fmpz_mat_window_with_perm_clear(Ap);
	       _fmpz_mat_window_with_perm_clear(Bp);
        }
    } else
        fmpz_zero(den);

    _perm_clear(perm);
    fmpz_mat_clear(LU);

    return result;
}
