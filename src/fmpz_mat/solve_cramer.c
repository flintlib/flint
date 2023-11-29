/*
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mat.h"

#define AA(i,j) fmpz_mat_entry(A, i, j)
#define BB(i,j) fmpz_mat_entry(B, i, j)
#define XX(i,j) fmpz_mat_entry(X, i, j)

int
_fmpz_mat_solve_cramer_3x3(fmpz_mat_t X, fmpz_t den,
    const fmpz_mat_t A, const fmpz_mat_t B)
{
    fmpz_t t15, t16, t17;
    int success;

    fmpz_init(t15);
    fmpz_init(t16);
    fmpz_init(t17);

    fmpz_fmms(t17, AA(1,0), AA(2,1), AA(1,1), AA(2,0));
    fmpz_fmms(t16, AA(1,2), AA(2,0), AA(1,0), AA(2,2));
    fmpz_fmms(t15, AA(1,1), AA(2,2), AA(1,2), AA(2,1));

    fmpz_mul   (den, t15, AA(0,0));
    fmpz_addmul(den, t16, AA(0,1));
    fmpz_addmul(den, t17, AA(0,2));

    success = !fmpz_is_zero(den);

    if (success)
    {
        fmpz_t t12, t13, t14, x0, x1, x2;
        slong i, n = fmpz_mat_ncols(B);

        fmpz_init(t12);
        fmpz_init(t13);
        fmpz_init(t14);
        fmpz_init(x0);
        fmpz_init(x1);
        fmpz_init(x2);

        for (i = 0; i < n; i++)
        {
            fmpz_fmms(t14, AA(2,0), BB(1,i), AA(1,0), BB(2,i));
            fmpz_fmms(t13, AA(2,1), BB(1,i), AA(1,1), BB(2,i));
            fmpz_fmms(t12, AA(2,2), BB(1,i), AA(1,2), BB(2,i));

            fmpz_mul   (x0, t15, BB(0,i));
            fmpz_addmul(x0, t13, AA(0,2));
            fmpz_submul(x0, t12, AA(0,1));

            fmpz_mul   (x1, t16, BB(0,i));
            fmpz_addmul(x1, t12, AA(0,0));
            fmpz_submul(x1, t14, AA(0,2));

            fmpz_mul   (x2, t17, BB(0,i));
            fmpz_addmul(x2, t14, AA(0,1));
            fmpz_submul(x2, t13, AA(0,0));

            fmpz_swap(XX(0,i), x0);
            fmpz_swap(XX(1,i), x1);
            fmpz_swap(XX(2,i), x2);
        }

        fmpz_clear(t12);
        fmpz_clear(t13);
        fmpz_clear(t14);
        fmpz_clear(x0);
        fmpz_clear(x1);
        fmpz_clear(x2);
    }

    fmpz_clear(t15);
    fmpz_clear(t16);
    fmpz_clear(t17);

    return success;
}

int
fmpz_mat_solve_cramer(fmpz_mat_t X, fmpz_t den,
                            const fmpz_mat_t A, const fmpz_mat_t B)
{
    slong i, dim = fmpz_mat_nrows(A);

    if (dim == 0 || fmpz_mat_ncols(B) == 0)
    {
        fmpz_one(den);
        return 1;
    }
    else if (dim == 1)
    {
        fmpz_set(den, fmpz_mat_entry(A, 0, 0));

        if (fmpz_is_zero(den))
            return 0;

        if (!fmpz_mat_is_empty(B))
            _fmpz_vec_set(X->rows[0], B->rows[0], fmpz_mat_ncols(B));
        return 1;
    }
    else if (dim == 2)
    {
        fmpz_t t, u;

        fmpz_fmms(den, fmpz_mat_entry(A, 0, 0), fmpz_mat_entry(A, 1, 1),
                       fmpz_mat_entry(A, 0, 1), fmpz_mat_entry(A, 1, 0));

        if (fmpz_is_zero(den))
            return 0;

        fmpz_init(t);
        fmpz_init(u);

        for (i = 0; i < fmpz_mat_ncols(B); i++)
        {
            fmpz_fmms(t, fmpz_mat_entry(A, 1, 1), fmpz_mat_entry(B, 0, i),
                         fmpz_mat_entry(A, 0, 1), fmpz_mat_entry(B, 1, i));
            fmpz_fmms(u, fmpz_mat_entry(A, 0, 0), fmpz_mat_entry(B, 1, i),
                         fmpz_mat_entry(A, 1, 0), fmpz_mat_entry(B, 0, i));
            fmpz_swap(fmpz_mat_entry(X, 0, i), t);
            fmpz_swap(fmpz_mat_entry(X, 1, i), u);
        }

        fmpz_clear(t);
        fmpz_clear(u);

        return 1;
    }
    else if (dim == 3)
    {
        if (X == A)
        {
            int success;
            fmpz_mat_t T;
            fmpz_mat_init(T, 3, 3);
            success = _fmpz_mat_solve_cramer_3x3(T, den, A, B);
            fmpz_mat_swap_entrywise(T, X);
            fmpz_mat_clear(T);
            return success;
        }
        else
        {
            return _fmpz_mat_solve_cramer_3x3(X, den, A, B);
        }
    }
    else
    {
        flint_throw(FLINT_ERROR, "Exception (fmpz_mat_solve_cramer). dim > 3 not implemented.");
    }
}
