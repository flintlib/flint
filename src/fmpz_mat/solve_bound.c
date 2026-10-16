/*
    Copyright (C) 2010, 2018, 2026 Fredrik Johansson
    Developed using Claude Fable 5

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mag.h"
#include "fmpz.h"
#include "fmpz_mat.h"

/*
    By Cramer's rule, the entries of the solution of A X = B, written over
    the common denominator det(A), are (up to sign) determinants of the
    matrix A with one column replaced by a column of B. Hadamard's bound
    for such a determinant is the product of the column norms of A with
    one factor replaced by the largest column norm of B, which is what we
    use for N; note that the row Hadamard bound, which is also valid for
    D = |det(A)|, may not be substituted here. All norms are
    accumulated in mag_t arithmetic (an upper bound with a few guard bits),
    so no multiprecision operations on huge entries are needed.
*/
void
fmpz_mat_solve_bound(fmpz_t N, fmpz_t D,
                        const fmpz_mat_t A, const fmpz_mat_t B)
{
    slong i, j, m, n;
    mag_t t, u, rbound, cbound;
    mag_struct * rnorms;

    m = A->r;
    n = B->c;

    if (m == 0)
    {
        fmpz_one(N);
        fmpz_one(D);
        return;
    }

    rnorms = _mag_vec_init(m);
    mag_init(t);
    mag_init(u);
    mag_init(rbound);
    mag_init(cbound);

    /* row norms of A, and column bound = product of column norms of A */
    mag_one(cbound);
    for (j = 0; j < m; j++)
    {
        mag_zero(u);
        for (i = 0; i < m; i++)
        {
            mag_set_fmpz(t, fmpz_mat_entry(A, i, j));
            mag_fast_addmul(rnorms + i, t, t);
            mag_fast_addmul(u, t, t);
        }
        mag_sqrt(u, u);
        mag_mul(cbound, cbound, u);
    }
    mag_one(rbound);
    for (i = 0; i < m; i++)
    {
        mag_sqrt(t, rnorms + i);
        mag_mul(rbound, rbound, t);
    }
    /* both Hadamard products bound |det(A)| */
    mag_min(t, rbound, cbound);
    mag_get_fmpz(D, t);

    /* largest column norm of B */
    mag_zero(u);
    for (j = 0; j < n; j++)
    {
        mag_zero(t);
        for (i = 0; i < m; i++)
        {
            mag_set_fmpz(rbound, fmpz_mat_entry(B, i, j));
            mag_fast_addmul(t, rbound, rbound);
        }
        mag_max(u, u, t);
    }
    mag_sqrt(u, u);
    /* Only the column bound may be used here: the numerator is a
       determinant of A with one COLUMN replaced by a column of B, and
       Hadamard applied to its columns replaces one column norm of A by
       the norm of that column of B. The row Hadamard product does not
       bound the corresponding cofactor vector (only sqrt(n) times it
       does), so min(row, col) would be unsound for N even though it is
       correct for D. */
    mag_mul(u, u, cbound);
    mag_get_fmpz(N, u);

    _mag_vec_clear(rnorms, m);
    mag_clear(t);
    mag_clear(u);
    mag_clear(rbound);
    mag_clear(cbound);
}
