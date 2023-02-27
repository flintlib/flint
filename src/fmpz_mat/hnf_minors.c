/*
    Copyright (C) 2014 Alex J. Best
    Copyright (C) 2017 Tommy Hofmann

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mat.h"

/*
  This is the algorithm of Kannan, Bachem, "Polynomial algorithms for computing
  the Smith and Hermite normal forms of an integer matrix", Siam J. Comput.,
  Vol. 8, No. 4, pp. 499-507.
*/
void
fmpz_mat_hnf_minors(fmpz_mat_t H, const fmpz_mat_t A)
{
    slong j, j2, i, k, l, m, n;
    fmpz_t u, v, d, r2d, r1d, q, b;

    m = fmpz_mat_nrows(A);
    n = fmpz_mat_ncols(A);

    fmpz_init(u);
    fmpz_init(v);
    fmpz_init(d);
    fmpz_init(r1d);
    fmpz_init(r2d);
    fmpz_init(q);
    fmpz_init(b);
    fmpz_mat_set(H, A);

    /* put the kth principal minor in HNF */
    for (k = 0, l = m - 1; k < n; k++)
    {
        for (j = 0; j < k; j++)
        {
            if (fmpz_is_zero(fmpz_mat_entry(H, k, j)))
            {
                continue;
            }

            fmpz_xgcd(d, u, v, fmpz_mat_entry(H, j, j),
                      fmpz_mat_entry(H, k, j));

            if (fmpz_cmpabs(d, fmpz_mat_entry(H, j, j)) == 0)
            {
                 fmpz_divexact(b, fmpz_mat_entry(H, k, j), fmpz_mat_entry(H, j, j));
                 for (j2 = j; j2 < n; j2++)
                 {
                     fmpz_submul(fmpz_mat_entry(H, k, j2), b, fmpz_mat_entry(H, j, j2));
                 }
                 continue;
            }
                

            fmpz_divexact(r1d, fmpz_mat_entry(H, j, j), d);
            fmpz_divexact(r2d, fmpz_mat_entry(H, k, j), d);
            for (j2 = j; j2 < n; j2++)
            {
                fmpz_mul(b, u, fmpz_mat_entry(H, j, j2));
                fmpz_addmul(b, v, fmpz_mat_entry(H, k, j2));
                fmpz_mul(fmpz_mat_entry(H, k, j2), r1d,
                         fmpz_mat_entry(H, k, j2));
                fmpz_submul(fmpz_mat_entry(H, k, j2), r2d,
                            fmpz_mat_entry(H, j, j2));
                fmpz_set(fmpz_mat_entry(H, j, j2), b);
            }
        }
        /* if H_k,k is zero we swap row k for some other row (starting with the
           last) */
        if (fmpz_is_zero(fmpz_mat_entry(H, k, k)))
        {
            fmpz_mat_swap_rows(H, NULL, k, l);
            l--;
            k--;
            continue;
        }
        /* ensure H_k,k is positive */
        if (fmpz_sgn(fmpz_mat_entry(H, k, k)) < 0)
        {
            for (j = k; j < n; j++)
            {
                fmpz_neg(fmpz_mat_entry(H, k, j), fmpz_mat_entry(H, k, j));
            }
        }
        /* reduce above diagonal elements of each row i */
        for (i = k - 1; i >= 0; i--)
        {
            for (j = i + 1; j <= k; j++)
            {
                fmpz_fdiv_q(q, fmpz_mat_entry(H, i, j),
                            fmpz_mat_entry(H, j, j));
                if (fmpz_is_zero(q))
                {
                    continue;
                }
                for (j2 = j; j2 < n; j2++)
                {
                    fmpz_submul(fmpz_mat_entry(H, i, j2), q,
                                fmpz_mat_entry(H, j, j2));
                }
            }
        }
        l = m - 1;
    }

    /* reduce final rows */
    for (k = n; k < m; k++)
    {
        for (j = 0; j < n; j++)
        {
            fmpz_xgcd(d, u, v, fmpz_mat_entry(H, j, j),
                      fmpz_mat_entry(H, k, j));
            
            if (fmpz_cmpabs(d, fmpz_mat_entry(H, j, j)) == 0)
            {
                 fmpz_divexact(b, fmpz_mat_entry(H, k, j), fmpz_mat_entry(H, j, j));
                 for (j2 = j; j2 < n; j2++)
                 {
                     fmpz_submul(fmpz_mat_entry(H, k, j2), b, fmpz_mat_entry(H, j, j2));
                 }
                 continue;
            }

            fmpz_divexact(r1d, fmpz_mat_entry(H, j, j), d);
            fmpz_divexact(r2d, fmpz_mat_entry(H, k, j), d);
            for (j2 = j; j2 < n; j2++)
            {
                fmpz_mul(b, u, fmpz_mat_entry(H, j, j2));
                fmpz_addmul(b, v, fmpz_mat_entry(H, k, j2));
                fmpz_mul(fmpz_mat_entry(H, k, j2), r1d,
                         fmpz_mat_entry(H, k, j2));
                fmpz_submul(fmpz_mat_entry(H, k, j2), r2d,
                            fmpz_mat_entry(H, j, j2));
                fmpz_set(fmpz_mat_entry(H, j, j2), b);
            }
        }
        /* reduce above diagonal elements of each row i */
        for (i = n - 1; i >= 0; i--)
        {
            for (j = i + 1; j < n; j++)
            {
                fmpz_fdiv_q(q, fmpz_mat_entry(H, i, j),
                            fmpz_mat_entry(H, j, j));

                if (fmpz_is_zero(q))
                {
                    continue;
                }
                for (j2 = j; j2 < n; j2++)
                {
                    fmpz_submul(fmpz_mat_entry(H, i, j2), q,
                                fmpz_mat_entry(H, j, j2));
                }
            }
        }
    }

    fmpz_clear(b);
    fmpz_clear(q);
    fmpz_clear(r2d);
    fmpz_clear(r1d);
    fmpz_clear(d);
    fmpz_clear(v);
    fmpz_clear(u);
}
