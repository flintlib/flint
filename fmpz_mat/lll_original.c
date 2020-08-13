/*
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mat.h"
#include "fmpq_mat.h"

void
fmpz_mat_lll_original(fmpz_mat_t A, const fmpq_t delta, const fmpq_t eta)
{
    slong i, j, k, l, m, n;
    fmpz_t r, one;
    fmpq_t chi, nu, xi, half, rat;
    fmpq_mat_t R, mu;

    if (A->r == 0)
    {
        return;
    }

    m = A->r;
    n = A->c;
    fmpq_mat_init(R, m, m);
    fmpq_mat_init(mu, m, m);
    fmpz_init(r);
    fmpz_init_set_ui(one, 1);
    fmpq_init(chi);
    fmpq_init(nu);
    fmpq_init(xi);
    fmpq_init(half);
    fmpq_init(rat);
    fmpq_set_si(half, 1, 2);

    /* compute the rational GSO */
    for (i = 0; i < m; i++)
    {
        _fmpz_vec_dot(fmpq_mat_entry_num(mu, i, i), A->rows[i], A->rows[i], n);
        for (j = 0; j <= i - 1; j++)
        {
            _fmpz_vec_dot(fmpq_mat_entry_num(R, i, j), A->rows[i], A->rows[j],
                          n);
            for (k = 0; k <= j - 1; k++)
            {
                fmpq_submul(fmpq_mat_entry(R, i, j), fmpq_mat_entry(mu, j, k),
                            fmpq_mat_entry(R, i, k));
            }
            fmpq_div(fmpq_mat_entry(mu, i, j), fmpq_mat_entry(R, i, j),
                     fmpq_mat_entry(mu, j, j));
            fmpq_submul(fmpq_mat_entry(mu, i, i), fmpq_mat_entry(mu, i, j),
                        fmpq_mat_entry(R, i, j));
        }
    }

    /* index k counts the current number of LLL-reduced rows */
    k = 1;
    while (k < m)
    {
        /* size reduce row k against row k - 1 */
        fmpq_abs(rat, fmpq_mat_entry(mu, k, k - 1));
        if (fmpq_cmp(rat, eta) > 0)
        {
            /* determine reduction coefficient */
            fmpq_sub(rat, fmpq_mat_entry(mu, k, k - 1), half);
            fmpz_cdiv_q(r, fmpq_numref(rat), fmpq_denref(rat));
            /* perform reduction */
            for (i = 0; i < n; i++)
                fmpz_submul(fmpz_mat_entry(A, k, i), r,
                            fmpz_mat_entry(A, k - 1, i));
            /* update mu */
            fmpq_set_fmpz_frac(rat, r, one);
            for (j = 0; j <= k - 2; j++)
                fmpq_submul(fmpq_mat_entry(mu, k, j), rat,
                            fmpq_mat_entry(mu, k - 1, j));
            fmpq_sub(fmpq_mat_entry(mu, k, k - 1),
                     fmpq_mat_entry(mu, k, k - 1), rat);
        }
        /* check exchange condition */
        fmpq_set(rat, delta);
        fmpq_submul(rat, fmpq_mat_entry(mu, k, k - 1),
                    fmpq_mat_entry(mu, k, k - 1));
        fmpq_mul(rat, rat, fmpq_mat_entry(mu, k - 1, k - 1));
        if (fmpq_cmp(fmpq_mat_entry(mu, k, k), rat) >= 0)
        {
            for (l = k - 2; l >= 0; l--)
            {
                /* size reduce row k against row l */
                fmpq_abs(rat, fmpq_mat_entry(mu, k, l));
                if (fmpq_cmp(rat, eta) > 0)
                {
                    fmpq_sub(rat, fmpq_mat_entry(mu, k, l), half);
                    fmpz_cdiv_q(r, fmpq_numref(rat), fmpq_denref(rat));
                    for (i = 0; i < n; i++)
                        fmpz_submul(fmpz_mat_entry(A, k, i), r,
                                    fmpz_mat_entry(A, l, i));
                    fmpq_set_fmpz_frac(rat, r, one);
                    for (j = 0; j <= l - 1; j++)
                        fmpq_submul(fmpq_mat_entry(mu, k, j), rat,
                                    fmpq_mat_entry(mu, l, j));
                    fmpq_sub(fmpq_mat_entry(mu, k, l),
                             fmpq_mat_entry(mu, k, l), rat);
                }
            }
            /* increment LLL-reduced index */
            k++;
        }
        else
        {
            fmpq_set(nu, fmpq_mat_entry(mu, k, k - 1));
            fmpq_mul(chi, fmpq_mat_entry(mu, k - 1, k - 1), nu);
            fmpq_mul(chi, chi, nu);
            fmpq_add(chi, chi, fmpq_mat_entry(mu, k, k));
            fmpq_mul(fmpq_mat_entry(mu, k, k - 1),
                     fmpq_mat_entry(mu, k, k - 1), fmpq_mat_entry(mu, k - 1,
                                                                  k - 1));
            fmpq_div(fmpq_mat_entry(mu, k, k - 1),
                     fmpq_mat_entry(mu, k, k - 1), chi);
            fmpq_mul(fmpq_mat_entry(mu, k, k), fmpq_mat_entry(mu, k, k),
                     fmpq_mat_entry(mu, k - 1, k - 1));
            fmpq_div(fmpq_mat_entry(mu, k, k), fmpq_mat_entry(mu, k, k), chi);
            fmpq_set(fmpq_mat_entry(mu, k - 1, k - 1), chi);
            /* swap row k - 1 and row k */
            fmpz_mat_swap_rows(A, NULL, k - 1, k);
            /* update mu */
            for (j = 0; j <= k - 2; j++)
            {
                fmpq_swap(fmpq_mat_entry(mu, k - 1, j),
                          fmpq_mat_entry(mu, k, j));
            }
            for (i = k + 1; i < m; i++)
            {
                fmpq_set(xi, fmpq_mat_entry(mu, i, k));
                fmpq_set(fmpq_mat_entry(mu, i, k),
                         fmpq_mat_entry(mu, i, k - 1));
                fmpq_submul(fmpq_mat_entry(mu, i, k), nu, xi);
                fmpq_mul(fmpq_mat_entry(mu, i, k - 1),
                         fmpq_mat_entry(mu, k, k - 1), fmpq_mat_entry(mu, i,
                                                                      k));
                fmpq_add(fmpq_mat_entry(mu, i, k - 1),
                         fmpq_mat_entry(mu, i, k - 1), xi);
            }
            /* decrement LLL-reduced index */
            if (k > 1)
                k--;
        }
    }

    fmpz_clear(r);
    fmpz_clear(one);
    fmpq_clear(chi);
    fmpq_clear(nu);
    fmpq_clear(xi);
    fmpq_clear(half);
    fmpq_clear(rat);
    fmpq_mat_clear(R);
    fmpq_mat_clear(mu);
}
