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

    Copyright (C) 2014 Abhinav Baid

******************************************************************************/

#include "fmpq_vec.h"
#include "fmpq_mat.h"

void
fmpq_mat_lll(fmpz_mat_t B, const fmpz_mat_t A, const fmpq_lll_t fl)
{
    slong i, j, k, l, m, n;
    fmpz_t r, one;
    fmpq_t delta, nu, xi, half, rat;
    fmpq *gstar, *temp;
    fmpq_mat_t Bstar, mu;

    if (B->r != A->r || B->c != A->c)
    {
        flint_printf("Exception (fmpq_mat_lll). Incompatible dimensions.\n");
        abort();
    }

    if (B == A)
    {
        fmpz_mat_t t;
        fmpz_mat_init(t, B->r, B->c);
        fmpq_mat_lll(t, A, fl);
        fmpz_mat_swap(B, t);
        fmpz_mat_clear(t);
        return;
    }

    if (!A->r)
    {
        return;
    }

    fmpz_mat_set(B, A);
    m = B->r;
    n = B->c;
    fmpq_mat_init(Bstar, m, n);
    fmpq_mat_init(mu, m, m);
    gstar = _fmpq_vec_init(m);
    temp = _fmpq_vec_init(n);
    fmpz_init(r);
    fmpz_init_set_ui(one, 1);
    fmpq_init(delta);
    fmpq_init(nu);
    fmpq_init(xi);
    fmpq_init(half);
    fmpq_init(rat);
    fmpq_set_si(half, 1, 2);

    for (i = 0; i < m; i++)
    {
        for (k = 0; k < n; k++)
        {
            fmpq_set_fmpz_frac(fmpq_mat_entry(Bstar, i, k),
                               fmpz_mat_entry(B, i, k), one);
        }
        for (j = 0; j <= i - 1; j++)
        {
            _fmpq_vec_set_fmpz_vec(temp, B->rows[i], n);
            _fmpq_vec_dot(fmpq_mat_entry(mu, i, j), temp, Bstar->rows[j], n);
            fmpq_div(fmpq_mat_entry(mu, i, j), fmpq_mat_entry(mu, i, j),
                     gstar + j);
            for (k = 0; k < n; k++)
            {
                fmpq_submul(fmpq_mat_entry(Bstar, i, k),
                            fmpq_mat_entry(Bstar, j, k), fmpq_mat_entry(mu, i,
                                                                        j));
            }
        }
        _fmpq_vec_dot(gstar + i, Bstar->rows[i], Bstar->rows[i], n);
    }

    /* index k counts the current number of LLL-reduced rows */
    k = 1;
    while (k < m)
    {
        /* size reduce row k against row k - 1 */
        fmpq_abs(rat, fmpq_mat_entry(mu, k, k - 1));
        if (fmpq_cmp(rat, fl->eta) > 0)
        {
            /* determine reduction coefficient */
            fmpq_set(rat, fmpq_mat_entry(mu, k, k - 1));
            fmpq_sub(rat, rat, half);
            fmpz_cdiv_q(r, fmpq_numref(rat), fmpq_denref(rat));
            /* perform reduction */
            for (i = 0; i < n; i++)
                fmpz_submul(fmpz_mat_entry(B, k, i), r,
                            fmpz_mat_entry(B, k - 1, i));
            /* update mu */
            fmpq_set_fmpz_frac(rat, r, one);
            for (j = 0; j <= k - 2; j++)
                fmpq_submul(fmpq_mat_entry(mu, k, j), rat,
                            fmpq_mat_entry(mu, k - 1, j));
            fmpq_sub(fmpq_mat_entry(mu, k, k - 1),
                     fmpq_mat_entry(mu, k, k - 1), rat);
        }
        /* check exchange condition */
        fmpq_set(rat, fl->delta);
        fmpq_submul(rat, fmpq_mat_entry(mu, k, k - 1),
                    fmpq_mat_entry(mu, k, k - 1));
        fmpq_mul(rat, rat, gstar + k - 1);
        if (fmpq_cmp(gstar + k, rat) >= 0)
        {
            for (l = k - 2; l >= 0; l--)
            {
                /* size reduce row k against row l */
                fmpq_abs(rat, fmpq_mat_entry(mu, k, l));
                if (fmpq_cmp(rat, fl->eta) > 0)
                {
                    fmpq_set(rat, fmpq_mat_entry(mu, k, l));
                    fmpq_sub(rat, rat, half);
                    fmpz_cdiv_q(r, fmpq_numref(rat), fmpq_denref(rat));
                    for (i = 0; i < n; i++)
                        fmpz_submul(fmpz_mat_entry(B, k, i), r,
                                    fmpz_mat_entry(B, l, i));
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
            fmpq_set(delta, gstar + k - 1);
            fmpq_mul(delta, delta, nu);
            fmpq_mul(delta, delta, nu);
            fmpq_add(delta, delta, gstar + k);
            fmpq_set(fmpq_mat_entry(mu, k, k - 1), nu);
            fmpq_mul(fmpq_mat_entry(mu, k, k - 1),
                     fmpq_mat_entry(mu, k, k - 1), gstar + k - 1);
            fmpq_div(fmpq_mat_entry(mu, k, k - 1),
                     fmpq_mat_entry(mu, k, k - 1), delta);
            fmpq_mul(gstar + k, gstar + k, gstar + k - 1);
            fmpq_div(gstar + k, gstar + k, delta);
            fmpq_set(gstar + k - 1, delta);
            /* swap row k - 1 and row k */
            fmpz_mat_swap_rows(B, 0, k - 1, k);
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
    fmpq_clear(delta);
    fmpq_clear(nu);
    fmpq_clear(xi);
    fmpq_clear(half);
    fmpq_clear(rat);
    _fmpq_vec_clear(gstar, m);
    _fmpq_vec_clear(temp, n);
    fmpq_mat_clear(Bstar);
    fmpq_mat_clear(mu);
}
