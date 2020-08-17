/*
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mat.h"

void
fmpz_mat_lll_storjohann(fmpz_mat_t A, const fmpq_t delta, const fmpq_t eta)
{
    slong n, np, i, j, k;
    double e;
    fmpz_t M, lhs, rhs;
    fmpz_mat_t T;
    fmpq_t max, gsn, half;

    if (A->r == 0)
    {
        return;
    }

    n = A->r;
    np = A->c;

    fmpz_init(M);
    fmpz_init(lhs);
    fmpz_init(rhs);

    fmpz_mat_init(T, n, n);
    fmpz_mat_gram(T, A);
    for (i = 0; i < n - 1; i++)
    {
        for (j = i + 1; j < n; j++)
        {
            for (k = i + 1; k < n; k++)
            {
                fmpz_mul(fmpz_mat_entry(T, j, k), fmpz_mat_entry(T, j, k),
                         fmpz_mat_entry(T, i, i));
                fmpz_submul(fmpz_mat_entry(T, j, k), fmpz_mat_entry(T, j, i),
                            fmpz_mat_entry(T, i, k));

                if (i > 0)
                {
                    fmpz_divexact(fmpz_mat_entry(T, j, k),
                                  fmpz_mat_entry(T, j, k), fmpz_mat_entry(T,
                                                                          i -
                                                                          1,
                                                                          i -
                                                                          1));
                }
            }

            fmpz_zero(fmpz_mat_entry(T, j, i));
        }
    }

    fmpq_init(max);
    fmpq_init(gsn);
    fmpq_init(half);
    fmpq_set_si(half, 1, 2);

    fmpz_set(fmpq_numref(max), fmpz_mat_entry(T, 0, 0));
    fmpz_one(fmpq_denref(max));
    for (i = 1; i < n; i++)
    {
        fmpq_set_fmpz_frac(gsn, fmpz_mat_entry(T, i, i),
                           fmpz_mat_entry(T, i - 1, i - 1));
        if (fmpq_cmp(gsn, max) > 0)
        {
            fmpq_set(max, gsn);
        }
    }
    fmpz_set_si(M, n);
    fmpq_mul_fmpz(max, max, M);
    e = 2 *
        ceil(sqrt
             (n * fmpz_get_d(fmpq_numref(max)) /
              fmpz_get_d(fmpq_denref(max)))) + 1;
    fmpz_set_d(M, e);
    k = 1;

    while (k < n)
    {
        fmpq_set_fmpz_frac(max, fmpz_mat_entry(T, k - 1, k),
                           fmpz_mat_entry(T, k - 1, k - 1));
        fmpq_abs(gsn, max);
        if (fmpq_cmp(gsn, eta) > 0)
        {
            fmpq_sub(max, max, half);
            fmpz_cdiv_q(lhs, fmpq_numref(max), fmpq_denref(max));
            _fmpz_vec_scalar_submul_fmpz(A->rows[k], A->rows[k - 1], np, lhs);
            for (i = 0; i < n; i++)
            {
                fmpz_submul(fmpz_mat_entry(T, i, k), lhs,
                            fmpz_mat_entry(T, i, k - 1));
                if (i <= k - 1)
                {
                    fmpz_mul(rhs, fmpz_mat_entry(T, i, i), M);
                    if (i > 0)
                    {
                        fmpz_mul(rhs, rhs, fmpz_mat_entry(T, i - 1, i - 1));
                    }
                    fmpz_smod(fmpz_mat_entry(T, i, k), fmpz_mat_entry(T, i, k),
                              rhs);
                }
            }
            for (j = 0; j < np; j++)
            {
                fmpz_smod(fmpz_mat_entry(A, k, j), fmpz_mat_entry(A, k, j), M);
            }
        }
        fmpq_set_fmpz_frac(max, fmpz_mat_entry(T, k, k),
                           fmpz_mat_entry(T, k - 1, k - 1));
        fmpq_div_fmpz(max, max, fmpz_mat_entry(T, k - 1, k - 1));
        if (k > 1)
        {
            fmpq_mul_fmpz(max, max, fmpz_mat_entry(T, k - 2, k - 2));
        }
        fmpq_set(gsn, delta);
        _fmpq_submul(fmpq_numref(gsn), fmpq_denref(gsn),
                     fmpz_mat_entry(T, k - 1, k), fmpz_mat_entry(T, k - 1,
                                                                 k - 1),
                     fmpz_mat_entry(T, k - 1, k), fmpz_mat_entry(T, k - 1,
                                                                 k - 1));
        if (fmpq_cmp(max, gsn) < 0)
        {
            fmpz_mat_swap_rows(A, NULL, k - 1, k);
            if (k > 1)
            {
                _fmpz_vec_scalar_mul_fmpz(T->rows[k], T->rows[k], n,
                                          fmpz_mat_entry(T, k - 2, k - 2));
            }
            _fmpz_vec_scalar_addmul_fmpz(T->rows[k], T->rows[k - 1], n,
                                         fmpz_mat_entry(T, k - 1, k));
            _fmpz_vec_scalar_divexact_fmpz(T->rows[k], T->rows[k], n,
                                           fmpz_mat_entry(T, k - 1, k - 1));
            fmpz_mat_swap_rows(T, NULL, k - 1, k);
            for (i = 0; i < n; i++)
            {
                fmpz_swap(fmpz_mat_entry(T, i, k - 1),
                          fmpz_mat_entry(T, i, k));
            }
            _fmpz_vec_scalar_mul_fmpz(T->rows[k], T->rows[k], n,
                                      fmpz_mat_entry(T, k - 1, k - 1));
            _fmpz_vec_scalar_submul_fmpz(T->rows[k], T->rows[k - 1], n,
                                         fmpz_mat_entry(T, k - 1, k));
            if (k > 1)
            {
                _fmpz_vec_scalar_divexact_fmpz(T->rows[k], T->rows[k], n,
                                               fmpz_mat_entry(T, k - 2,
                                                              k - 2));
            }
            for (i = 0; i <= k - 2; i++)
            {
                fmpz_mul(rhs, fmpz_mat_entry(T, i, i), M);
                if (i > 0)
                {
                    fmpz_mul(rhs, rhs, fmpz_mat_entry(T, i - 1, i - 1));
                }
                fmpz_smod(fmpz_mat_entry(T, i, k - 1),
                          fmpz_mat_entry(T, i, k - 1), rhs);
                fmpz_smod(fmpz_mat_entry(T, i, k), fmpz_mat_entry(T, i, k),
                          rhs);
            }
            fmpz_mul(rhs, fmpz_mat_entry(T, k - 1, k - 1), M);
            fmpz_mul(lhs, rhs, fmpz_mat_entry(T, k, k));
            if (k > 1)
            {
                fmpz_mul(rhs, rhs, fmpz_mat_entry(T, k - 2, k - 2));
            }
            fmpz_smod(fmpz_mat_entry(T, k - 1, k), fmpz_mat_entry(T, k - 1, k),
                      rhs);
            for (j = k + 1; j < n; j++)
            {
                fmpz_smod(fmpz_mat_entry(T, k - 1, j),
                          fmpz_mat_entry(T, k - 1, j), rhs);
                fmpz_smod(fmpz_mat_entry(T, k, j), fmpz_mat_entry(T, k, j),
                          lhs);
            }
            if (k > 1)
                k--;
        }
        else
        {
            k++;
        }
    }

    for (k = 1; k < n; k++)
    {
        for (j = k - 1; j >= 0; j--)
        {
            fmpq_set_fmpz_frac(max, fmpz_mat_entry(T, j, k),
                               fmpz_mat_entry(T, j, j));
            fmpq_abs(gsn, max);
            if (fmpq_cmp(gsn, eta) > 0)
            {
                fmpq_sub(max, max, half);
                fmpz_cdiv_q(lhs, fmpq_numref(max), fmpq_denref(max));
                _fmpz_vec_scalar_submul_fmpz(A->rows[k], A->rows[j], np, lhs);
                for (i = 0; i < n; i++)
                {
                    fmpz_submul(fmpz_mat_entry(T, i, k), lhs,
                                fmpz_mat_entry(T, i, j));
                    if (i <= k - 1)
                    {
                        fmpz_mul(rhs, fmpz_mat_entry(T, i, i), M);
                        if (i > 0)
                        {
                            fmpz_mul(rhs, rhs,
                                     fmpz_mat_entry(T, i - 1, i - 1));
                        }
                        fmpz_smod(fmpz_mat_entry(T, i, k),
                                  fmpz_mat_entry(T, i, k), rhs);
                    }
                }
                for (i = 0; i < np; i++)
                {
                    fmpz_smod(fmpz_mat_entry(A, k, i), fmpz_mat_entry(A, k, i),
                              M);
                }
            }
        }
    }

    fmpz_clear(M);
    fmpz_clear(lhs);
    fmpz_clear(rhs);

    fmpz_mat_clear(T);

    fmpq_clear(max);
    fmpq_clear(gsn);
    fmpq_clear(half);
}
