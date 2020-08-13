/*
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mat.h"
#include "fmpq_vec.h"
#include "fmpq_mat.h"

int
fmpz_mat_is_reduced(const fmpz_mat_t A, double delta, double eta)
{
    slong i, j, k, d = A->r, n = A->c;
    fmpq_mat_t Aq, Bq, mu;
    mpq_t deltax, etax;
    fmpq_t deltaq, etaq, tmp;

    if (d == 0 || d == 1)
        return 1;

    fmpq_mat_init(Aq, d, n);
    fmpq_mat_init(Bq, d, n);
    fmpq_mat_init(mu, d, d);

    mpq_init(deltax);
    mpq_init(etax);

    fmpq_init(deltaq);
    fmpq_init(etaq);
    fmpq_init(tmp);

    fmpq_mat_set_fmpz_mat(Aq, A);

    mpq_set_d(deltax, delta);
    mpq_set_d(etax, eta);
    fmpq_set_mpq(deltaq, deltax);
    fmpq_set_mpq(etaq, etax);
    mpq_clears(deltax, etax, NULL);

    for (j = 0; j < n; j++)
    {
        fmpq_set(fmpq_mat_entry(Bq, 0, j), fmpq_mat_entry(Aq, 0, j));
    }
    /* diagonal of mu stores the squared GS norms */
    _fmpq_vec_dot(fmpq_mat_entry(mu, 0, 0), Bq->rows[0], Bq->rows[0], n);

    for (i = 1; i < d; i++)
    {
        for (j = 0; j < n; j++)
        {
            fmpq_set(fmpq_mat_entry(Bq, i, j), fmpq_mat_entry(Aq, i, j));
        }

        for (j = 0; j < i; j++)
        {
            _fmpq_vec_dot(tmp, Aq->rows[i], Bq->rows[j], n);

            fmpq_div(fmpq_mat_entry(mu, i, j), tmp, fmpq_mat_entry(mu, j, j));

            for (k = 0; k < n; k++)
            {
                fmpq_submul(fmpq_mat_entry(Bq, i, k),
                            fmpq_mat_entry(mu, i, j), fmpq_mat_entry(Bq, j,
                                                                     k));
            }
            fmpq_abs(tmp, fmpq_mat_entry(mu, i, j));
            if (fmpq_cmp(tmp, etaq) > 0)    /* check size reduction */
            {
                fmpq_mat_clear(Aq);
                fmpq_mat_clear(Bq);
                fmpq_mat_clear(mu);
                fmpq_clear(deltaq);
                fmpq_clear(etaq);
                fmpq_clear(tmp);
                return 0;
            }
        }
        fmpq_set(tmp, deltaq);
        fmpq_submul(tmp, fmpq_mat_entry(mu, i, i - 1),
                    fmpq_mat_entry(mu, i, i - 1));
        fmpq_mul(tmp, tmp, fmpq_mat_entry(mu, i - 1, i - 1));
        _fmpq_vec_dot(fmpq_mat_entry(mu, i, i), Bq->rows[i], Bq->rows[i], n);
        if (fmpq_cmp(tmp, fmpq_mat_entry(mu, i, i)) > 0)    /* check Lovasz condition */
        {
            fmpq_mat_clear(Aq);
            fmpq_mat_clear(Bq);
            fmpq_mat_clear(mu);
            fmpq_clear(deltaq);
            fmpq_clear(etaq);
            fmpq_clear(tmp);
            return 0;
        }
    }
    fmpq_mat_clear(Aq);
    fmpq_mat_clear(Bq);
    fmpq_mat_clear(mu);
    fmpq_clear(deltaq);
    fmpq_clear(etaq);
    fmpq_clear(tmp);
    return 1;
}
