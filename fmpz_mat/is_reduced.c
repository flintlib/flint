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

#include "fmpz_mat.h"
#include "fmpq_vec.h"
#include "fmpq_mat.h"

int
fmpz_mat_is_reduced(const fmpz_mat_t A, double delta, double eta)
{
    slong i, j, k, d = A->r, n = A->c;
    fmpq_mat_t Aq, Bq, mu;
    mpq_t deltax, etax;
    fmpq_t deltaq, etaq, tmp, rtmp;

    if (d == 1)
        return 1;

    fmpq_mat_init(Aq, d, n);
    fmpq_mat_init(Bq, d, n);
    fmpq_mat_init(mu, d, d);

    mpq_init(deltax);
    mpq_init(etax);

    fmpq_init(deltaq);
    fmpq_init(etaq);
    fmpq_init(tmp);
    fmpq_init(rtmp);

    fmpq_mat_set_fmpz_mat(Aq, A);

    mpq_set_d(deltax, delta);
    mpq_set_d(etax, eta);
    fmpq_set_mpq(deltaq, deltax);
    fmpq_set_mpq(etaq, etax);
    mpq_clears(deltax, etax, '\0');

    for (i = 0; i < d; i++)
    {
        for (j = 0; j < n; j++)
        {
            fmpq_set(fmpq_mat_entry(Bq, i, j), fmpq_mat_entry(Aq, i, j));
        }

        for (j = 0; j < i; j++)
        {
            _fmpq_vec_dot(tmp, Aq->rows[i], Bq->rows[j], n);
            _fmpq_vec_dot(rtmp, Bq->rows[j], Bq->rows[j], n);

            fmpq_div(fmpq_mat_entry(mu, i, j), tmp, rtmp);

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
                fmpq_clear(rtmp);
                return 0;
            }
        }
        if (i > 0)
        {
            fmpq_set(tmp, deltaq);
            fmpq_submul(tmp, fmpq_mat_entry(mu, i, i - 1),
                        fmpq_mat_entry(mu, i, i - 1));
            fmpq_mul(tmp, tmp, rtmp);
            _fmpq_vec_dot(rtmp, Bq->rows[i], Bq->rows[i], n);
            if (fmpq_cmp(tmp, rtmp) > 0)    /* check Lovasz condition */
            {
                fmpq_mat_clear(Aq);
                fmpq_mat_clear(Bq);
                fmpq_mat_clear(mu);
                fmpq_clear(deltaq);
                fmpq_clear(etaq);
                fmpq_clear(tmp);
                fmpq_clear(rtmp);
                return 0;
            }
        }
    }
    fmpq_mat_clear(Aq);
    fmpq_mat_clear(Bq);
    fmpq_mat_clear(mu);
    fmpq_clear(deltaq);
    fmpq_clear(etaq);
    fmpq_clear(tmp);
    fmpq_clear(rtmp);
    return 1;
}
