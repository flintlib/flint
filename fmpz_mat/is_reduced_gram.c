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
#include "fmpq_mat.h"

int
fmpz_mat_is_reduced_gram(const fmpz_mat_t A, double delta, double eta)
{
    slong i, j, k, d = A->r;
    fmpq_mat_t mu;
    mpq_t deltax, etax;
    fmpq_t deltaq, etaq, tmp, rtmp;
    fmpz_t one;

    if (d == 1)
        return 1;

    fmpq_mat_init(mu, d, d);

    mpq_init(deltax);
    mpq_init(etax);

    fmpq_init(deltaq);
    fmpq_init(etaq);
    fmpq_init(tmp);
    fmpq_init(rtmp);

    fmpz_init_set_ui(one, 1);

    mpq_set_d(deltax, delta);
    mpq_set_d(etax, eta);
    fmpq_set_mpq(deltaq, deltax);
    fmpq_set_mpq(etaq, etax);
    mpq_clears(deltax, etax, '\0');

    for (i = 0; i < d; i++)
    {
        fmpq_set_fmpz_frac(fmpq_mat_entry(mu, i, i), fmpz_mat_entry(A, i, i),
                           one);
        for (j = 0; j <= i - 1; j++)
        {
            fmpq_set_fmpz_frac(fmpq_mat_entry(mu, i, j),
                               fmpz_mat_entry(A, i, j), one);
            for (k = 0; k <= j - 1; k++)
            {
                fmpq_mul(tmp, fmpq_mat_entry(mu, i, k),
                         fmpq_mat_entry(mu, k, k));
                fmpq_submul(fmpq_mat_entry(mu, i, j), fmpq_mat_entry(mu, j, k),
                            tmp);
            }
            fmpq_set(tmp, fmpq_mat_entry(mu, i, j));
            fmpq_div(fmpq_mat_entry(mu, i, j), fmpq_mat_entry(mu, i, j),
                     fmpq_mat_entry(mu, j, j));
            fmpq_submul(fmpq_mat_entry(mu, i, i), fmpq_mat_entry(mu, i, j),
                        tmp);
            fmpq_abs(rtmp, fmpq_mat_entry(mu, i, j));
            if (fmpq_cmp(rtmp, etaq) > 0)   /* check size reduction */
            {
                fmpq_mat_clear(mu);
                fmpq_clear(deltaq);
                fmpq_clear(etaq);
                fmpq_clear(tmp);
                fmpq_clear(rtmp);
                fmpz_clear(one);
                return 0;
            }
        }
        if (i > 0)
        {
            fmpq_set(rtmp, fmpq_mat_entry(mu, i, i));
            fmpq_addmul(rtmp, fmpq_mat_entry(mu, i, i - 1), tmp);
            fmpq_mul(tmp, deltaq, fmpq_mat_entry(mu, i - 1, i - 1));
            if (fmpq_cmp(tmp, rtmp) > 0)    /* check Lovasz condition */
            {
                fmpq_mat_clear(mu);
                fmpq_clear(deltaq);
                fmpq_clear(etaq);
                fmpq_clear(tmp);
                fmpq_clear(rtmp);
                fmpz_clear(one);
                return 0;
            }
        }
    }
    fmpq_mat_clear(mu);
    fmpq_clear(deltaq);
    fmpq_clear(etaq);
    fmpq_clear(tmp);
    fmpq_clear(rtmp);
    fmpz_clear(one);
    return 1;
}
