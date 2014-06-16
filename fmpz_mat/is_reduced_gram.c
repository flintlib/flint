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
    fmpq_mat_t r, mu;
    fmpq *s;
    mpq_t deltax, etax;
    fmpq_t deltaq, etaq, tmp;
    fmpz_t one;

    if (d == 1)
        return 1;

    fmpq_mat_init(r, d, d);
    fmpq_mat_init(mu, d, d);

    s = _fmpq_vec_init(d);

    mpq_init(deltax);
    mpq_init(etax);

    fmpq_init(deltaq);
    fmpq_init(etaq);
    fmpq_init(tmp);

    fmpz_init_set_ui(one, 1);

    mpq_set_d(deltax, delta);
    mpq_set_d(etax, eta);
    fmpq_set_mpq(deltaq, deltax);
    fmpq_set_mpq(etaq, etax);
    mpq_clears(deltax, etax, '\0');

    for (i = 0; i < d; i++)
    {
        for (j = 0; j <= i - 1; j++)
        {
            fmpq_set_fmpz_frac(fmpq_mat_entry(r, i, j),
                               fmpz_mat_entry(A, i, j), one);
            for (k = 0; k <= j - 1; k++)
            {
                fmpq_submul(fmpq_mat_entry(r, i, j), fmpq_mat_entry(mu, j, k),
                            fmpq_mat_entry(r, i, k));
            }
            fmpq_div(fmpq_mat_entry(mu, i, j), fmpq_mat_entry(r, i, j),
                     fmpq_mat_entry(r, j, j));
            fmpq_abs(tmp, fmpq_mat_entry(mu, i, j));
            if (fmpq_cmp(tmp, etaq) > 0)    /* check size reduction */
            {
                fmpq_mat_clear(r);
                fmpq_mat_clear(mu);
                fmpq_clear(deltaq);
                fmpq_clear(etaq);
                fmpq_clear(tmp);
                _fmpq_vec_clear(s, d);
                fmpz_clear(one);
                return 0;
            }
        }
        fmpq_set_fmpz_frac(s, fmpz_mat_entry(A, i, i), one);
        for (j = 1; j <= i; j++)
        {
            fmpq_set(s + j, s + j - 1);
            fmpq_submul(s + j, fmpq_mat_entry(mu, i, j - 1),
                        fmpq_mat_entry(r, i, j - 1));
        }
        fmpq_set(fmpq_mat_entry(r, i, i), s + i);
        if (i > 0)
        {
            fmpq_mul(tmp, deltaq, fmpq_mat_entry(r, i - 1, i - 1));
            if (fmpq_cmp(tmp, s + i - 1) > 0)   /* check Lovasz condition */
            {
                fmpq_mat_clear(r);
                fmpq_mat_clear(mu);
                fmpq_clear(deltaq);
                fmpq_clear(etaq);
                fmpq_clear(tmp);
                _fmpq_vec_clear(s, d);
                fmpz_clear(one);
                return 0;
            }
        }
    }
    fmpq_mat_clear(r);
    fmpq_mat_clear(mu);
    fmpq_clear(deltaq);
    fmpq_clear(etaq);
    fmpq_clear(tmp);
    _fmpq_vec_clear(s, d);
    fmpz_clear(one);
    return 1;
}
