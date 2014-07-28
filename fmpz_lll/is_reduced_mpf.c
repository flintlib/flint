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

#include "fmpz_lll.h"

int
fmpz_lll_is_reduced_mpf(const fmpz_mat_t B, const fmpz_lll_t fl,
                        mp_bitcnt_t prec)
{
    if (fl->rt == Z_BASIS)
    {
        slong i, j, k, d = B->r, n = B->c;
        mpf_mat_t appB, Q, mu;
        mpf_t tmp, rtmp;

        if (d == 0 || d == 1)
            return 1;

        mpf_mat_init(appB, d, n, prec);
        mpf_mat_init(Q, d, n, prec);
        mpf_mat_init(mu, d, d, prec);

        mpf_init2(tmp, prec);
        mpf_init2(rtmp, prec);

        fmpz_mat_get_mpf_mat(appB, B);

        for (j = 0; j < n; j++)
        {
            mpf_set(mpf_mat_entry(Q, 0, j), mpf_mat_entry(appB, 0, j));
        }
        /* diagonal of mu stores the squared GS norms */
        _mpf_vec_norm2(mpf_mat_entry(mu, 0, 0), Q->rows[0], n, prec);

        for (i = 1; i < d; i++)
        {
            for (j = 0; j < n; j++)
            {
                mpf_set(mpf_mat_entry(Q, i, j), mpf_mat_entry(appB, i, j));
            }

            for (j = 0; j < i; j++)
            {
                _mpf_vec_dot2(tmp, appB->rows[i], Q->rows[j], n, prec);
                mpf_div(mpf_mat_entry(mu, i, j), tmp, mpf_mat_entry(mu, j, j));

                for (k = 0; k < n; k++)
                {
                    mpf_mul(tmp, mpf_mat_entry(mu, i, j),
                            mpf_mat_entry(Q, j, k));
                    mpf_sub(mpf_mat_entry(Q, i, k), mpf_mat_entry(Q, i, k),
                            tmp);
                }

                mpf_abs(tmp, mpf_mat_entry(mu, i, j));
                if (mpf_cmp_d(tmp, fl->eta) > 0)    /* check size reduction */
                {
                    mpf_mat_clear(appB);
                    mpf_mat_clear(Q);
                    mpf_mat_clear(mu);
                    mpf_clears(tmp, rtmp, '\0');
                    return 0;
                }
            }

            _mpf_vec_norm2(mpf_mat_entry(mu, i, i), Q->rows[i], n, prec);
            mpf_set_d(rtmp, fl->delta);
            mpf_pow_ui(tmp, tmp, 2);
            mpf_sub(rtmp, rtmp, tmp);
            mpf_mul(tmp, rtmp, mpf_mat_entry(mu, i - 1, i - 1));
            if (mpf_cmp(tmp, mpf_mat_entry(mu, i, i)) > 0)  /* check Lovasz condition */
            {
                mpf_mat_clear(appB);
                mpf_mat_clear(Q);
                mpf_mat_clear(mu);
                mpf_clears(tmp, rtmp, '\0');
                return 0;
            }
        }
        mpf_mat_clear(appB);
        mpf_mat_clear(Q);
        mpf_mat_clear(mu);
        mpf_clears(tmp, rtmp, '\0');
    }
    else
    {
        slong i, j, k, d = B->r;
        mpf_mat_t r, mu;
        mpf *s;
        mpf_t tmp;

        if (d == 0 || d == 1)
            return 1;

        mpf_mat_init(r, d, d, prec);
        mpf_mat_init(mu, d, d, prec);

        s = _mpf_vec_init(d, prec);

        mpf_init2(tmp, prec);

        fmpz_get_mpf(mpf_mat_entry(r, 0, 0), fmpz_mat_entry(B, 0, 0));

        for (i = 1; i < d; i++)
        {
            fmpz_get_mpf(s, fmpz_mat_entry(B, i, i));
            for (j = 0; j <= i - 1; j++)
            {
                fmpz_get_mpf(mpf_mat_entry(r, i, j), fmpz_mat_entry(B, i, j));
                for (k = 0; k <= j - 1; k++)
                {
                    mpf_mul(tmp, mpf_mat_entry(mu, j, k),
                            mpf_mat_entry(r, i, k));
                    mpf_sub(mpf_mat_entry(r, i, j), mpf_mat_entry(r, i, j),
                            tmp);
                }
                mpf_div(mpf_mat_entry(mu, i, j), mpf_mat_entry(r, i, j),
                        mpf_mat_entry(r, j, j));
                mpf_abs(tmp, mpf_mat_entry(mu, i, j));
                if (mpf_cmp_d(tmp, fl->eta) > 0)    /* check size reduction */
                {
                    mpf_mat_clear(r);
                    mpf_mat_clear(mu);
                    _mpf_vec_clear(s, d);
                    mpf_clear(tmp);
                    return 0;
                }
                mpf_mul(tmp, mpf_mat_entry(mu, i, j), mpf_mat_entry(r, i, j));
                mpf_sub(s + j + 1, s + j, tmp);
            }
            mpf_set(mpf_mat_entry(r, i, i), s + i);
            mpf_set_d(tmp, fl->delta);
            mpf_mul(tmp, tmp, mpf_mat_entry(r, i - 1, i - 1));
            if (mpf_cmp(tmp, s + i - 1) > 0)    /* check Lovasz condition */
            {
                mpf_mat_clear(r);
                mpf_mat_clear(mu);
                _mpf_vec_clear(s, d);
                mpf_clear(tmp);
                return 0;
            }
        }
        mpf_mat_clear(r);
        mpf_mat_clear(mu);
        _mpf_vec_clear(s, d);
        mpf_clear(tmp);
    }
    return 1;
}
