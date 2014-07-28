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
fmpz_lll_is_reduced_d_with_removal(const fmpz_mat_t B, const fmpz_lll_t fl,
                                   const fmpz_t gs_B, int newd)
{
    if (fl->rt == Z_BASIS)
    {
        slong i, j, k, d = B->r, n = B->c;
        d_mat_t appB, Q, mu;
        double d_gs_B;

        if (d == 0 || d == 1)
            return 1;

        d_mat_init(appB, d, n);
        d_mat_init(Q, d, n);
        d_mat_init(mu, d, d);

        if (fmpz_mat_get_d_mat(appB, B) == -1)
        {
            d_mat_clear(appB);
            d_mat_clear(Q);
            d_mat_clear(mu);
            return 0;
        }

        d_gs_B = fmpz_get_d(gs_B);

        for (j = 0; j < n; j++)
        {
            d_mat_entry(Q, 0, j) = d_mat_entry(appB, 0, j);
        }
        /* diagonal of mu stores the squared GS norms */
        d_mat_entry(mu, 0, 0) = _d_vec_norm(Q->rows[0], n);
        if (newd == 0 && d_mat_entry(mu, 0, 0) < d_gs_B)
        {
            d_mat_clear(appB);
            d_mat_clear(Q);
            d_mat_clear(mu);
            return 0;
        }

        for (i = 1; i < d; i++)
        {
            for (j = 0; j < n; j++)
            {
                d_mat_entry(Q, i, j) = d_mat_entry(appB, i, j);
            }

            for (j = 0; j < i; j++)
            {
                d_mat_entry(mu, i, j) =
                    _d_vec_dot_heuristic(appB->rows[i], Q->rows[j], n,
                                         NULL) / d_mat_entry(mu, j, j);

                for (k = 0; k < n; k++)
                {
                    d_mat_entry(Q, i, k) -=
                        d_mat_entry(mu, i, j) * d_mat_entry(Q, j, k);
                }

                if (i < newd)
                {
                    if (fabs(d_mat_entry(mu, i, j)) > fl->eta)  /* check size reduction */
                    {
                        d_mat_clear(appB);
                        d_mat_clear(Q);
                        d_mat_clear(mu);
                        return 0;
                    }
                }
            }

            d_mat_entry(mu, i, i) = _d_vec_norm(Q->rows[i], n);
            if (i >= newd && d_mat_entry(mu, i, i) < d_gs_B)    /* check removals */
            {
                d_mat_clear(appB);
                d_mat_clear(Q);
                d_mat_clear(mu);
                return 0;
            }
            if (i < newd)
            {
                if ((fl->delta - pow(d_mat_entry(mu, i, i - 1), 2)) * d_mat_entry(mu, i - 1, i - 1) > d_mat_entry(mu, i, i))    /* check Lovasz condition */
                {
                    d_mat_clear(appB);
                    d_mat_clear(Q);
                    d_mat_clear(mu);
                    return 0;
                }
            }
        }
        d_mat_clear(appB);
        d_mat_clear(Q);
        d_mat_clear(mu);
    }
    else
    {
        slong i, j, k, d = B->r;
        d_mat_t r, mu;
        double *s;
        double d_gs_B;
        fmpz_t dmax;

        if (d == 0 || d == 1)
            return 1;

        d_mat_init(r, d, d);
        d_mat_init(mu, d, d);

        s = _d_vec_init(d);

        fmpz_init(dmax);
        fmpz_set_d(dmax, DBL_MAX);

        if (fmpz_cmpabs(fmpz_mat_entry(B, 0, 0), dmax) > 0)
        {
            d_mat_clear(r);
            d_mat_clear(mu);
            _d_vec_clear(s);
            fmpz_clear(dmax);
            return 0;
        }
        d_mat_entry(r, 0, 0) = fmpz_get_d(fmpz_mat_entry(B, 0, 0));
        if (newd == 0 && fmpz_cmp(fmpz_mat_entry(B, 0, 0), gs_B) < 0)
        {
            d_mat_clear(r);
            d_mat_clear(mu);
            _d_vec_clear(s);
            fmpz_clear(dmax);
            return 0;
        }

        d_gs_B = fmpz_get_d(gs_B);
        for (i = 1; i < d; i++)
        {
            if (fmpz_cmpabs(fmpz_mat_entry(B, i, i), dmax) > 0)
            {
                d_mat_clear(r);
                d_mat_clear(mu);
                _d_vec_clear(s);
                fmpz_clear(dmax);
                return 0;
            }
            s[0] = fmpz_get_d(fmpz_mat_entry(B, i, i));
            for (j = 0; j <= i - 1; j++)
            {
                if (fmpz_cmpabs(fmpz_mat_entry(B, i, j), dmax) > 0)
                {
                    d_mat_clear(r);
                    d_mat_clear(mu);
                    _d_vec_clear(s);
                    fmpz_clear(dmax);
                    return 0;
                }
                d_mat_entry(r, i, j) = fmpz_get_d(fmpz_mat_entry(B, i, j));
                for (k = 0; k <= j - 1; k++)
                {
                    d_mat_entry(r, i, j) -=
                        d_mat_entry(mu, j, k) * d_mat_entry(r, i, k);
                }
                d_mat_entry(mu, i, j) =
                    d_mat_entry(r, i, j) / d_mat_entry(r, j, j);
                if (i < newd)
                {
                    if (fabs(d_mat_entry(mu, i, j)) > fl->eta)  /* check size reduction */
                    {
                        d_mat_clear(r);
                        d_mat_clear(mu);
                        _d_vec_clear(s);
                        fmpz_clear(dmax);
                        return 0;
                    }
                }
                s[j + 1] = s[j] - d_mat_entry(mu, i, j) * d_mat_entry(r, i, j);
            }
            d_mat_entry(r, i, i) = s[i];
            if (i >= newd && d_mat_entry(r, i, i) < d_gs_B) /* check removals */
            {
                d_mat_clear(r);
                d_mat_clear(mu);
                _d_vec_clear(s);
                fmpz_clear(dmax);
                return 0;
            }
            if (i < newd)
            {
                if ((fl->delta * d_mat_entry(r, i - 1, i - 1)) > s[i - 1])  /* check Lovasz condition */
                {
                    d_mat_clear(r);
                    d_mat_clear(mu);
                    _d_vec_clear(s);
                    fmpz_clear(dmax);
                    return 0;
                }
            }
        }
        d_mat_clear(r);
        d_mat_clear(mu);
        _d_vec_clear(s);
        fmpz_clear(dmax);
    }
    return 1;
}
