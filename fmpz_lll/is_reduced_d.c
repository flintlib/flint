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

#include <fenv.h>
#include "fmpz_lll.h"

int
fmpz_lll_is_reduced_d(const fmpz_mat_t B, const fmpz_lll_t fl)
{
    if (fl->rt == Z_BASIS)
    {
        slong i, j, k, m, n;
        d_mat_t A, Q, R, V, Wu, Wd, bound, bound2, bound3, boundt, mm, rm, mn,
            rn, absR;
        double *du, *dd;
        double s, norm = 0, ti, tj;
        int rounding_direction = fegetround();

        if (B->r == 0 || B->r == 1)
            return 1;

        m = B->c;
        n = B->r;

        d_mat_init(A, m, n);
        d_mat_init(Q, m, n);
        d_mat_init(R, n, n);
        d_mat_init(V, n, n);
        d_mat_zero(R);
        d_mat_zero(V);

        if (fmpz_mat_get_d_mat_transpose(A, B) == -1)
        {
            d_mat_clear(A);
            d_mat_clear(Q);
            d_mat_clear(R);
            d_mat_clear(V);
            return 0;
        }

        for (k = 0; k < n; k++)
        {
            for (j = 0; j < m; j++)
            {
                d_mat_entry(Q, j, k) = d_mat_entry(A, j, k);
            }
            for (i = 0; i < k; i++)
            {
                s = 0;
                for (j = 0; j < m; j++)
                {
                    s += d_mat_entry(Q, j, i) * d_mat_entry(Q, j, k);
                }
                d_mat_entry(R, i, k) = s;
                for (j = 0; j < m; j++)
                {
                    d_mat_entry(Q, j, k) -= s * d_mat_entry(Q, j, i);
                }
            }
            s = 0;
            for (j = 0; j < m; j++)
            {
                s += d_mat_entry(Q, j, k) * d_mat_entry(Q, j, k);
            }
            d_mat_entry(R, k, k) = s = sqrt(s);
            if (s != 0)
            {
                s = 1 / s;
                for (j = 0; j < m; j++)
                {
                    d_mat_entry(Q, j, k) *= s;
                }
            }
        }
        d_mat_clear(Q);

        for (j = n - 1; j >= 0; j--)
        {
            d_mat_entry(V, j, j) = 1.0 / d_mat_entry(R, j, j);
            for (i = j + 1; i < n; i++)
            {
                for (k = j + 1; k < n; k++)
                {
                    d_mat_entry(V, j, i) +=
                        d_mat_entry(V, k, i) * d_mat_entry(R, j, k);
                }
                d_mat_entry(V, j, i) *= -d_mat_entry(V, j, j);
            }
        }

        d_mat_init(Wu, n, n);
        d_mat_init(Wd, n, n);
        du = _d_vec_init(n);
        dd = _d_vec_init(n);

        fesetround(FE_DOWNWARD);
        d_mat_mul(Wd, R, V);
        for (i = 0; i < n; i++)
        {
            dd[i] = d_mat_entry(Wd, i, i) - 1;
        }
        fesetround(FE_UPWARD);
        d_mat_mul(Wu, R, V);
        for (i = 0; i < n; i++)
        {
            du[i] = d_mat_entry(Wu, i, i) - 1;
        }
        for (i = 0; i < n; i++)
        {
            s = 0;
            for (j = 0; j < n; j++)
            {
                if (i != j)
                    s += FLINT_MAX(fabs(d_mat_entry(Wd, i, j)),
                                   fabs(d_mat_entry(Wu, i, j)));
                else
                    s += FLINT_MAX(fabs(dd[i]), fabs(du[i]));
            }
            norm = FLINT_MAX(norm, s);
        }
        if (norm >= 1)
        {
            d_mat_clear(A);
            d_mat_clear(R);
            d_mat_clear(V);
            d_mat_clear(Wu);
            d_mat_clear(Wd);
            _d_vec_clear(du);
            _d_vec_clear(dd);
            fesetround(rounding_direction);
            return 0;
        }

        d_mat_init(bound, n, n);

        fesetround(FE_DOWNWARD);
        for (i = 0; i < n; i++)
        {
            dd[i] = d_mat_entry(Wd, i, i) - 2;
        }
        fesetround(FE_UPWARD);
        for (i = 0; i < n; i++)
        {
            du[i] = d_mat_entry(Wu, i, i) - 2;
        }
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                if (j > i)
                {
                    d_mat_entry(bound, i, j) =
                        FLINT_MAX(fabs(d_mat_entry(Wd, i, j)),
                                  fabs(d_mat_entry(Wu, i, j))) +
                        norm * norm / (1.0 - norm);
                }
                else if (j < i)
                {
                    d_mat_entry(bound, i, j) =
                        FLINT_MAX(fabs(d_mat_entry(Wd, i, j)),
                                  fabs(d_mat_entry(Wu, i, j)));
                }
                else
                {
                    d_mat_entry(bound, i, j) =
                        FLINT_MAX(fabs(dd[i]),
                                  fabs(du[i])) + norm * norm / (1.0 - norm);
                }
            }
        }
        _d_vec_clear(dd);
        _d_vec_clear(du);

        d_mat_init(mm, n, n);
        d_mat_init(rm, n, n);
        d_mat_init(mn, n, n);
        d_mat_init(rn, n, n);
        d_mat_init(bound2, n, n);

        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                d_mat_entry(mm, j, i) =
                    (d_mat_entry(Wu, i, j) + d_mat_entry(Wd, i, j)) / 2;
                d_mat_entry(rm, j, i) =
                    d_mat_entry(mm, j, i) - d_mat_entry(Wd, i, j);
                d_mat_entry(mn, i, j) =
                    (d_mat_entry(Wu, i, j) + d_mat_entry(Wd, i, j)) / 2;
                d_mat_entry(rn, i, j) =
                    d_mat_entry(mn, i, j) - d_mat_entry(Wd, i, j);
            }
        }
        fesetround(FE_DOWNWARD);
        d_mat_mul(Wd, mm, mn);
        for (i = 0; i < n; i++)
        {
            d_mat_entry(Wd, i, i) -= 1;
        }
        fesetround(FE_UPWARD);
        d_mat_mul(Wu, mm, mn);
        for (i = 0; i < n; i++)
        {
            d_mat_entry(Wu, i, i) -= 1;
            for (j = 0; j < n; j++)
            {
                d_mat_entry(Wu, i, j) =
                    FLINT_MAX(fabs(d_mat_entry(Wd, i, j)),
                              fabs(d_mat_entry(Wu, i, j)));
                d_mat_entry(mm, i, j) = fabs(d_mat_entry(mm, i, j));
                d_mat_entry(mn, i, j) = fabs(d_mat_entry(mn, i, j));
            }
        }
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                d_mat_entry(bound2, i, j) =
                    d_mat_entry(mn, i, j) + d_mat_entry(rn, i, j);
            }
        }
        d_mat_mul(bound2, rm, bound2);
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                d_mat_entry(bound2, i, j) += d_mat_entry(Wu, i, j);
            }
        }
        d_mat_mul(Wu, mm, rn);
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                d_mat_entry(bound2, i, j) += d_mat_entry(Wu, i, j);
            }
        }

        d_mat_clear(Wu);
        d_mat_clear(Wd);
        d_mat_clear(mm);
        d_mat_clear(mn);
        d_mat_clear(rm);
        d_mat_clear(rn);

        d_mat_init(Wu, m, n);
        d_mat_init(Wd, m, n);
        d_mat_init(mm, n, m);
        d_mat_init(mn, m, n);
        d_mat_init(rm, n, m);
        d_mat_init(rn, m, n);

        fesetround(FE_DOWNWARD);
        d_mat_mul(Wd, A, V);
        fesetround(FE_UPWARD);
        d_mat_mul(Wu, A, V);

        d_mat_clear(A);
        d_mat_clear(V);

        d_mat_init(bound3, n, n);

        for (i = 0; i < m; i++)
        {
            for (j = 0; j < n; j++)
            {
                d_mat_entry(mm, j, i) =
                    (d_mat_entry(Wu, i, j) + d_mat_entry(Wd, i, j)) / 2;
                d_mat_entry(rm, j, i) =
                    d_mat_entry(mm, j, i) - d_mat_entry(Wd, i, j);
                d_mat_entry(mn, i, j) =
                    (d_mat_entry(Wu, i, j) + d_mat_entry(Wd, i, j)) / 2;
                d_mat_entry(rn, i, j) =
                    d_mat_entry(mn, i, j) - d_mat_entry(Wd, i, j);
            }
        }

        d_mat_clear(Wd);
        d_mat_clear(Wu);

        d_mat_init(Wd, n, n);
        d_mat_init(Wu, n, n);

        fesetround(FE_DOWNWARD);
        d_mat_mul(Wd, mm, mn);
        for (i = 0; i < n; i++)
        {
            d_mat_entry(Wd, i, i) -= 1;
        }
        fesetround(FE_UPWARD);
        d_mat_mul(Wu, mm, mn);
        for (i = 0; i < n; i++)
        {
            d_mat_entry(Wu, i, i) -= 1;
            for (j = 0; j < m; j++)
            {
                if (j < n)
                {
                    d_mat_entry(Wu, i, j) =
                        FLINT_MAX(fabs(d_mat_entry(Wd, i, j)),
                                  fabs(d_mat_entry(Wu, i, j)));
                }
                d_mat_entry(mm, i, j) = fabs(d_mat_entry(mm, i, j));
            }
        }

        d_mat_clear(Wd);
        d_mat_init(Wd, m, n);

        for (i = 0; i < m; i++)
        {
            for (j = 0; j < n; j++)
            {
                d_mat_entry(Wd, i, j) =
                    fabs(d_mat_entry(mn, i, j)) + d_mat_entry(rn, i, j);
            }
        }
        d_mat_mul(bound3, rm, Wd);
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                d_mat_entry(bound3, i, j) += d_mat_entry(Wu, i, j);
            }
        }
        d_mat_mul(Wu, mm, rn);
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                d_mat_entry(bound3, i, j) += d_mat_entry(Wu, i, j);
            }
        }

        d_mat_clear(Wu);
        d_mat_clear(Wd);
        d_mat_clear(mm);
        d_mat_clear(mn);
        d_mat_clear(rm);
        d_mat_clear(rn);

        d_mat_init(boundt, n, n);

        d_mat_transpose(boundt, bound);
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                d_mat_entry(bound2, i, j) =
                    fabs(d_mat_entry(bound2, i, j)) +
                    fabs(d_mat_entry(bound3, i, j));
            }
        }
        d_mat_mul(bound, bound2, bound);
        d_mat_mul(bound, boundt, bound);

        d_mat_clear(bound2);
        d_mat_clear(bound3);
        d_mat_clear(boundt);

        norm = 0;
        for (i = 0; i < n; i++)
        {
            s = 0;
            for (j = 0; j < n; j++)
            {
                s += fabs(d_mat_entry(bound, i, j));
            }
            norm = FLINT_MAX(norm, s);
        }
        if (norm >= 1)
        {
            d_mat_clear(R);
            d_mat_clear(bound);
            fesetround(rounding_direction);
            return 0;
        }

        d_mat_init(absR, n, n);
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                if (j >= i)
                {
                    d_mat_entry(bound, i, j) += norm * norm / (1.0 - norm);
                }
                else
                {
                    d_mat_entry(bound, i, j) = 0;
                }
                d_mat_entry(absR, i, j) = fabs(d_mat_entry(R, i, j));
            }
        }
        d_mat_mul(bound, bound, absR);

        d_mat_clear(absR);

        for (i = 0; i < n - 1; i++)
        {
            fesetround(FE_DOWNWARD);
            ti = (d_mat_entry(R, i, i) - d_mat_entry(bound, i, i)) * fl->eta;
            fesetround(FE_UPWARD);
            for (j = i + 1; j < n; j++)
            {
                tj = fabs(d_mat_entry(R, i, j)) + d_mat_entry(bound, i, j);
                if (tj > ti)
                {
                    d_mat_clear(R);
                    d_mat_clear(bound);
                    return 0;
                }
            }
            ti = d_mat_entry(R, i, i) + d_mat_entry(bound, i, i);
            fesetround(FE_DOWNWARD);
            tj = d_mat_entry(R, i + 1, i + 1) - d_mat_entry(bound, i + 1,
                                                            i + 1);
            s = ((fabs(d_mat_entry(R, i, i + 1)) -
                  d_mat_entry(bound, i,
                              i + 1)) / ti) * ((fabs(d_mat_entry(R, i,
                                                                 i + 1)) -
                                                d_mat_entry(bound, i,
                                                            i + 1)) / ti) -
                fl->delta;
            s = -s;
            fesetround(FE_UPWARD);
            s = sqrt(s) * ti;
            if (s > tj)
            {
                d_mat_clear(R);
                d_mat_clear(bound);
                return 0;
            }
        }

        d_mat_clear(R);
        d_mat_clear(bound);
        fesetround(rounding_direction);

#if 0
        slong i, j, k, m, n;
        fmpz_mat_t Btrans;
        d_mat_t A;
        double *d;
        double s, fak, alpha, up, deltap, etap, thetap;
        double c1, c2, c3, c4, c5, c6, c7;

        alpha = 1 / sqrt(fl->deltsa - (fl->eta * fl->eta));

        if (B->r == 0 || B->r == 1)
            return 1;

        fmpz_mat_init(Btrans, B->c, B->r);
        fmpz_mat_transpose(Btrans, B);
        m = Btrans->r;
        n = Btrans->c;

        c1 = (sqrt(6) + sqrt(3)) * m * sqrt(n - 1);
        c2 = 2 * sqrt(1 + (n - 2) * fl->eta * fl->eta) / (1 + fl->eta);
        c3 = (fabs(1 - fl->eta) * alpha +
              1) * m * sqrt(n) / (((1 + fl->eta) * alpha - 1) * (sqrt(1.5) -
                                                                 1));
        c4 = FLINT_MAX(c3, sqrt(2) * c1 * c2);
        c5 = 4 * (6 * m + 63);
        c6 = 0.5 * n * c5;
        c7 = c4 * c6;
        up = c7 * pow((1 + fl->eta), n) * pow(alpha, n) * D_EPS;
        if (up >= 1)
        {
            fmpz_mat_clear(Btrans);
            return 0;
        }

        d_mat_init(A, m, n);

        d = _d_vec_init(n);

        if (fmpz_mat_get_d_mat(A, Btrans) == -1)
        {
            fmpz_mat_clear(Btrans);
            d_mat_clear(A);
            _d_vec_clear(d);
            return 0;
        }
        fmpz_mat_clear(Btrans);

        for (j = 0; j < n; j++)
        {
            s = 0.0;
            for (i = j; i < m; i++)
            {
                s += d_mat_entry(A, i, j) * d_mat_entry(A, i, j);
            }
            s = sqrt(s);
            d[j] = (d_mat_entry(A, j, j) > 0) ? (-s) : s;
            fak = sqrt(s * (s + fabs(d_mat_entry(A, j, j))));
            d_mat_entry(A, j, j) -= d[j];
            if (fak != 0.0)
            {
                for (k = j; k < m; k++)
                {
                    d_mat_entry(A, k, j) /= fak;
                }
                for (i = j + 1; i < n; i++)
                {
                    s = 0.0;
                    for (k = j; k < m; k++)
                    {
                        s += d_mat_entry(A, k, j) * d_mat_entry(A, k, i);
                    }
                    for (k = j; k < m; k++)
                    {
                        d_mat_entry(A, k, i) -= d_mat_entry(A, k, j) * s;
                    }
                }
            }
        }

        deltap =
            fl->delta * (1 - up) * (1 -
                                    up) / ((1 + up) * (1 + up) * (1 +
                                                                  2 * up *
                                                                  fl->eta *
                                                                  alpha));
        etap = fl->eta / (1 - up);
        thetap = up / (1 - up);

        for (i = 0; i < n - 1; i++)
        {
            for (j = i + 1; j < n; j++)
            {
                if (fabs(d_mat_entry(A, i, j)) >
                    (etap * fabs(d[i]) + thetap * fabs(d[j])))
                {
                    d_mat_clear(A);
                    _d_vec_clear(d);
                    return 0;
                }
            }
            if (deltap * d[i] * d[i] >
                (d_mat_entry(A, i, i + 1) * d_mat_entry(A, i, i + 1) +
                 d[i + 1] * d[i + 1]))
            {
                d_mat_clear(A);
                _d_vec_clear(d);
                return 0;
            }
        }

        d_mat_clear(A);
        _d_vec_clear(d);


        slong i, j, k, d = B->r, n = B->c;
        d_mat_t appB, Q, mu;

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

        for (j = 0; j < n; j++)
        {
            d_mat_entry(Q, 0, j) = d_mat_entry(appB, 0, j);
        }
        /* diagonal of mu stores the squared GS norms */
        d_mat_entry(mu, 0, 0) = _d_vec_norm(Q->rows[0], n);

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

                if (fabs(d_mat_entry(mu, i, j)) > fl->eta)  /* check size reduction */
                {
                    d_mat_clear(appB);
                    d_mat_clear(Q);
                    d_mat_clear(mu);
                    return 0;
                }
            }

            d_mat_entry(mu, i, i) = _d_vec_norm(Q->rows[i], n);
            if ((fl->delta - pow(d_mat_entry(mu, i, i - 1), 2)) * d_mat_entry(mu, i - 1, i - 1) > d_mat_entry(mu, i, i))    /* check Lovasz condition */
            {
                d_mat_clear(appB);
                d_mat_clear(Q);
                d_mat_clear(mu);
                return 0;
            }
        }
        d_mat_clear(appB);
        d_mat_clear(Q);
        d_mat_clear(mu);
#endif
    }
    else
    {
        slong i, j, k, d = B->r;
        d_mat_t r, mu;
        double *s;
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
                if (fabs(d_mat_entry(mu, i, j)) > fl->eta)  /* check size reduction */
                {
                    d_mat_clear(r);
                    d_mat_clear(mu);
                    _d_vec_clear(s);
                    fmpz_clear(dmax);
                    return 0;
                }
                s[j + 1] = s[j] - d_mat_entry(mu, i, j) * d_mat_entry(r, i, j);
            }
            d_mat_entry(r, i, i) = s[i];
            if ((fl->delta * d_mat_entry(r, i - 1, i - 1)) > s[i - 1])  /* check Lovasz condition */
            {
                d_mat_clear(r);
                d_mat_clear(mu);
                _d_vec_clear(s);
                fmpz_clear(dmax);
                return 0;
            }
        }
        d_mat_clear(r);
        d_mat_clear(mu);
        _d_vec_clear(s);
        fmpz_clear(dmax);
    }
    return 1;
}
