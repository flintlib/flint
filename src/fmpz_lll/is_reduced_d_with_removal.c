/*
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_lll.h"
#if FLINT_USES_FENV
#include <fenv.h>
#endif

int
fmpz_lll_is_reduced_d_with_removal(const fmpz_mat_t B, const fmpz_lll_t fl,
                                   const fmpz_t gs_B, int newd)
{
#if FLINT_USES_FENV
    if (fl->rt == Z_BASIS)
    {
        /* NOTE: this algorithm should *not* be changed */
        slong i, j, k, m, n;
        d_mat_t A, Q, R, V, Wu, Wd, bound, bound2, bound3, boundt, mm, rm, mn,
            rn, absR;
        double *du, *dd;
        double s, norm = 0, ti, tj, d_gs_B;
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

        d_gs_B = fmpz_get_d(gs_B);

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
        d_mat_mul_classical(Wd, R, V);
        for (i = 0; i < n; i++)
        {
            dd[i] = d_mat_entry(Wd, i, i) - 1;
        }
        fesetround(FE_UPWARD);
        d_mat_mul_classical(Wu, R, V);
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
        if (!(norm < 1))
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
        d_mat_mul_classical(Wd, mm, mn);
        for (i = 0; i < n; i++)
        {
            d_mat_entry(Wd, i, i) -= 1;
        }
        fesetround(FE_UPWARD);
        d_mat_mul_classical(Wu, mm, mn);
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
        d_mat_mul_classical(bound2, rm, bound2);
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                d_mat_entry(bound2, i, j) += d_mat_entry(Wu, i, j);
            }
        }
        d_mat_mul_classical(Wu, mm, rn);
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
        d_mat_mul_classical(Wd, A, V);
        fesetround(FE_UPWARD);
        d_mat_mul_classical(Wu, A, V);

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
        d_mat_mul_classical(Wd, mm, mn);
        for (i = 0; i < n; i++)
        {
            d_mat_entry(Wd, i, i) -= 1;
        }
        fesetround(FE_UPWARD);
        d_mat_mul_classical(Wu, mm, mn);
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
        d_mat_mul_classical(bound3, rm, Wd);
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                d_mat_entry(bound3, i, j) += d_mat_entry(Wu, i, j);
            }
        }
        d_mat_mul_classical(Wu, mm, rn);
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
        d_mat_mul_classical(bound, bound2, bound);
        d_mat_mul_classical(bound, boundt, bound);

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
        if (!(norm < 1))
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
        d_mat_mul_classical(bound, bound, absR);

        d_mat_clear(absR);

        for (i = 0; i < n - 1; i++)
        {
            fesetround(FE_DOWNWARD);
            ti = (s =
                  (d_mat_entry(R, i, i) - d_mat_entry(bound, i, i))) * fl->eta;
            if (i >= newd && !(s*s >= d_gs_B))
            {
                d_mat_clear(R);
                d_mat_clear(bound);
                fesetround(rounding_direction);
                return 0;
            }
            fesetround(FE_UPWARD);
            for (j = i + 1; j < n; j++)
            {
                tj = fabs(d_mat_entry(R, i, j)) + d_mat_entry(bound, i, j);
                if (i < newd && !(tj <= ti))
                {
                    d_mat_clear(R);
                    d_mat_clear(bound);
                    fesetround(rounding_direction);
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
            if (i < newd && !(s <= tj))
            {
                d_mat_clear(R);
                d_mat_clear(bound);
                fesetround(rounding_direction);
                return 0;
            }
        }
        fesetround(FE_DOWNWARD);
        s = (d_mat_entry(R, i, i) - d_mat_entry(bound, i, i));
        if (i >= newd && !(s*s >= d_gs_B))
        {
            d_mat_clear(R);
            d_mat_clear(bound);
            fesetround(rounding_direction);
            return 0;
        }

        d_mat_clear(R);
        d_mat_clear(bound);
        fesetround(rounding_direction);
    }
    else
    {
        /* NOTE: this algorithm should *not* be changed */
        slong i, j, k, m, n;
        d_mat_t A, R, V, Wu, Wd, bound, bound2, bound3, boundt, mm, rm, mn,
            rn, absR;
        double *du, *dd;
        double s, norm = 0, ti, tj, d_gs_B;
        int rounding_direction = fegetround();

        if (B->r == 0 || B->r == 1)
            return 1;

        m = B->c;
        n = B->r;

        d_mat_init(A, m, n);
        d_mat_init(R, n, n);
        d_mat_init(V, n, n);
        d_mat_zero(R);
        d_mat_zero(V);

        if (fmpz_mat_get_d_mat_transpose(A, B) == -1)
        {
            d_mat_clear(A);
            d_mat_clear(R);
            d_mat_clear(V);
            return 0;
        }

        for (j = 0; j < n; j++)
        {
            d_mat_entry(R, j, j) = d_mat_entry(A, j, j);
            for (i = 0; i <= j - 1; i++)
            {
                d_mat_entry(R, i, j) = d_mat_entry(A, j, i);
                for (k = 0; k <= i - 1; k++)
                {
                    d_mat_entry(R, i, j) -=
                        d_mat_entry(R, k, i) * d_mat_entry(R, k, j);
                }
                if (d_mat_entry(R, i, i) != 0)
                {
                    d_mat_entry(R, i, j) /= d_mat_entry(R, i, i);
                    d_mat_entry(R, j, j) -=
                        d_mat_entry(R, i, j) * d_mat_entry(R, i, j);
                }
            }

            if (!(d_mat_entry(R, j, j) > 0))
            {
                /* going to take sqrt and then divide by it */
                d_mat_clear(A);
                d_mat_clear(R);
                d_mat_clear(V);
                return 0;
            }

            d_mat_entry(R, j, j) = sqrt(d_mat_entry(R, j, j));
        }

        d_gs_B = fmpz_get_d(gs_B);

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
        d_mat_mul_classical(Wd, R, V);
        for (i = 0; i < n; i++)
        {
            dd[i] = d_mat_entry(Wd, i, i) - 1;
        }
        fesetround(FE_UPWARD);
        d_mat_mul_classical(Wu, R, V);
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
        if (!(norm < 1))
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
        d_mat_mul_classical(Wd, mm, mn);
        for (i = 0; i < n; i++)
        {
            d_mat_entry(Wd, i, i) -= 1;
        }
        fesetround(FE_UPWARD);
        d_mat_mul_classical(Wu, mm, mn);
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
        d_mat_mul_classical(bound2, rm, bound2);
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                d_mat_entry(bound2, i, j) += d_mat_entry(Wu, i, j);
            }
        }
        d_mat_mul_classical(Wu, mm, rn);
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

        d_mat_transpose(mm, V);
        fesetround(FE_DOWNWARD);
        d_mat_mul_classical(Wd, mm, A);
        fesetround(FE_UPWARD);
        d_mat_mul_classical(Wu, mm, A);

        d_mat_clear(A);

        d_mat_init(bound3, n, n);

        fesetround(FE_DOWNWARD);
        d_mat_mul_classical(mm, Wd, V);
        for (i = 0; i < n; i++)
        {
            d_mat_entry(mm, i, i) -= 1;
        }
        fesetround(FE_UPWARD);
        d_mat_mul_classical(rm, Wd, V);
        for (i = 0; i < n; i++)
        {
            d_mat_entry(rm, i, i) -= 1;
        }

        fesetround(FE_DOWNWARD);
        d_mat_mul_classical(mn, Wu, V);
        for (i = 0; i < n; i++)
        {
            d_mat_entry(mn, i, i) -= 1;
        }
        fesetround(FE_UPWARD);
        d_mat_mul_classical(rn, Wu, V);
        for (i = 0; i < n; i++)
        {
            d_mat_entry(rn, i, i) -= 1;
        }

        d_mat_clear(Wd);
        d_mat_clear(Wu);
        d_mat_clear(V);

        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                d_mat_entry(bound3, i, j) =
                    FLINT_MAX(fabs(d_mat_entry(mm, i, j)),
                              fabs(d_mat_entry(mn, i, j)));
                d_mat_entry(bound3, i, j) =
                    FLINT_MAX(fabs(d_mat_entry(bound3, i, j)),
                              fabs(d_mat_entry(rm, i, j)));
                d_mat_entry(bound3, i, j) =
                    FLINT_MAX(fabs(d_mat_entry(bound3, i, j)),
                              fabs(d_mat_entry(rn, i, j)));
            }
        }

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
        d_mat_mul_classical(bound, bound2, bound);
        d_mat_mul_classical(bound, boundt, bound);

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
        if (!(norm < 1))
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
        d_mat_mul_classical(bound, bound, absR);

        d_mat_clear(absR);

        for (i = 0; i < n - 1; i++)
        {
            fesetround(FE_DOWNWARD);
            ti = (s =
                  (d_mat_entry(R, i, i) - d_mat_entry(bound, i, i))) * fl->eta;
            if (i >= newd && !(s*s >= d_gs_B))
            {
                d_mat_clear(R);
                d_mat_clear(bound);
                fesetround(rounding_direction);
                return 0;
            }
            fesetround(FE_UPWARD);
            for (j = i + 1; j < n; j++)
            {
                tj = fabs(d_mat_entry(R, i, j)) + d_mat_entry(bound, i, j);
                if (i < newd && !(tj <= ti))
                {
                    d_mat_clear(R);
                    d_mat_clear(bound);
                    fesetround(rounding_direction);
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
            if (i < newd && !(s <= tj))
            {
                d_mat_clear(R);
                d_mat_clear(bound);
                fesetround(rounding_direction);
                return 0;
            }
        }
        fesetround(FE_DOWNWARD);
        s = (d_mat_entry(R, i, i) - d_mat_entry(bound, i, i));
        if (i >= newd && !(s*s >= d_gs_B))
        {
            d_mat_clear(R);
            d_mat_clear(bound);
            fesetround(rounding_direction);
            return 0;
        }

        d_mat_clear(R);
        d_mat_clear(bound);
        fesetround(rounding_direction);
    }

    FLINT_ASSERT((fl->rt == Z_BASIS
          ? fmpz_mat_is_reduced_with_removal(B, fl->delta, fl->eta, gs_B, newd)
          : fmpz_mat_is_reduced_gram_with_removal(B, fl->delta, fl->eta, gs_B, newd)));

    return 1;
#else
    return 0;
#endif
}
