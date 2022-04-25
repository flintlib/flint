/*
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_lll.h"
#include "mpfr_vec.h"
#include "mpfr_mat.h"

int
fmpz_lll_is_reduced_mpfr_with_removal(const fmpz_mat_t B, const fmpz_lll_t fl,
                                      const fmpz_t gs_B, int newd,
                                      flint_bitcnt_t prec)
{
    if (fl->rt == Z_BASIS)
    {
        /* NOTE: this algorithm should *not* be changed */
        slong i, j, k, m, n;
        mpfr_mat_t A, Q, R, V, Wu, Wd, bound, bound2, bound3, boundt, mm, rm,
            mn, rn, absR;
        flint_mpfr *du, *dd;
        mpfr_t s, norm, ti, tj, tmp, mpfr_gs_B;

        if (B->r == 0 || B->r == 1)
            return 1;

        m = B->c;
        n = B->r;

        mpfr_mat_init(A, m, n, prec);
        mpfr_mat_init(Q, m, n, prec);
        mpfr_mat_init(R, n, n, prec);
        mpfr_mat_init(V, n, n, prec);
        mpfr_mat_zero(R);
        mpfr_mat_zero(V);

        mpfr_inits2(prec, s, norm, ti, tj, tmp, NULL);

        for (i = 0; i < n; i++)
        {
            for (j = 0; j < m; j++)
            {
                fmpz_get_mpfr(mpfr_mat_entry(A, j, i), fmpz_mat_entry(B, i, j),
                              MPFR_RNDN);
            }
        }

        for (k = 0; k < n; k++)
        {
            for (j = 0; j < m; j++)
            {
                mpfr_set(mpfr_mat_entry(Q, j, k), mpfr_mat_entry(A, j, k),
                         MPFR_RNDN);
            }
            for (i = 0; i < k; i++)
            {
                mpfr_set_zero(s, 1);
                for (j = 0; j < m; j++)
                {
                    mpfr_mul(norm, mpfr_mat_entry(Q, j, i),
                             mpfr_mat_entry(Q, j, k), MPFR_RNDN);
                    mpfr_add(s, s, norm, MPFR_RNDN);
                }
                mpfr_set(mpfr_mat_entry(R, i, k), s, MPFR_RNDN);
                for (j = 0; j < m; j++)
                {
                    mpfr_mul(norm, s, mpfr_mat_entry(Q, j, i), MPFR_RNDN);
                    mpfr_sub(mpfr_mat_entry(Q, j, k), mpfr_mat_entry(Q, j, k),
                             norm, MPFR_RNDN);
                }
            }
            mpfr_set_zero(s, 1);
            for (j = 0; j < m; j++)
            {
                mpfr_mul(norm, mpfr_mat_entry(Q, j, k),
                         mpfr_mat_entry(Q, j, k), MPFR_RNDN);
                mpfr_add(s, s, norm, MPFR_RNDN);
            }
            mpfr_sqrt(s, s, MPFR_RNDN);
            mpfr_set(mpfr_mat_entry(R, k, k), s, MPFR_RNDN);
            if (!mpfr_zero_p(s))
            {
                mpfr_ui_div(s, 1, s, MPFR_RNDN);
                for (j = 0; j < m; j++)
                {
                    mpfr_mul(mpfr_mat_entry(Q, j, k), mpfr_mat_entry(Q, j, k),
                             s, MPFR_RNDN);
                }
            }
        }
        mpfr_mat_clear(Q);

        for (j = n - 1; j >= 0; j--)
        {
            mpfr_ui_div(mpfr_mat_entry(V, j, j), 1, mpfr_mat_entry(R, j, j),
                        MPFR_RNDN);
            for (i = j + 1; i < n; i++)
            {
                for (k = j + 1; k < n; k++)
                {
                    mpfr_mul(norm, mpfr_mat_entry(V, k, i),
                             mpfr_mat_entry(R, j, k), MPFR_RNDN);
                    mpfr_add(mpfr_mat_entry(V, j, i), mpfr_mat_entry(V, j, i),
                             norm, MPFR_RNDN);
                }
                mpfr_mul_si(mpfr_mat_entry(V, j, i), mpfr_mat_entry(V, j, i),
                            -WORD(1), MPFR_RNDN);
                mpfr_mul(mpfr_mat_entry(V, j, i), mpfr_mat_entry(V, j, i),
                         mpfr_mat_entry(V, j, j), MPFR_RNDN);
            }
        }

        mpfr_mat_init(Wu, n, n, prec);
        mpfr_mat_init(Wd, n, n, prec);
        du = _mpfr_vec_init(n, prec);
        dd = _mpfr_vec_init(n, prec);

        mpfr_mat_mul_classical(Wd, R, V, MPFR_RNDD);
        for (i = 0; i < n; i++)
        {
            mpfr_sub_ui(dd + i, mpfr_mat_entry(Wd, i, i), 1, MPFR_RNDD);
        }
        mpfr_mat_mul_classical(Wu, R, V, MPFR_RNDU);
        for (i = 0; i < n; i++)
        {
            mpfr_sub_ui(du + i, mpfr_mat_entry(Wu, i, i), 1, MPFR_RNDU);
        }
        mpfr_set_zero(norm, 1);
        for (i = 0; i < n; i++)
        {
            mpfr_set_zero(s, 1);
            for (j = 0; j < n; j++)
            {
                if (i != j)
                {
                    mpfr_abs(ti, mpfr_mat_entry(Wd, i, j), MPFR_RNDU);
                    mpfr_abs(tj, mpfr_mat_entry(Wu, i, j), MPFR_RNDU);
                    mpfr_max(tmp, ti, tj, MPFR_RNDU);
                    mpfr_add(s, s, tmp, MPFR_RNDU);
                }
                else
                {
                    mpfr_abs(ti, dd + i, MPFR_RNDU);
                    mpfr_abs(tj, du + i, MPFR_RNDU);
                    mpfr_max(tmp, ti, tj, MPFR_RNDU);
                    mpfr_add(s, s, tmp, MPFR_RNDU);
                }
            }
            mpfr_max(norm, norm, s, MPFR_RNDU);
        }
        if (mpfr_cmp_ui(norm, 1) >= 0)
        {
            mpfr_mat_clear(A);
            mpfr_mat_clear(R);
            mpfr_mat_clear(V);
            mpfr_mat_clear(Wu);
            mpfr_mat_clear(Wd);
            _mpfr_vec_clear(du, n);
            _mpfr_vec_clear(dd, n);
            mpfr_clears(s, norm, ti, tj, tmp, NULL);
            return 0;
        }

        mpfr_mat_init(bound, n, n, prec);

        for (i = 0; i < n; i++)
        {
            mpfr_sub_ui(dd + i, mpfr_mat_entry(Wd, i, i), 2, MPFR_RNDD);
        }
        for (i = 0; i < n; i++)
        {
            mpfr_sub_ui(du + i, mpfr_mat_entry(Wu, i, i), 2, MPFR_RNDU);
        }
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                if (j > i)
                {
                    mpfr_abs(ti, mpfr_mat_entry(Wd, i, j), MPFR_RNDU);
                    mpfr_abs(tj, mpfr_mat_entry(Wu, i, j), MPFR_RNDU);
                    mpfr_max(mpfr_mat_entry(bound, i, j), ti, tj, MPFR_RNDU);
                    mpfr_mul(ti, norm, norm, MPFR_RNDU);
                    mpfr_ui_sub(tj, 1, norm, MPFR_RNDU);
                    mpfr_div(tmp, ti, tj, MPFR_RNDU);
                    mpfr_add(mpfr_mat_entry(bound, i, j),
                             mpfr_mat_entry(bound, i, j), tmp, MPFR_RNDU);
                }
                else if (j < i)
                {
                    mpfr_abs(ti, mpfr_mat_entry(Wd, i, j), MPFR_RNDU);
                    mpfr_abs(tj, mpfr_mat_entry(Wu, i, j), MPFR_RNDU);
                    mpfr_max(mpfr_mat_entry(bound, i, j), ti, tj, MPFR_RNDU);
                }
                else
                {
                    mpfr_abs(ti, dd + i, MPFR_RNDU);
                    mpfr_abs(tj, du + i, MPFR_RNDU);
                    mpfr_max(mpfr_mat_entry(bound, i, j), ti, tj, MPFR_RNDU);
                    mpfr_mul(ti, norm, norm, MPFR_RNDU);
                    mpfr_ui_sub(tj, 1, norm, MPFR_RNDU);
                    mpfr_div(tmp, ti, tj, MPFR_RNDU);
                    mpfr_add(mpfr_mat_entry(bound, i, j),
                             mpfr_mat_entry(bound, i, j), tmp, MPFR_RNDU);
                }
            }
        }
        _mpfr_vec_clear(dd, n);
        _mpfr_vec_clear(du, n);

        mpfr_mat_init(mm, n, n, prec);
        mpfr_mat_init(rm, n, n, prec);
        mpfr_mat_init(mn, n, n, prec);
        mpfr_mat_init(rn, n, n, prec);
        mpfr_mat_init(bound2, n, n, prec);

        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                mpfr_add(tmp, mpfr_mat_entry(Wu, i, j),
                         mpfr_mat_entry(Wd, i, j), MPFR_RNDU);
                mpfr_div_ui(mpfr_mat_entry(mm, j, i), tmp, 2, MPFR_RNDU);
                mpfr_sub(mpfr_mat_entry(rm, j, i), mpfr_mat_entry(mm, j, i),
                         mpfr_mat_entry(Wd, i, j), MPFR_RNDU);
                mpfr_div_ui(mpfr_mat_entry(mn, i, j), tmp, 2, MPFR_RNDU);
                mpfr_sub(mpfr_mat_entry(rn, i, j), mpfr_mat_entry(mn, i, j),
                         mpfr_mat_entry(Wd, i, j), MPFR_RNDU);
            }
        }
        mpfr_mat_mul_classical(Wd, mm, mn, MPFR_RNDD);
        for (i = 0; i < n; i++)
        {
            mpfr_sub_ui(mpfr_mat_entry(Wd, i, i), mpfr_mat_entry(Wd, i, i), 1,
                        MPFR_RNDD);
        }
        mpfr_mat_mul_classical(Wu, mm, mn, MPFR_RNDU);
        for (i = 0; i < n; i++)
        {
            mpfr_sub_ui(mpfr_mat_entry(Wu, i, i), mpfr_mat_entry(Wu, i, i), 1,
                        MPFR_RNDU);
            for (j = 0; j < n; j++)
            {
                mpfr_abs(ti, mpfr_mat_entry(Wd, i, j), MPFR_RNDU);
                mpfr_abs(tj, mpfr_mat_entry(Wu, i, j), MPFR_RNDU);
                mpfr_max(mpfr_mat_entry(Wu, i, j), ti, tj, MPFR_RNDU);
                mpfr_abs(mpfr_mat_entry(mm, i, j), mpfr_mat_entry(mm, i, j),
                         MPFR_RNDU);
                mpfr_abs(mpfr_mat_entry(mn, i, j), mpfr_mat_entry(mn, i, j),
                         MPFR_RNDU);
            }
        }
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                mpfr_add(mpfr_mat_entry(bound2, i, j),
                         mpfr_mat_entry(mn, i, j), mpfr_mat_entry(rn, i, j),
                         MPFR_RNDU);
            }
        }
        mpfr_mat_mul_classical(bound2, rm, bound2, MPFR_RNDU);
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                mpfr_add(mpfr_mat_entry(bound2, i, j),
                         mpfr_mat_entry(bound2, i, j), mpfr_mat_entry(Wu, i,
                                                                      j),
                         MPFR_RNDU);
            }
        }
        mpfr_mat_mul_classical(Wu, mm, rn, MPFR_RNDU);
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                mpfr_add(mpfr_mat_entry(bound2, i, j),
                         mpfr_mat_entry(bound2, i, j), mpfr_mat_entry(Wu, i,
                                                                      j),
                         MPFR_RNDU);
            }
        }

        mpfr_mat_clear(Wu);
        mpfr_mat_clear(Wd);
        mpfr_mat_clear(mm);
        mpfr_mat_clear(mn);
        mpfr_mat_clear(rm);
        mpfr_mat_clear(rn);

        mpfr_mat_init(Wu, m, n, prec);
        mpfr_mat_init(Wd, m, n, prec);
        mpfr_mat_init(mm, n, m, prec);
        mpfr_mat_init(mn, m, n, prec);
        mpfr_mat_init(rm, n, m, prec);
        mpfr_mat_init(rn, m, n, prec);

        mpfr_mat_mul_classical(Wd, A, V, MPFR_RNDD);
        mpfr_mat_mul_classical(Wu, A, V, MPFR_RNDU);

        mpfr_mat_clear(A);
        mpfr_mat_clear(V);

        mpfr_mat_init(bound3, n, n, prec);

        for (i = 0; i < m; i++)
        {
            for (j = 0; j < n; j++)
            {
                mpfr_add(tmp, mpfr_mat_entry(Wu, i, j),
                         mpfr_mat_entry(Wd, i, j), MPFR_RNDU);
                mpfr_div_ui(mpfr_mat_entry(mm, j, i), tmp, 2, MPFR_RNDU);
                mpfr_sub(mpfr_mat_entry(rm, j, i), mpfr_mat_entry(mm, j, i),
                         mpfr_mat_entry(Wd, i, j), MPFR_RNDU);
                mpfr_div_ui(mpfr_mat_entry(mn, i, j), tmp, 2, MPFR_RNDU);
                mpfr_sub(mpfr_mat_entry(rn, i, j), mpfr_mat_entry(mn, i, j),
                         mpfr_mat_entry(Wd, i, j), MPFR_RNDU);
            }
        }

        mpfr_mat_clear(Wd);
        mpfr_mat_clear(Wu);

        mpfr_mat_init(Wd, n, n, prec);
        mpfr_mat_init(Wu, n, n, prec);

        mpfr_mat_mul_classical(Wd, mm, mn, MPFR_RNDD);
        for (i = 0; i < n; i++)
        {
            mpfr_sub_ui(mpfr_mat_entry(Wd, i, i), mpfr_mat_entry(Wd, i, i), 1,
                        MPFR_RNDD);
        }
        mpfr_mat_mul_classical(Wu, mm, mn, MPFR_RNDU);
        for (i = 0; i < n; i++)
        {
            mpfr_sub_ui(mpfr_mat_entry(Wu, i, i), mpfr_mat_entry(Wu, i, i), 1,
                        MPFR_RNDU);
            for (j = 0; j < m; j++)
            {
                if (j < n)
                {
                    mpfr_abs(ti, mpfr_mat_entry(Wd, i, j), MPFR_RNDU);
                    mpfr_abs(tj, mpfr_mat_entry(Wu, i, j), MPFR_RNDU);
                    mpfr_max(mpfr_mat_entry(Wu, i, j), ti, tj, MPFR_RNDU);
                }
                mpfr_abs(mpfr_mat_entry(mm, i, j), mpfr_mat_entry(mm, i, j),
                         MPFR_RNDU);
            }
        }

        mpfr_mat_clear(Wd);
        mpfr_mat_init(Wd, m, n, prec);

        for (i = 0; i < m; i++)
        {
            for (j = 0; j < n; j++)
            {
                mpfr_abs(mpfr_mat_entry(mn, i, j), mpfr_mat_entry(mn, i, j),
                         MPFR_RNDU);
                mpfr_add(mpfr_mat_entry(Wd, i, j), mpfr_mat_entry(mn, i, j),
                         mpfr_mat_entry(rn, i, j), MPFR_RNDU);
            }
        }
        mpfr_mat_mul_classical(bound3, rm, Wd, MPFR_RNDU);
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                mpfr_add(mpfr_mat_entry(bound3, i, j),
                         mpfr_mat_entry(bound3, i, j), mpfr_mat_entry(Wu, i,
                                                                      j),
                         MPFR_RNDU);
            }
        }
        mpfr_mat_mul_classical(Wu, mm, rn, MPFR_RNDU);
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                mpfr_add(mpfr_mat_entry(bound3, i, j),
                         mpfr_mat_entry(bound3, i, j), mpfr_mat_entry(Wu, i,
                                                                      j),
                         MPFR_RNDU);
            }
        }

        mpfr_mat_clear(Wu);
        mpfr_mat_clear(Wd);
        mpfr_mat_clear(mm);
        mpfr_mat_clear(mn);
        mpfr_mat_clear(rm);
        mpfr_mat_clear(rn);

        mpfr_mat_init(boundt, n, n, prec);

        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                mpfr_set(mpfr_mat_entry(boundt, j, i),
                         mpfr_mat_entry(bound, i, j), MPFR_RNDU);
                mpfr_set(ti, mpfr_mat_entry(bound2, i, j), MPFR_RNDU);
                mpfr_set(tj, mpfr_mat_entry(bound3, i, j), MPFR_RNDU);
                mpfr_add(mpfr_mat_entry(bound2, i, j), ti, tj, MPFR_RNDU);
            }
        }
        mpfr_mat_mul_classical(bound, bound2, bound, MPFR_RNDU);
        mpfr_mat_mul_classical(bound, boundt, bound, MPFR_RNDU);

        mpfr_mat_clear(bound2);
        mpfr_mat_clear(bound3);
        mpfr_mat_clear(boundt);

        mpfr_set_zero(norm, 1);
        for (i = 0; i < n; i++)
        {
            mpfr_set_zero(s, 1);
            for (j = 0; j < n; j++)
            {
                mpfr_abs(tmp, mpfr_mat_entry(bound, i, j), MPFR_RNDU);
                mpfr_add(s, s, tmp, MPFR_RNDU);
            }
            mpfr_max(norm, norm, s, MPFR_RNDU);
        }
        if (mpfr_cmp_ui(norm, 1) >= 0)
        {
            mpfr_mat_clear(R);
            mpfr_mat_clear(bound);
            mpfr_clears(s, norm, ti, tj, tmp, NULL);
            return 0;
        }

        mpfr_mat_init(absR, n, n, prec);
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                if (j >= i)
                {
                    mpfr_mul(ti, norm, norm, MPFR_RNDU);
                    mpfr_ui_sub(tj, 1, norm, MPFR_RNDU);
                    mpfr_div(tmp, ti, tj, MPFR_RNDU);
                    mpfr_add(mpfr_mat_entry(bound, i, j),
                             mpfr_mat_entry(bound, i, j), tmp, MPFR_RNDU);
                }
                else
                {
                    mpfr_set_zero(mpfr_mat_entry(bound, i, j), 1);
                }
                mpfr_abs(mpfr_mat_entry(absR, i, j), mpfr_mat_entry(R, i, j),
                         MPFR_RNDU);
            }
        }
        mpfr_mat_mul_classical(bound, bound, absR, MPFR_RNDU);

        mpfr_mat_clear(absR);

        mpfr_init2(mpfr_gs_B, prec);
        fmpz_get_mpfr(mpfr_gs_B, gs_B, MPFR_RNDN);
        for (i = 0; i < n - 1; i++)
        {
            mpfr_sub(tmp, mpfr_mat_entry(R, i, i), mpfr_mat_entry(bound, i, i),
                     MPFR_RNDD);
            mpfr_mul_d(ti, tmp, fl->eta, MPFR_RNDD);
            mpfr_sqr(tmp, tmp, MPFR_RNDD);
            if (i >= newd && !mpfr_greaterequal_p(tmp, mpfr_gs_B))
            {
                mpfr_mat_clear(R);
                mpfr_mat_clear(bound);
                mpfr_clears(s, norm, ti, tj, tmp, mpfr_gs_B, NULL);
                return 0;
            }
            for (j = i + 1; j < n; j++)
            {
                mpfr_abs(tmp, mpfr_mat_entry(R, i, j), MPFR_RNDU);
                mpfr_add(tj, tmp, mpfr_mat_entry(bound, i, j), MPFR_RNDU);
                if (i < newd && !mpfr_lessequal_p(tj, ti))
                {
                    mpfr_mat_clear(R);
                    mpfr_mat_clear(bound);
                    mpfr_clears(s, norm, ti, tj, tmp, mpfr_gs_B, NULL);
                    return 0;
                }
            }
            mpfr_add(ti, mpfr_mat_entry(R, i, i), mpfr_mat_entry(bound, i, i),
                     MPFR_RNDU);
            mpfr_sub(tj, mpfr_mat_entry(R, i + 1, i + 1),
                     mpfr_mat_entry(bound, i + 1, i + 1), MPFR_RNDD);
            mpfr_abs(tmp, mpfr_mat_entry(R, i, i + 1), MPFR_RNDD);
            mpfr_sub(norm, tmp, mpfr_mat_entry(bound, i, i + 1), MPFR_RNDD);
            mpfr_div(tmp, norm, ti, MPFR_RNDD);
            mpfr_sqr(norm, tmp, MPFR_RNDD);
            mpfr_sub_d(s, norm, fl->delta, MPFR_RNDD);
            mpfr_neg(s, s, MPFR_RNDD);
            mpfr_sqrt(tmp, s, MPFR_RNDU);
            mpfr_mul(s, tmp, ti, MPFR_RNDU);
            if (i < newd && !mpfr_lessequal_p(s, tj))
            {
                mpfr_mat_clear(R);
                mpfr_mat_clear(bound);
                mpfr_clears(s, norm, ti, tj, tmp, mpfr_gs_B, NULL);
                return 0;
            }
        }
        mpfr_sub(tmp, mpfr_mat_entry(R, i, i), mpfr_mat_entry(bound, i, i),
                 MPFR_RNDD);
        mpfr_sqr(tmp, tmp, MPFR_RNDD);
        if (i >= newd && !mpfr_greaterequal_p(tmp, mpfr_gs_B))
        {
            mpfr_mat_clear(R);
            mpfr_mat_clear(bound);
            mpfr_clears(s, norm, ti, tj, tmp, mpfr_gs_B, NULL);
            return 0;
        }

        mpfr_mat_clear(R);
        mpfr_mat_clear(bound);
        mpfr_clears(s, norm, ti, tj, tmp, mpfr_gs_B, NULL);
    }
    else
    {
        slong i, j, k, m, n;
        mpfr_mat_t A, R, V, Wu, Wd, bound, bound2, bound3, boundt, mm, rm,
            mn, rn, absR;
        flint_mpfr *du, *dd;
        mpfr_t s, norm, ti, tj, tmp, mpfr_gs_B;

        if (B->r == 0 || B->r == 1)
            return 1;

        m = B->c;
        n = B->r;

        mpfr_mat_init(A, m, n, prec);
        mpfr_mat_init(R, n, n, prec);
        mpfr_mat_init(V, n, n, prec);
        mpfr_mat_zero(R);
        mpfr_mat_zero(V);

        mpfr_inits2(prec, s, norm, ti, tj, tmp, NULL);

        for (i = 0; i < n; i++)
        {
            for (j = 0; j < m; j++)
            {
                fmpz_get_mpfr(mpfr_mat_entry(A, j, i), fmpz_mat_entry(B, i, j),
                              MPFR_RNDN);
            }
        }

        for (j = 0; j < n; j++)
        {
            mpfr_set(mpfr_mat_entry(R, j, j), mpfr_mat_entry(A, j, j),
                     MPFR_RNDN);
            for (i = 0; i <= j - 1; i++)
            {
                mpfr_set(mpfr_mat_entry(R, i, j), mpfr_mat_entry(A, j, i),
                         MPFR_RNDN);
                for (k = 0; k <= i - 1; k++)
                {
                    mpfr_mul(tmp, mpfr_mat_entry(R, k, i),
                             mpfr_mat_entry(R, k, j), MPFR_RNDN);
                    mpfr_sub(mpfr_mat_entry(R, i, j), mpfr_mat_entry(R, i, j),
                             tmp, MPFR_RNDN);
                }
                if (!mpfr_zero_p(mpfr_mat_entry(R, i, i)))
                {
                    mpfr_div(mpfr_mat_entry(R, i, j), mpfr_mat_entry(R, i, j),
                             mpfr_mat_entry(R, i, i), MPFR_RNDN);
                    mpfr_mul(tmp, mpfr_mat_entry(R, i, j),
                             mpfr_mat_entry(R, i, j), MPFR_RNDN);
                    mpfr_sub(mpfr_mat_entry(R, j, j), mpfr_mat_entry(R, j, j),
                             tmp, MPFR_RNDN);
                }
            }

            if (mpfr_sgn(mpfr_mat_entry(R, j, j)) <= 0)
            {
                /* going to take sqrt and then divide by it */
                mpfr_mat_clear(A);
                mpfr_mat_clear(R);
                mpfr_mat_clear(V);
                return 0;
            }

            mpfr_sqrt(mpfr_mat_entry(R, j, j), mpfr_mat_entry(R, j, j),
                      MPFR_RNDN);
        }

        for (j = n - 1; j >= 0; j--)
        {
            mpfr_ui_div(mpfr_mat_entry(V, j, j), 1, mpfr_mat_entry(R, j, j),
                        MPFR_RNDN);
            for (i = j + 1; i < n; i++)
            {
                for (k = j + 1; k < n; k++)
                {
                    mpfr_mul(norm, mpfr_mat_entry(V, k, i),
                             mpfr_mat_entry(R, j, k), MPFR_RNDN);
                    mpfr_add(mpfr_mat_entry(V, j, i), mpfr_mat_entry(V, j, i),
                             norm, MPFR_RNDN);
                }
                mpfr_mul_si(mpfr_mat_entry(V, j, i), mpfr_mat_entry(V, j, i),
                            -WORD(1), MPFR_RNDN);
                mpfr_mul(mpfr_mat_entry(V, j, i), mpfr_mat_entry(V, j, i),
                         mpfr_mat_entry(V, j, j), MPFR_RNDN);
            }
        }

        mpfr_mat_init(Wu, n, n, prec);
        mpfr_mat_init(Wd, n, n, prec);
        du = _mpfr_vec_init(n, prec);
        dd = _mpfr_vec_init(n, prec);

        mpfr_mat_mul_classical(Wd, R, V, MPFR_RNDD);
        for (i = 0; i < n; i++)
        {
            mpfr_sub_ui(dd + i, mpfr_mat_entry(Wd, i, i), 1, MPFR_RNDD);
        }
        mpfr_mat_mul_classical(Wu, R, V, MPFR_RNDU);
        for (i = 0; i < n; i++)
        {
            mpfr_sub_ui(du + i, mpfr_mat_entry(Wu, i, i), 1, MPFR_RNDU);
        }
        mpfr_set_zero(norm, 1);
        for (i = 0; i < n; i++)
        {
            mpfr_set_zero(s, 1);
            for (j = 0; j < n; j++)
            {
                if (i != j)
                {
                    mpfr_abs(ti, mpfr_mat_entry(Wd, i, j), MPFR_RNDU);
                    mpfr_abs(tj, mpfr_mat_entry(Wu, i, j), MPFR_RNDU);
                    mpfr_max(tmp, ti, tj, MPFR_RNDU);
                    mpfr_add(s, s, tmp, MPFR_RNDU);
                }
                else
                {
                    mpfr_abs(ti, dd + i, MPFR_RNDU);
                    mpfr_abs(tj, du + i, MPFR_RNDU);
                    mpfr_max(tmp, ti, tj, MPFR_RNDU);
                    mpfr_add(s, s, tmp, MPFR_RNDU);
                }
            }
            mpfr_max(norm, norm, s, MPFR_RNDU);
        }
        if (mpfr_cmp_ui(norm, 1) >= 0)
        {
            mpfr_mat_clear(A);
            mpfr_mat_clear(R);
            mpfr_mat_clear(V);
            mpfr_mat_clear(Wu);
            mpfr_mat_clear(Wd);
            _mpfr_vec_clear(du, n);
            _mpfr_vec_clear(dd, n);
            mpfr_clears(s, norm, ti, tj, tmp, NULL);
            return 0;
        }

        mpfr_mat_init(bound, n, n, prec);

        for (i = 0; i < n; i++)
        {
            mpfr_sub_ui(dd + i, mpfr_mat_entry(Wd, i, i), 2, MPFR_RNDD);
        }
        for (i = 0; i < n; i++)
        {
            mpfr_sub_ui(du + i, mpfr_mat_entry(Wu, i, i), 2, MPFR_RNDU);
        }
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                if (j > i)
                {
                    mpfr_abs(ti, mpfr_mat_entry(Wd, i, j), MPFR_RNDU);
                    mpfr_abs(tj, mpfr_mat_entry(Wu, i, j), MPFR_RNDU);
                    mpfr_max(mpfr_mat_entry(bound, i, j), ti, tj, MPFR_RNDU);
                    mpfr_mul(ti, norm, norm, MPFR_RNDU);
                    mpfr_ui_sub(tj, 1, norm, MPFR_RNDU);
                    mpfr_div(tmp, ti, tj, MPFR_RNDU);
                    mpfr_add(mpfr_mat_entry(bound, i, j),
                             mpfr_mat_entry(bound, i, j), tmp, MPFR_RNDU);
                }
                else if (j < i)
                {
                    mpfr_abs(ti, mpfr_mat_entry(Wd, i, j), MPFR_RNDU);
                    mpfr_abs(tj, mpfr_mat_entry(Wu, i, j), MPFR_RNDU);
                    mpfr_max(mpfr_mat_entry(bound, i, j), ti, tj, MPFR_RNDU);
                }
                else
                {
                    mpfr_abs(ti, dd + i, MPFR_RNDU);
                    mpfr_abs(tj, du + i, MPFR_RNDU);
                    mpfr_max(mpfr_mat_entry(bound, i, j), ti, tj, MPFR_RNDU);
                    mpfr_mul(ti, norm, norm, MPFR_RNDU);
                    mpfr_ui_sub(tj, 1, norm, MPFR_RNDU);
                    mpfr_div(tmp, ti, tj, MPFR_RNDU);
                    mpfr_add(mpfr_mat_entry(bound, i, j),
                             mpfr_mat_entry(bound, i, j), tmp, MPFR_RNDU);
                }
            }
        }
        _mpfr_vec_clear(dd, n);
        _mpfr_vec_clear(du, n);

        mpfr_mat_init(mm, n, n, prec);
        mpfr_mat_init(rm, n, n, prec);
        mpfr_mat_init(mn, n, n, prec);
        mpfr_mat_init(rn, n, n, prec);
        mpfr_mat_init(bound2, n, n, prec);

        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                mpfr_add(tmp, mpfr_mat_entry(Wu, i, j),
                         mpfr_mat_entry(Wd, i, j), MPFR_RNDU);
                mpfr_div_ui(mpfr_mat_entry(mm, j, i), tmp, 2, MPFR_RNDU);
                mpfr_sub(mpfr_mat_entry(rm, j, i), mpfr_mat_entry(mm, j, i),
                         mpfr_mat_entry(Wd, i, j), MPFR_RNDU);
                mpfr_div_ui(mpfr_mat_entry(mn, i, j), tmp, 2, MPFR_RNDU);
                mpfr_sub(mpfr_mat_entry(rn, i, j), mpfr_mat_entry(mn, i, j),
                         mpfr_mat_entry(Wd, i, j), MPFR_RNDU);
            }
        }
        mpfr_mat_mul_classical(Wd, mm, mn, MPFR_RNDD);
        for (i = 0; i < n; i++)
        {
            mpfr_sub_ui(mpfr_mat_entry(Wd, i, i), mpfr_mat_entry(Wd, i, i), 1,
                        MPFR_RNDD);
        }
        mpfr_mat_mul_classical(Wu, mm, mn, MPFR_RNDU);
        for (i = 0; i < n; i++)
        {
            mpfr_sub_ui(mpfr_mat_entry(Wu, i, i), mpfr_mat_entry(Wu, i, i), 1,
                        MPFR_RNDU);
            for (j = 0; j < n; j++)
            {
                mpfr_abs(ti, mpfr_mat_entry(Wd, i, j), MPFR_RNDU);
                mpfr_abs(tj, mpfr_mat_entry(Wu, i, j), MPFR_RNDU);
                mpfr_max(mpfr_mat_entry(Wu, i, j), ti, tj, MPFR_RNDU);
                mpfr_abs(mpfr_mat_entry(mm, i, j), mpfr_mat_entry(mm, i, j),
                         MPFR_RNDU);
                mpfr_abs(mpfr_mat_entry(mn, i, j), mpfr_mat_entry(mn, i, j),
                         MPFR_RNDU);
            }
        }
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                mpfr_add(mpfr_mat_entry(bound2, i, j),
                         mpfr_mat_entry(mn, i, j), mpfr_mat_entry(rn, i, j),
                         MPFR_RNDU);
            }
        }
        mpfr_mat_mul_classical(bound2, rm, bound2, MPFR_RNDU);
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                mpfr_add(mpfr_mat_entry(bound2, i, j),
                         mpfr_mat_entry(bound2, i, j), mpfr_mat_entry(Wu, i,
                                                                      j),
                         MPFR_RNDU);
            }
        }
        mpfr_mat_mul_classical(Wu, mm, rn, MPFR_RNDU);
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                mpfr_add(mpfr_mat_entry(bound2, i, j),
                         mpfr_mat_entry(bound2, i, j), mpfr_mat_entry(Wu, i,
                                                                      j),
                         MPFR_RNDU);
            }
        }

        mpfr_mat_clear(Wu);
        mpfr_mat_clear(Wd);
        mpfr_mat_clear(mm);
        mpfr_mat_clear(mn);
        mpfr_mat_clear(rm);
        mpfr_mat_clear(rn);

        mpfr_mat_init(Wu, m, n, prec);
        mpfr_mat_init(Wd, m, n, prec);
        mpfr_mat_init(mm, n, m, prec);
        mpfr_mat_init(mn, m, n, prec);
        mpfr_mat_init(rm, n, m, prec);
        mpfr_mat_init(rn, m, n, prec);

        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                mpfr_set(mpfr_mat_entry(mm, j, i), mpfr_mat_entry(V, i, j),
                         MPFR_RNDU);
            }
        }
        mpfr_mat_mul_classical(Wd, mm, A, MPFR_RNDD);
        mpfr_mat_mul_classical(Wu, mm, A, MPFR_RNDU);

        mpfr_mat_clear(A);

        mpfr_mat_init(bound3, n, n, prec);

        mpfr_mat_mul_classical(mm, Wd, V, MPFR_RNDD);
        for (i = 0; i < n; i++)
        {
            mpfr_sub_ui(mpfr_mat_entry(mm, i, i), mpfr_mat_entry(mm, i, i), 1,
                        MPFR_RNDD);
        }
        mpfr_mat_mul_classical(rm, Wd, V, MPFR_RNDU);
        for (i = 0; i < n; i++)
        {
            mpfr_sub_ui(mpfr_mat_entry(rm, i, i), mpfr_mat_entry(rm, i, i), 1,
                        MPFR_RNDU);
        }

        mpfr_mat_mul_classical(mn, Wu, V, MPFR_RNDD);
        for (i = 0; i < n; i++)
        {
            mpfr_sub_ui(mpfr_mat_entry(mn, i, i), mpfr_mat_entry(mn, i, i), 1,
                        MPFR_RNDD);
        }
        mpfr_mat_mul_classical(rn, Wu, V, MPFR_RNDU);
        for (i = 0; i < n; i++)
        {
            mpfr_sub_ui(mpfr_mat_entry(rn, i, i), mpfr_mat_entry(rn, i, i), 1,
                        MPFR_RNDU);
        }

        mpfr_mat_clear(Wd);
        mpfr_mat_clear(Wu);
        mpfr_mat_clear(V);

        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                mpfr_abs(ti, mpfr_mat_entry(mm, i, j), MPFR_RNDU);
                mpfr_abs(tj, mpfr_mat_entry(mn, i, j), MPFR_RNDU);
                mpfr_max(mpfr_mat_entry(bound3, i, j), ti, tj, MPFR_RNDU);
                mpfr_abs(tmp, mpfr_mat_entry(rm, i, j), MPFR_RNDU);
                mpfr_max(mpfr_mat_entry(bound3, i, j),
                         mpfr_mat_entry(bound3, i, j), tmp, MPFR_RNDU);
                mpfr_abs(tmp, mpfr_mat_entry(rn, i, j), MPFR_RNDU);
                mpfr_max(mpfr_mat_entry(bound3, i, j),
                         mpfr_mat_entry(bound3, i, j), tmp, MPFR_RNDU);
            }
        }

        mpfr_mat_clear(mm);
        mpfr_mat_clear(mn);
        mpfr_mat_clear(rm);
        mpfr_mat_clear(rn);

        mpfr_mat_init(boundt, n, n, prec);

        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                mpfr_set(mpfr_mat_entry(boundt, j, i),
                         mpfr_mat_entry(bound, i, j), MPFR_RNDU);
                mpfr_set(ti, mpfr_mat_entry(bound2, i, j), MPFR_RNDU);
                mpfr_set(tj, mpfr_mat_entry(bound3, i, j), MPFR_RNDU);
                mpfr_add(mpfr_mat_entry(bound2, i, j), ti, tj, MPFR_RNDU);
            }
        }
        mpfr_mat_mul_classical(bound, bound2, bound, MPFR_RNDU);
        mpfr_mat_mul_classical(bound, boundt, bound, MPFR_RNDU);

        mpfr_mat_clear(bound2);
        mpfr_mat_clear(bound3);
        mpfr_mat_clear(boundt);

        mpfr_set_zero(norm, 1);
        for (i = 0; i < n; i++)
        {
            mpfr_set_zero(s, 1);
            for (j = 0; j < n; j++)
            {
                mpfr_abs(tmp, mpfr_mat_entry(bound, i, j), MPFR_RNDU);
                mpfr_add(s, s, tmp, MPFR_RNDU);
            }
            mpfr_max(norm, norm, s, MPFR_RNDU);
        }
        if (mpfr_cmp_ui(norm, 1) >= 0)
        {
            mpfr_mat_clear(R);
            mpfr_mat_clear(bound);
            mpfr_clears(s, norm, ti, tj, tmp, NULL);
            return 0;
        }

        mpfr_mat_init(absR, n, n, prec);
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                if (j >= i)
                {
                    mpfr_mul(ti, norm, norm, MPFR_RNDU);
                    mpfr_ui_sub(tj, 1, norm, MPFR_RNDU);
                    mpfr_div(tmp, ti, tj, MPFR_RNDU);
                    mpfr_add(mpfr_mat_entry(bound, i, j),
                             mpfr_mat_entry(bound, i, j), tmp, MPFR_RNDU);
                }
                else
                {
                    mpfr_set_zero(mpfr_mat_entry(bound, i, j), 1);
                }
                mpfr_abs(mpfr_mat_entry(absR, i, j), mpfr_mat_entry(R, i, j),
                         MPFR_RNDU);
            }
        }
        mpfr_mat_mul_classical(bound, bound, absR, MPFR_RNDU);

        mpfr_mat_clear(absR);

        mpfr_init2(mpfr_gs_B, prec);
        fmpz_get_mpfr(mpfr_gs_B, gs_B, MPFR_RNDN);
        for (i = 0; i < n - 1; i++)
        {
            mpfr_sub(tmp, mpfr_mat_entry(R, i, i), mpfr_mat_entry(bound, i, i),
                     MPFR_RNDD);
            mpfr_mul_d(ti, tmp, fl->eta, MPFR_RNDD);
            mpfr_sqr(tmp, tmp, MPFR_RNDD);
            if (i >= newd && !mpfr_greaterequal_p(tmp, mpfr_gs_B))
            {
                mpfr_mat_clear(R);
                mpfr_mat_clear(bound);
                mpfr_clears(s, norm, ti, tj, tmp, mpfr_gs_B, NULL);
                return 0;
            }
            for (j = i + 1; j < n; j++)
            {
                mpfr_abs(tmp, mpfr_mat_entry(R, i, j), MPFR_RNDU);
                mpfr_add(tj, tmp, mpfr_mat_entry(bound, i, j), MPFR_RNDU);
                if (i < newd && !mpfr_lessequal_p(tj, ti))
                {
                    mpfr_mat_clear(R);
                    mpfr_mat_clear(bound);
                    mpfr_clears(s, norm, ti, tj, tmp, mpfr_gs_B, NULL);
                    return 0;
                }
            }
            mpfr_add(ti, mpfr_mat_entry(R, i, i), mpfr_mat_entry(bound, i, i),
                     MPFR_RNDU);
            mpfr_sub(tj, mpfr_mat_entry(R, i + 1, i + 1),
                     mpfr_mat_entry(bound, i + 1, i + 1), MPFR_RNDD);
            mpfr_abs(tmp, mpfr_mat_entry(R, i, i + 1), MPFR_RNDD);
            mpfr_sub(norm, tmp, mpfr_mat_entry(bound, i, i + 1), MPFR_RNDD);
            mpfr_div(tmp, norm, ti, MPFR_RNDD);
            mpfr_sqr(norm, tmp, MPFR_RNDD);
            mpfr_sub_d(s, norm, fl->delta, MPFR_RNDD);
            mpfr_neg(s, s, MPFR_RNDD);
            mpfr_sqrt(tmp, s, MPFR_RNDU);
            mpfr_mul(s, tmp, ti, MPFR_RNDU);
            if (i < newd && !mpfr_lessequal_p(s, tj))
            {
                mpfr_mat_clear(R);
                mpfr_mat_clear(bound);
                mpfr_clears(s, norm, ti, tj, tmp, mpfr_gs_B, NULL);
                return 0;
            }
        }
        mpfr_sub(tmp, mpfr_mat_entry(R, i, i), mpfr_mat_entry(bound, i, i),
                 MPFR_RNDD);
        mpfr_sqr(tmp, tmp, MPFR_RNDD);
        if (i >= newd && !mpfr_greaterequal_p(tmp, mpfr_gs_B))
        {
            mpfr_mat_clear(R);
            mpfr_mat_clear(bound);
            mpfr_clears(s, norm, ti, tj, tmp, mpfr_gs_B, NULL);
            return 0;
        }

        mpfr_mat_clear(R);
        mpfr_mat_clear(bound);
        mpfr_clears(s, norm, ti, tj, tmp, mpfr_gs_B, NULL);
    }

    FLINT_ASSERT((fl->rt == Z_BASIS
          ? fmpz_mat_is_reduced_with_removal(B, fl->delta, fl->eta, gs_B, newd)
          : fmpz_mat_is_reduced_gram_with_removal(B, fl->delta, fl->eta, gs_B, newd)));

    return 1;
}
