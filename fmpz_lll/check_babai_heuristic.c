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

    Copyright (C) 2009, 2010 William Hart
    Copyright (C) 2009, 2010 Andy Novocin
    Copyright (C) 2014 Abhinav Baid

******************************************************************************/

#include "fmpz_lll.h"

int
fmpz_lll_check_babai_heuristic(int kappa, fmpz_mat_t B, mpf_mat_t mu,
                               mpf_mat_t r, mpf * s, mpf_mat_t appB,
                               fmpz_gram_t A, int a, int zeros, int kappamax,
                               int n, mpf_t tmp, mpf_t rtmp, mp_bitcnt_t prec,
                               const fmpz_lll_t fl)
{
    if (fl->rt == Z_BASIS)
    {
        if (fl->gt == APPROX)
        {
            int i, j, k, test, aa;
            slong xx, exponent;
            fmpz_t ztmp;
            double halfplus, onedothalfplus;
            ulong loops;

            fmpz_init(ztmp);

            aa = (a > zeros) ? a : zeros + 1;

            halfplus = (4 * fl->eta + 0.5) / 5;
            onedothalfplus = 1.0 + halfplus;

            loops = 0;

            do
            {
                test = 0;

                loops++;
                if (loops > 20)
                {
                    return -1;
                }

                /* ************************************** */
                /* Step2: compute the GSO for stage kappa */
                /* ************************************** */

                for (j = aa; j < kappa; j++)
                {
                    if (mpf_cmp_d(mpf_mat_entry(A->appSP2, kappa, j), DBL_MIN)
                        == 0)
                    {
                        if (!
                            (_mpf_vec_dot2
                             (mpf_mat_entry(A->appSP2, kappa, j),
                              appB->rows[kappa], appB->rows[j], n, prec)))
                        {
/* In this case a heuristic told us that some cancelation probably happened so we just compute the scalar product at full precision */
                            _fmpz_vec_dot(ztmp, B->rows[kappa], B->rows[j], n);
                            fmpz_get_mpf(mpf_mat_entry(A->appSP2, kappa, j),
                                         ztmp);
                        }
                    }

                    if (j > zeros + 2)
                    {
                        mpf_mul(tmp, mpf_mat_entry(mu, j, zeros + 1),
                                mpf_mat_entry(r, kappa, zeros + 1));
                        mpf_sub(rtmp, mpf_mat_entry(A->appSP2, kappa, j), tmp);

                        for (k = zeros + 2; k < j - 1; k++)
                        {
                            mpf_mul(tmp, mpf_mat_entry(mu, j, k),
                                    mpf_mat_entry(r, kappa, k));
                            mpf_sub(rtmp, rtmp, tmp);
                        }

                        mpf_mul(tmp, mpf_mat_entry(mu, j, j - 1),
                                mpf_mat_entry(r, kappa, j - 1));
                        mpf_sub(mpf_mat_entry(r, kappa, j), rtmp, tmp);
                    }
                    else if (j == zeros + 2)
                    {
                        mpf_mul(tmp, mpf_mat_entry(mu, j, zeros + 1),
                                mpf_mat_entry(r, kappa, zeros + 1));
                        mpf_sub(mpf_mat_entry(r, kappa, j),
                                mpf_mat_entry(A->appSP2, kappa, j), tmp);
                    }
                    else
                        mpf_set(mpf_mat_entry(r, kappa, j),
                                mpf_mat_entry(A->appSP2, kappa, j));

                    mpf_div(mpf_mat_entry(mu, kappa, j),
                            mpf_mat_entry(r, kappa, j), mpf_mat_entry(r, j,
                                                                      j));
                }

                /* **************************** */
                /* Step3--5: compute the X_j's  */
                /* **************************** */

                for (j = kappa - 1; j > zeros; j--)
                {
                    /* test of the relaxed size-reduction condition */
                    mpf_abs(tmp, mpf_mat_entry(mu, kappa, j));

                    if (mpf_cmp_d(tmp, halfplus) > 0)
                    {
                        test = 1;

                        /* we consider separately the cases X = +-1 */
                        if (mpf_cmp_d(tmp, onedothalfplus) <= 0)
                        {
                            int sgn = mpf_sgn(mpf_mat_entry(mu, kappa, j));
                            if (sgn >= 0)   /* in this case, X is 1 */
                            {
                                for (k = zeros + 1; k < j; k++)
                                {
                                    mpf_sub(mpf_mat_entry(mu, kappa, k),
                                            mpf_mat_entry(mu, kappa, k),
                                            mpf_mat_entry(mu, j, k));
                                }
                                _fmpz_vec_sub(B->rows[kappa], B->rows[kappa],
                                              B->rows[j], n);
                            }
                            else    /* otherwise X is -1 */
                            {
                                for (k = zeros + 1; k < j; k++)
                                {
                                    mpf_add(mpf_mat_entry(mu, kappa, k),
                                            mpf_mat_entry(mu, kappa, k),
                                            mpf_mat_entry(mu, j, k));
                                }
                                _fmpz_vec_add(B->rows[kappa], B->rows[kappa],
                                              B->rows[j], n);
                            }
                        }
                        else    /* we must have |X| >= 2 */
                        {
                            mpf_set(tmp, mpf_mat_entry(mu, kappa, j));
                            mpf_set_d(rtmp, 0.5);
                            if (mpf_cmp_ui(tmp, 0) < 0)
                            {
                                mpf_sub(tmp, tmp, rtmp);
                                mpf_ceil(tmp, tmp);
                            }
                            else
                            {
                                mpf_add(tmp, tmp, rtmp);
                                mpf_floor(tmp, tmp);
                            }
                            for (k = zeros + 1; k < j; k++)
                            {
                                mpf_mul(rtmp, tmp, mpf_mat_entry(mu, j, k));
                                mpf_sub(mpf_mat_entry(mu, kappa, k),
                                        mpf_mat_entry(mu, kappa, k), rtmp);
                            }

                            mpf_get_d_2exp(&exponent, tmp);
                            if (exponent < CPU_SIZE_1 - 2)
                            {
                                /* X is stored in an slong */
                                xx = mpf_get_si(tmp);
                                _fmpz_vec_scalar_submul_si(B->rows[kappa],
                                                           B->rows[j], n, xx);
                            }
                            else
                            {
                                fmpz_set_mpf(ztmp, tmp);
                                _fmpz_vec_scalar_submul_fmpz(B->rows[kappa],
                                                             B->rows[j], n,
                                                             ztmp);
                            }
                        }
                    }
                }

                if (test)       /* Anything happened? */
                {
                    _fmpz_vec_get_mpf_vec(appB->rows[kappa],
                                          B->rows[kappa], n);
                    aa = zeros + 1;

                    for (i = zeros + 1; i <= kappa; i++)
                        mpf_set_d(mpf_mat_entry(A->appSP2, kappa, i), DBL_MIN);

                    for (i = kappa + 1; i <= kappamax; i++)
                        mpf_set_d(mpf_mat_entry(A->appSP2, i, kappa), DBL_MIN);
                }
            } while (test);

            if (mpf_cmp_d(mpf_mat_entry(A->appSP2, kappa, kappa), DBL_MIN) ==
                0)
            {
                _mpf_vec_norm2(mpf_mat_entry(A->appSP2, kappa, kappa),
                               appB->rows[kappa], n, prec);
            }

            mpf_set(s + zeros + 1, mpf_mat_entry(A->appSP2, kappa, kappa));

            for (k = zeros + 1; k < kappa; k++)
            {
                mpf_mul(tmp, mpf_mat_entry(mu, kappa, k),
                        mpf_mat_entry(r, kappa, k));
                mpf_sub(s + k + 1, s + k, tmp);
            }
            fmpz_clear(ztmp);
        }
        else
        {
            int i, j, k, test, aa;
            slong xx, exponent;
            fmpz_t ztmp;
            double halfplus, onedothalfplus;
            ulong loops;

            fmpz_init(ztmp);

            aa = (a > zeros) ? a : zeros + 1;

            halfplus = (4 * fl->eta + 0.5) / 5;
            onedothalfplus = 1.0 + halfplus;

            loops = 0;

            do
            {
                test = 0;

                loops++;
                if (loops > 20)
                {
                    return -1;
                }

                /* ************************************** */
                /* Step2: compute the GSO for stage kappa */
                /* ************************************** */

                for (j = aa; j < kappa; j++)
                {
                    if (fmpz_cmp_si
                        (fmpz_mat_entry(A->exactSP, kappa, j), WORD_MIN) == 0)
                    {
                        _fmpz_vec_dot(fmpz_mat_entry(A->exactSP, kappa, j),
                                      B->rows[kappa], B->rows[j], n);
                    }

                    if (j > zeros + 2)
                    {
                        mpf_mul(tmp, mpf_mat_entry(mu, j, zeros + 1),
                                mpf_mat_entry(r, kappa, zeros + 1));
                        fmpz_get_mpf(rtmp,
                                     fmpz_mat_entry(A->exactSP, kappa, j));
                        mpf_sub(rtmp, rtmp, tmp);

                        for (k = zeros + 2; k < j - 1; k++)
                        {
                            mpf_mul(tmp, mpf_mat_entry(mu, j, k),
                                    mpf_mat_entry(r, kappa, k));
                            mpf_sub(rtmp, rtmp, tmp);
                        }

                        mpf_mul(tmp, mpf_mat_entry(mu, j, j - 1),
                                mpf_mat_entry(r, kappa, j - 1));
                        mpf_sub(mpf_mat_entry(r, kappa, j), rtmp, tmp);
                    }
                    else if (j == zeros + 2)
                    {
                        mpf_mul(tmp, mpf_mat_entry(mu, j, zeros + 1),
                                mpf_mat_entry(r, kappa, zeros + 1));
                        fmpz_get_mpf(mpf_mat_entry(r, kappa, j),
                                     fmpz_mat_entry(A->exactSP, kappa, j));
                        mpf_sub(mpf_mat_entry(r, kappa, j),
                                mpf_mat_entry(r, kappa, j), tmp);
                    }
                    else
                        fmpz_get_mpf(mpf_mat_entry(r, kappa, j),
                                     fmpz_mat_entry(A->exactSP, kappa, j));

                    mpf_div(mpf_mat_entry(mu, kappa, j),
                            mpf_mat_entry(r, kappa, j), mpf_mat_entry(r, j,
                                                                      j));
                }

                /* **************************** */
                /* Step3--5: compute the X_j's  */
                /* **************************** */

                for (j = kappa - 1; j > zeros; j--)
                {
                    /* test of the relaxed size-reduction condition */
                    mpf_abs(tmp, mpf_mat_entry(mu, kappa, j));

                    if (mpf_cmp_d(tmp, halfplus) > 0)
                    {
                        test = 1;

                        /* we consider separately the cases X = +-1 */
                        if (mpf_cmp_d(tmp, onedothalfplus) <= 0)
                        {
                            int sgn = mpf_sgn(mpf_mat_entry(mu, kappa, j));
                            if (sgn >= 0)   /* in this case, X is 1 */
                            {
                                for (k = zeros + 1; k < j; k++)
                                {
                                    mpf_sub(mpf_mat_entry(mu, kappa, k),
                                            mpf_mat_entry(mu, kappa, k),
                                            mpf_mat_entry(mu, j, k));
                                }
                                _fmpz_vec_sub(B->rows[kappa], B->rows[kappa],
                                              B->rows[j], n);
                            }
                            else    /* otherwise X is -1 */
                            {
                                for (k = zeros + 1; k < j; k++)
                                {
                                    mpf_add(mpf_mat_entry(mu, kappa, k),
                                            mpf_mat_entry(mu, kappa, k),
                                            mpf_mat_entry(mu, j, k));
                                }
                                _fmpz_vec_add(B->rows[kappa], B->rows[kappa],
                                              B->rows[j], n);
                            }
                        }
                        else    /* we must have |X| >= 2 */
                        {
                            mpf_set(tmp, mpf_mat_entry(mu, kappa, j));
                            mpf_set_d(rtmp, 0.5);
                            if (mpf_cmp_ui(tmp, 0) < 0)
                            {
                                mpf_sub(tmp, tmp, rtmp);
                                mpf_ceil(tmp, tmp);
                            }
                            else
                            {
                                mpf_add(tmp, tmp, rtmp);
                                mpf_floor(tmp, tmp);
                            }
                            for (k = zeros + 1; k < j; k++)
                            {
                                mpf_mul(rtmp, tmp, mpf_mat_entry(mu, j, k));
                                mpf_sub(mpf_mat_entry(mu, kappa, k),
                                        mpf_mat_entry(mu, kappa, k), rtmp);
                            }

                            mpf_get_d_2exp(&exponent, tmp);
                            if (exponent < CPU_SIZE_1 - 2)
                            {
                                /* X is stored in an slong */
                                xx = mpf_get_si(tmp);
                                _fmpz_vec_scalar_submul_si(B->rows[kappa],
                                                           B->rows[j], n, xx);
                            }
                            else
                            {
                                fmpz_set_mpf(ztmp, tmp);
                                _fmpz_vec_scalar_submul_fmpz(B->rows[kappa],
                                                             B->rows[j], n,
                                                             ztmp);
                            }
                        }
                    }
                }

                if (test)       /* Anything happened? */
                {
                    _fmpz_vec_get_mpf_vec(appB->rows[kappa],
                                          B->rows[kappa], n);
                    aa = zeros + 1;

                    for (i = zeros + 1; i <= kappa; i++)
                    {
                        fmpz_set_si(fmpz_mat_entry(A->exactSP, kappa, i),
                                    WORD_MIN);
                    }

                    for (i = kappa + 1; i <= kappamax; i++)
                    {
                        fmpz_set_si(fmpz_mat_entry(A->exactSP, i, kappa),
                                    WORD_MIN);
                    }
                }
            } while (test);

            if (fmpz_cmp_si(fmpz_mat_entry(A->exactSP, kappa, kappa), WORD_MIN)
                == 0)
            {
                _fmpz_vec_dot(fmpz_mat_entry(A->exactSP, kappa, kappa),
                              B->rows[kappa], B->rows[kappa], n);
            }

            fmpz_get_mpf(s + zeros + 1,
                         fmpz_mat_entry(A->exactSP, kappa, kappa));

            for (k = zeros + 1; k < kappa; k++)
            {
                mpf_mul(tmp, mpf_mat_entry(mu, kappa, k),
                        mpf_mat_entry(r, kappa, k));
                mpf_sub(s + k + 1, s + k, tmp);
            }
            fmpz_clear(ztmp);
        }
    }
    return 0;
}
