/*
    Copyright (C) 2009, 2010 William Hart
    Copyright (C) 2009, 2010 Andy Novocin
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_lll.h"
#ifdef GM
#undef GM
#endif
#define GM ((fl->rt == Z_BASIS) ? A->exactSP : B)

int
fmpz_lll_check_babai_heuristic(int kappa, fmpz_mat_t B, fmpz_mat_t U,
                               mpf_mat_t mu, mpf_mat_t r, mpf * s,
                               mpf_mat_t appB, fmpz_gram_t A, int a, int zeros,
                               int kappamax, int n, mpf_t tmp, mpf_t rtmp,
                               flint_bitcnt_t prec, const fmpz_lll_t fl)
{
    if (fl->rt == Z_BASIS && fl->gt == APPROX)
    {
        int i, j, k, test, aa;
        slong xx, exponent, max_expo = WORD_MAX;
        fmpz_t ztmp;
        double halfplus, onedothalfplus;
        ulong loops;

        fmpz_init(ztmp);

        aa = (a > zeros) ? a : zeros + 1;

        halfplus = (fl->eta + 0.5) / 2;
        onedothalfplus = 1.0 + halfplus;

        loops = 0;

        do
        {
            test = 0;

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
                        fmpz_get_mpf(mpf_mat_entry(A->appSP2, kappa, j), ztmp);
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
                        mpf_mat_entry(r, kappa, j), mpf_mat_entry(r, j, j));
            }

            if (loops >= 20)
            {
                slong new_max_expo = WORD_MIN;
                for (j = 0; j < kappa; j++)
                {
                    slong expo2;
                    flint_mpf_get_d_2exp(&expo2, mpf_mat_entry(mu, kappa, j));
                    new_max_expo = FLINT_MAX(new_max_expo, expo2);
                }
                if (new_max_expo > max_expo - SIZE_RED_FAILURE_THRESH)
                {
                    fmpz_clear(ztmp);
                    return -1;
                }
                max_expo = new_max_expo;
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
                            if (U != NULL)
                            {
                                _fmpz_vec_sub(U->rows[kappa],
                                              U->rows[kappa], U->rows[j],
                                              U->c);
                            }
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
                            if (U != NULL)
                            {
                                _fmpz_vec_add(U->rows[kappa],
                                              U->rows[kappa], U->rows[j],
                                              U->c);
                            }
                        }
                    }
                    else        /* we must have |X| >= 2 */
                    {
                        mpf_set(tmp, mpf_mat_entry(mu, kappa, j));
                        mpf_set_d(rtmp, 0.5);
                        if (flint_mpf_cmp_ui(tmp, 0) < 0)
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

                        flint_mpf_get_d_2exp(&exponent, tmp);
                        if (exponent < CPU_SIZE_1 - 2)
                        {
                            /* X is stored in an slong */
                            xx = flint_mpf_get_si(tmp);
                            _fmpz_vec_scalar_submul_si(B->rows[kappa],
                                                       B->rows[j], n, xx);
                            if (U != NULL)
                            {
                                _fmpz_vec_scalar_submul_si(U->rows[kappa],
                                                           U->rows[j],
                                                           U->c, xx);
                            }
                        }
                        else
                        {
                            fmpz_set_mpf(ztmp, tmp);
                            _fmpz_vec_scalar_submul_fmpz(B->rows[kappa],
                                                         B->rows[j], n, ztmp);
                            if (U != NULL)
                            {
                                _fmpz_vec_scalar_submul_fmpz(U->rows
                                                             [kappa],
                                                             U->rows[j],
                                                             U->c, ztmp);
                            }
                        }
                    }
                }
            }

            if (test)           /* Anything happened? */
            {
                _fmpz_vec_get_mpf_vec(appB->rows[kappa], B->rows[kappa], n);
                aa = zeros + 1;

                for (i = zeros + 1; i <= kappa; i++)
                    mpf_set_d(mpf_mat_entry(A->appSP2, kappa, i), DBL_MIN);

                for (i = kappa + 1; i <= kappamax; i++)
                    mpf_set_d(mpf_mat_entry(A->appSP2, i, kappa), DBL_MIN);
            }
            loops++;
        } while (test);

        if (mpf_cmp_d(mpf_mat_entry(A->appSP2, kappa, kappa), DBL_MIN) == 0)
        {
            _mpf_vec_norm2(mpf_mat_entry(A->appSP2, kappa, kappa),
                           appB->rows[kappa], n, prec);
        }

        mpf_set(s + zeros + 1, mpf_mat_entry(A->appSP2, kappa, kappa));

        for (k = zeros + 1; k < kappa - 1; k++)
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
        slong xx, exponent, max_expo = WORD_MAX;
        fmpz_t ztmp;
        double halfplus, onedothalfplus;
        ulong loops;

        fmpz_init(ztmp);

        aa = (a > zeros) ? a : zeros + 1;

        halfplus = (fl->eta + 0.5) / 2;
        onedothalfplus = 1.0 + halfplus;

        loops = 0;

        do
        {
            fmpz *x;

            test = 0;

            /* ************************************** */
            /* Step2: compute the GSO for stage kappa */
            /* ************************************** */

            for (j = aa; j < kappa; j++)
            {
                if (j > zeros + 2)
                {
                    mpf_mul(tmp, mpf_mat_entry(mu, j, zeros + 1),
                            mpf_mat_entry(r, kappa, zeros + 1));
                    fmpz_get_mpf(rtmp, fmpz_mat_entry(GM, kappa, j));
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
                                 fmpz_mat_entry(GM, kappa, j));
                    mpf_sub(mpf_mat_entry(r, kappa, j),
                            mpf_mat_entry(r, kappa, j), tmp);
                }
                else
                    fmpz_get_mpf(mpf_mat_entry(r, kappa, j),
                                 fmpz_mat_entry(GM, kappa, j));

                mpf_div(mpf_mat_entry(mu, kappa, j),
                        mpf_mat_entry(r, kappa, j), mpf_mat_entry(r, j, j));
            }

            if (loops >= 20)
            {
                slong new_max_expo = WORD_MIN;
                for (j = 0; j < kappa; j++)
                {
                    slong expo2;
                    flint_mpf_get_d_2exp(&expo2, mpf_mat_entry(mu, kappa, j));
                    new_max_expo = FLINT_MAX(new_max_expo, expo2);
                }
                if (new_max_expo > max_expo - SIZE_RED_FAILURE_THRESH)
                {
                    fmpz_clear(ztmp);
                    return -1;
                }
                max_expo = new_max_expo;
            }

            /* **************************** */
            /* Step3--5: compute the X_j's  */
            /* **************************** */

            x = _fmpz_vec_init(kappa - 1 - zeros);
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
                            fmpz_set_ui(x + j, 1);
                            for (k = zeros + 1; k < j; k++)
                            {
                                mpf_sub(mpf_mat_entry(mu, kappa, k),
                                        mpf_mat_entry(mu, kappa, k),
                                        mpf_mat_entry(mu, j, k));
                            }
                            if (fl->rt == Z_BASIS && B != NULL)
                            {
                                _fmpz_vec_sub(B->rows[kappa], B->rows[kappa],
                                              B->rows[j], n);
                            }
                            if (U != NULL)
                            {
                                _fmpz_vec_sub(U->rows[kappa],
                                              U->rows[kappa], U->rows[j],
                                              U->c);
                            }
                        }
                        else    /* otherwise X is -1 */
                        {
                            fmpz_set_si(x + j, -WORD(1));
                            for (k = zeros + 1; k < j; k++)
                            {
                                mpf_add(mpf_mat_entry(mu, kappa, k),
                                        mpf_mat_entry(mu, kappa, k),
                                        mpf_mat_entry(mu, j, k));
                            }
                            if (fl->rt == Z_BASIS && B != NULL)
                            {
                                _fmpz_vec_add(B->rows[kappa],
                                              B->rows[kappa], B->rows[j], n);
                            }
                            if (U != NULL)
                            {
                                _fmpz_vec_add(U->rows[kappa],
                                              U->rows[kappa], U->rows[j],
                                              U->c);
                            }
                        }
                    }
                    else        /* we must have |X| >= 2 */
                    {
                        mpf_set(tmp, mpf_mat_entry(mu, kappa, j));
                        mpf_set_d(rtmp, 0.5);
                        if (flint_mpf_cmp_ui(tmp, 0) < 0)
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

                        flint_mpf_get_d_2exp(&exponent, tmp);
                        if (exponent < CPU_SIZE_1 - 2)
                        {
                            /* X is stored in an slong */
                            xx = flint_mpf_get_si(tmp);
                            fmpz_set_si(x + j, xx);
                            if (fl->rt == Z_BASIS && B != NULL)
                            {
                                _fmpz_vec_scalar_submul_si(B->rows[kappa],
                                                           B->rows[j], n, xx);
                            }
                            if (U != NULL)
                            {
                                _fmpz_vec_scalar_submul_si(U->rows[kappa],
                                                           U->rows[j],
                                                           U->c, xx);
                            }
                        }
                        else
                        {
                            fmpz_set_mpf(x + j, tmp);
                            if (fl->rt == Z_BASIS && B != NULL)
                            {
                                _fmpz_vec_scalar_submul_fmpz(B->rows[kappa],
                                                             B->rows[j], n,
                                                             x + j);
                            }
                            if (U != NULL)
                            {
                                _fmpz_vec_scalar_submul_fmpz(U->rows
                                                             [kappa],
                                                             U->rows[j],
                                                             U->c, x + j);
                            }
                        }
                    }
                }
            }

            if (test)           /* Anything happened? */
            {
                aa = zeros + 1;

                for (j = zeros + 1; j < kappa; j++)
                {
                    fmpz_pow_ui(ztmp, x + j, 2);
                    fmpz_addmul(fmpz_mat_entry(GM, kappa, kappa),
                                ztmp, fmpz_mat_entry(GM, j, j));

                    fmpz_mul(ztmp, x + j, fmpz_mat_entry(GM, kappa, j));
                    fmpz_mul_2exp(ztmp, ztmp, 1);
                    fmpz_sub(fmpz_mat_entry(GM, kappa, kappa),
                             fmpz_mat_entry(GM, kappa, kappa), ztmp);

                    for (i = zeros + 1; i < j; i++)
                    {
                        fmpz_mul(ztmp, x + i, x + j);
                        fmpz_mul(ztmp, ztmp, fmpz_mat_entry(GM, j, i));
                        fmpz_mul_2exp(ztmp, ztmp, 1);
                        fmpz_add(fmpz_mat_entry(GM, kappa, kappa),
                                 fmpz_mat_entry(GM, kappa, kappa), ztmp);
                    }
                }

                for (i = zeros + 1; i < kappa; i++)
                {
                    for (j = zeros + 1; j <= i; j++)
                        fmpz_submul(fmpz_mat_entry(GM, kappa, i),
                                    x + j, fmpz_mat_entry(GM, i, j));
                    for (j = i + 1; j < kappa; j++)
                        fmpz_submul(fmpz_mat_entry(GM, kappa, i),
                                    x + j, fmpz_mat_entry(GM, j, i));
                }

                for (i = kappa + 1; i < GM->r; i++)
                {
                    for (j = zeros + 1; j < kappa; j++)
                        fmpz_submul(fmpz_mat_entry(GM, i, kappa),
                                    x + j, fmpz_mat_entry(GM, i, j));
                }
            }

            _fmpz_vec_clear(x, kappa - 1 - zeros);
            loops++;
        } while (test);

        fmpz_get_mpf(s + zeros + 1, fmpz_mat_entry(GM, kappa, kappa));

        for (k = zeros + 1; k < kappa - 1; k++)
        {
            mpf_mul(tmp, mpf_mat_entry(mu, kappa, k),
                    mpf_mat_entry(r, kappa, k));
            mpf_sub(s + k + 1, s + k, tmp);
        }

        fmpz_clear(ztmp);
    }
    return 0;
}

#undef GM
