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
fmpz_lll_mpf2(fmpz_mat_t B, mp_bitcnt_t prec, const fmpz_lll_t fl)
{
    if (fl->rt == Z_BASIS)
    {
        if (fl->gt == APPROX)
        {
            int kappa, kappa2, d, n, i, j, zeros, kappamax;
            mpf_mat_t mu, r, appB;
            fmpz_gram_t A;
            mpf *s, *mutmp, *appBtmp, *appSPtmp;
            mpf_t ctt, tmp, rtmp;
            int *alpha;
            fmpz *Btmp;

            n = B->c;
            d = B->r;

            mpf_init_set_d(ctt, (4 * fl->delta + 1) / 5);

            alpha = (int *) flint_malloc(d * sizeof(int));

            mpf_init2(tmp, prec);
            mpf_init2(rtmp, prec);

            mpf_mat_init(mu, d, d, prec);
            mpf_mat_init(r, d, d, prec);
            mpf_mat_init(appB, d, n, prec);
            mpf_mat_init(A->appSP2, d, d, prec);

            s = _mpf_vec_init(d, prec);
            appSPtmp = _mpf_vec_init(d, prec);

            for (i = 0; i < d; i++)
            {
                for (j = 0; j < d; j++)
                {
                    mpf_set_d(mpf_mat_entry(A->appSP2, i, j), DBL_MIN);
                }
            }

            /* ************************** */
            /* Step1: Initialization Step */
            /* ************************** */

            for (i = 0; i < d; i++)
                _fmpz_vec_get_mpf_vec(appB->rows[i], B->rows[i], n);

            /* ********************************* */
            /* Step2: Initializing the main loop */
            /* ********************************* */

            kappamax = 0;
            i = 0;

            do
            {
                _mpf_vec_norm2(mpf_mat_entry(A->appSP2, i, i), appB->rows[i],
                               n, prec);
            } while ((mpf_sgn(mpf_mat_entry(A->appSP2, i, i)) == 0.0)
                     && (++i < d));

            zeros = i - 1;      /* all vectors B[i] with i <= zeros are zero vectors */
            kappa = i + 1;
            kappamax = kappa;

            if (zeros < d - 1)
            {
                mpf_set(mpf_mat_entry(r, i, i),
                        mpf_mat_entry(A->appSP2, i, i));
            }

            for (i = zeros + 1; i < d; i++)
                alpha[i] = 0;

            while (kappa < d)
            {
                int babai_ok = 0;

                if (kappa > kappamax)
                    kappamax = kappa;

                /* ********************************** */
                /* Step3: Call to the Babai algorithm */
                /* ********************************** */
                babai_ok =
                    fmpz_lll_check_babai_heuristic(kappa, B, mu, r, s, appB, A,
                                                   alpha[kappa], zeros,
                                                   kappamax, n, tmp, rtmp,
                                                   prec, fl);

                if (babai_ok == -1)
                {
                    flint_free(alpha);
                    mpf_clears(ctt, tmp, rtmp, '\0');
                    mpf_mat_clear(mu);
                    mpf_mat_clear(r);
                    mpf_mat_clear(appB);
                    mpf_mat_clear(A->appSP2);
                    _mpf_vec_clear(s, d);
                    _mpf_vec_clear(appSPtmp, d);
                    /* Need to switch to mpf / arb */
                    return -1;
                }

                /* ************************************ */
                /* Step4: Success of Lovasz's condition */
                /* ************************************ */

                mpf_mul(tmp, mpf_mat_entry(r, kappa - 1, kappa - 1), ctt);

                if (mpf_cmp(tmp, s + kappa - 1) <= 0)
                {
                    alpha[kappa] = kappa;
                    mpf_mul(tmp, mpf_mat_entry(mu, kappa, kappa - 1),
                            mpf_mat_entry(r, kappa, kappa - 1));
                    mpf_sub(mpf_mat_entry(r, kappa, kappa), s + kappa - 1,
                            tmp);
                    kappa++;
                }
                else
                {

                    /* ******************************************* */
                    /* Step5: Find the right insertion index kappa */
                    /* kappa2 remains the initial kappa            */
                    /* ******************************************* */

                    kappa2 = kappa;
                    do
                    {
                        kappa--;
                        if (kappa > zeros + 1)
                        {
                            mpf_mul(tmp,
                                    mpf_mat_entry(r, kappa - 1, kappa - 1),
                                    ctt);
                        }
                    } while ((kappa >= zeros + 2)
                             && (mpf_cmp(s + kappa - 1, tmp) <= 0));

                    for (i = kappa; i < kappa2; i++)
                        if (kappa <= alpha[i])
                            alpha[i] = kappa;

                    for (i = kappa2; i > kappa; i--)
                        alpha[i] = alpha[i - 1];

                    for (i = kappa2 + 1; i <= kappamax; i++)
                        if (kappa < alpha[i])
                            alpha[i] = kappa;

                    alpha[kappa] = kappa;

                    /* ****************************** */
                    /* Step6: Update the mu's and r's */
                    /* ****************************** */

                    mutmp = mu->rows[kappa2];
                    for (i = kappa2; i > kappa; i--)
                        mu->rows[i] = mu->rows[i - 1];
                    mu->rows[kappa] = mutmp;

                    mutmp = r->rows[kappa2];
                    for (i = kappa2; i > kappa; i--)
                        r->rows[i] = r->rows[i - 1];
                    r->rows[kappa] = mutmp;

                    mpf_set(mpf_mat_entry(r, kappa, kappa), s + kappa);

                    /* ************************ */
                    /* Step7: Update B and appB */
                    /* ************************ */

                    Btmp = B->rows[kappa2];
                    for (i = kappa2; i > kappa; i--)
                        B->rows[i] = B->rows[i - 1];
                    B->rows[kappa] = Btmp;

                    appBtmp = appB->rows[kappa2];
                    for (i = kappa2; i > kappa; i--)
                        appB->rows[i] = appB->rows[i - 1];
                    appB->rows[kappa] = appBtmp;

                    /* *************************** */
                    /* Step8: Update appSP: tricky */
                    /* *************************** */

                    for (i = 0; i <= kappa2; i++)
                        mpf_set(appSPtmp + i,
                                mpf_mat_entry(A->appSP2, kappa2, i));

                    for (i = kappa2 + 1; i <= kappamax; i++)
                        mpf_set(appSPtmp + i,
                                mpf_mat_entry(A->appSP2, i, kappa2));

                    for (i = kappa2; i > kappa; i--)
                    {
                        for (j = 0; j < kappa; j++)
                            mpf_set(mpf_mat_entry(A->appSP2, i, j),
                                    mpf_mat_entry(A->appSP2, i - 1, j));
                        mpf_set(mpf_mat_entry(A->appSP2, i, kappa),
                                appSPtmp + i - 1);

                        for (j = kappa + 1; j <= i; j++)
                            mpf_set(mpf_mat_entry(A->appSP2, i, j),
                                    mpf_mat_entry(A->appSP2, i - 1, j - 1));

                        for (j = kappa2 + 1; j <= kappamax; j++)
                            mpf_set(mpf_mat_entry(A->appSP2, j, i),
                                    mpf_mat_entry(A->appSP2, j, i - 1));
                    }

                    for (i = 0; i < kappa; i++)
                        mpf_set(mpf_mat_entry(A->appSP2, kappa, i),
                                appSPtmp + i);
                    mpf_set(mpf_mat_entry(A->appSP2, kappa, kappa),
                            appSPtmp + kappa2);

                    for (i = kappa2 + 1; i <= kappamax; i++)
                        mpf_set(mpf_mat_entry(A->appSP2, i, kappa),
                                appSPtmp + i);

                    if (mpf_sgn(mpf_mat_entry(r, kappa, kappa)) <= 0.0)
                    {
                        zeros++;
                        kappa++;
                        _mpf_vec_norm2(mpf_mat_entry(A->appSP2, kappa, kappa),
                                       appB->rows[kappa], n, prec);
                        mpf_set(mpf_mat_entry(r, kappa, kappa),
                                mpf_mat_entry(A->appSP2, kappa, kappa));
                    }

                    kappa++;
                }
            }

            flint_free(alpha);
            mpf_clears(ctt, tmp, rtmp, '\0');
            mpf_mat_clear(mu);
            mpf_mat_clear(r);
            mpf_mat_clear(appB);
            mpf_mat_clear(A->appSP2);
            _mpf_vec_clear(s, d);
            _mpf_vec_clear(appSPtmp, d);
        }
        else
        {
            int kappa, kappa2, d, n, i, j, zeros, kappamax;
            mpf_mat_t mu, r, appB;
            fmpz_gram_t A;
            mpf *s, *mutmp, *appBtmp;
            fmpz *exactSPtmp;
            mpf_t ctt, tmp, rtmp;
            int *alpha;
            fmpz *Btmp;

            n = B->c;
            d = B->r;

            mpf_init_set_d(ctt, (4 * fl->delta + 1) / 5);

            alpha = (int *) flint_malloc(d * sizeof(int));

            mpf_init2(tmp, prec);
            mpf_init2(rtmp, prec);

            mpf_mat_init(mu, d, d, prec);
            mpf_mat_init(r, d, d, prec);
            mpf_mat_init(appB, d, n, prec);
            fmpz_mat_init(A->exactSP, d, d);

            s = _mpf_vec_init(d, prec);
            exactSPtmp = _fmpz_vec_init(d);

            for (i = 0; i < d; i++)
            {
                for (j = 0; j < d; j++)
                {
                    fmpz_set_si(fmpz_mat_entry(A->exactSP, i, j), WORD_MIN);
                }
            }

            /* ************************** */
            /* Step1: Initialization Step */
            /* ************************** */

            for (i = 0; i < d; i++)
                _fmpz_vec_get_mpf_vec(appB->rows[i], B->rows[i], n);

            /* ********************************* */
            /* Step2: Initializing the main loop */
            /* ********************************* */

            kappamax = 0;
            i = 0;

            do
                _fmpz_vec_dot(fmpz_mat_entry(A->exactSP, i, i), B->rows[i],
                              B->rows[i], n);
            while ((fmpz_cmp_ui(fmpz_mat_entry(A->exactSP, i, i), 0) <= 0)
                   && (++i < d));

            zeros = i - 1;      /* all vectors B[i] with i <= zeros are zero vectors */
            kappa = i + 1;
            kappamax = kappa;

            if (zeros < d - 1)
            {
                fmpz_get_mpf(mpf_mat_entry(r, i, i),
                             fmpz_mat_entry(A->exactSP, i, i));
            }

            for (i = zeros + 1; i < d; i++)
                alpha[i] = 0;

            while (kappa < d)
            {
                int babai_ok = 0;

                if (kappa > kappamax)
                    kappamax = kappa;

                /* ********************************** */
                /* Step3: Call to the Babai algorithm */
                /* ********************************** */
                babai_ok =
                    fmpz_lll_check_babai_heuristic(kappa, B, mu, r, s, appB, A,
                                                   alpha[kappa], zeros,
                                                   kappamax, n, tmp, rtmp,
                                                   prec, fl);

                if (babai_ok == -1)
                {
                    flint_free(alpha);
                    mpf_clears(ctt, tmp, rtmp, '\0');
                    mpf_mat_clear(mu);
                    mpf_mat_clear(r);
                    mpf_mat_clear(appB);
                    fmpz_mat_clear(A->exactSP);
                    _mpf_vec_clear(s, d);
                    _fmpz_vec_clear(exactSPtmp, d);
                    /* Need to switch to mpf / arb */
                    return -1;
                }

                /* ************************************ */
                /* Step4: Success of Lovasz's condition */
                /* ************************************ */

                mpf_mul(tmp, mpf_mat_entry(r, kappa - 1, kappa - 1), ctt);

                if (mpf_cmp(tmp, s + kappa - 1) <= 0)
                {
                    alpha[kappa] = kappa;
                    mpf_mul(tmp, mpf_mat_entry(mu, kappa, kappa - 1),
                            mpf_mat_entry(r, kappa, kappa - 1));
                    mpf_sub(mpf_mat_entry(r, kappa, kappa), s + kappa - 1,
                            tmp);
                    kappa++;
                }
                else
                {

                    /* ******************************************* */
                    /* Step5: Find the right insertion index kappa */
                    /* kappa2 remains the initial kappa            */
                    /* ******************************************* */

                    kappa2 = kappa;
                    do
                    {
                        kappa--;
                        if (kappa > zeros + 1)
                        {
                            mpf_mul(tmp,
                                    mpf_mat_entry(r, kappa - 1, kappa - 1),
                                    ctt);
                        }
                    } while ((kappa >= zeros + 2)
                             && (mpf_cmp(s + kappa - 1, tmp) <= 0));

                    for (i = kappa; i < kappa2; i++)
                        if (kappa <= alpha[i])
                            alpha[i] = kappa;

                    for (i = kappa2; i > kappa; i--)
                        alpha[i] = alpha[i - 1];

                    for (i = kappa2 + 1; i <= kappamax; i++)
                        if (kappa < alpha[i])
                            alpha[i] = kappa;

                    alpha[kappa] = kappa;

                    /* ****************************** */
                    /* Step6: Update the mu's and r's */
                    /* ****************************** */

                    mutmp = mu->rows[kappa2];
                    for (i = kappa2; i > kappa; i--)
                        mu->rows[i] = mu->rows[i - 1];
                    mu->rows[kappa] = mutmp;

                    mutmp = r->rows[kappa2];
                    for (i = kappa2; i > kappa; i--)
                        r->rows[i] = r->rows[i - 1];
                    r->rows[kappa] = mutmp;

                    mpf_set(mpf_mat_entry(r, kappa, kappa), s + kappa);

                    /* ************************ */
                    /* Step7: Update B and appB */
                    /* ************************ */

                    Btmp = B->rows[kappa2];
                    for (i = kappa2; i > kappa; i--)
                        B->rows[i] = B->rows[i - 1];
                    B->rows[kappa] = Btmp;

                    appBtmp = appB->rows[kappa2];
                    for (i = kappa2; i > kappa; i--)
                        appB->rows[i] = appB->rows[i - 1];
                    appB->rows[kappa] = appBtmp;

                    /* ********************* */
                    /* Step8: Update exactSP */
                    /* ********************* */

                    for (i = 0; i <= kappa2; i++)
                        fmpz_set(exactSPtmp + i,
                                 fmpz_mat_entry(A->exactSP, kappa2, i));

                    for (i = kappa2 + 1; i <= kappamax; i++)
                        fmpz_set(exactSPtmp + i,
                                 fmpz_mat_entry(A->exactSP, i, kappa2));

                    for (i = kappa2; i > kappa; i--)
                    {
                        for (j = 0; j < kappa; j++)
                            fmpz_set(fmpz_mat_entry(A->exactSP, i, j),
                                     fmpz_mat_entry(A->exactSP, i - 1, j));
                        fmpz_set(fmpz_mat_entry(A->exactSP, i, kappa),
                                 exactSPtmp + (i - 1));

                        for (j = kappa + 1; j <= i; j++)
                            fmpz_set(fmpz_mat_entry(A->exactSP, i, j),
                                     fmpz_mat_entry(A->exactSP, i - 1, j - 1));

                        for (j = kappa2 + 1; j <= kappamax; j++)
                            fmpz_set(fmpz_mat_entry(A->exactSP, j, i),
                                     fmpz_mat_entry(A->exactSP, j, i - 1));
                    }

                    for (i = 0; i < kappa; i++)
                        fmpz_set(fmpz_mat_entry(A->exactSP, kappa, i),
                                 exactSPtmp + i);
                    fmpz_set(fmpz_mat_entry(A->exactSP, kappa, kappa),
                             exactSPtmp + kappa2);

                    for (i = kappa2 + 1; i <= kappamax; i++)
                        fmpz_set(fmpz_mat_entry(A->exactSP, i, kappa),
                                 exactSPtmp + i);

                    if (mpf_sgn(mpf_mat_entry(r, kappa, kappa)) <= 0.0)
                    {
                        zeros++;
                        kappa++;
                        _fmpz_vec_dot(fmpz_mat_entry(A->exactSP, kappa, kappa),
                                      B->rows[kappa], B->rows[kappa], n);
                        fmpz_get_mpf(mpf_mat_entry(r, kappa, kappa),
                                     fmpz_mat_entry(A->exactSP, kappa, kappa));
                    }

                    kappa++;
                }
            }

            flint_free(alpha);
            mpf_clears(ctt, tmp, rtmp, '\0');
            mpf_mat_clear(mu);
            mpf_mat_clear(r);
            mpf_mat_clear(appB);
            fmpz_mat_clear(A->exactSP);
            _mpf_vec_clear(s, d);
            _fmpz_vec_clear(exactSPtmp, d);
        }
    }
    else
    {
        int kappa, d, i, j, test;
        mpf_mat_t mu, r;
        mpf *s;
        double ctt, halfplus, onedothalfplus;
        mpf_t tmp, rtmp;
        fmpz *x;
        fmpz_t t;
        ulong loops;

        d = B->r;

        mpf_mat_init(mu, d, d, prec);
        mpf_mat_init(r, d, d, prec);

        s = _mpf_vec_init(d, prec);

        ctt = (4 * fl->delta + 1) / 5;
        halfplus = (4 * fl->eta + 0.5) / 5;
        onedothalfplus = 1.0 + halfplus;

        mpf_init2(tmp, prec);
        mpf_init2(rtmp, prec);
        fmpz_init(t);

        i = 0;
        do
            ;
        while ((fmpz_cmp_ui(fmpz_mat_entry(B, i, i), 0) <= 0) && (++i < d));

        kappa = i + 1;

        fmpz_get_mpf(mpf_mat_entry(r, i, i), fmpz_mat_entry(B, i, i));

        loops = 0;
        while (kappa < d)
        {
            do
            {
                test = 0;

                loops++;
                if (loops > 2000)
                {
                    mpf_mat_clear(mu);
                    mpf_mat_clear(r);
                    _mpf_vec_clear(s, d);
                    mpf_clears(tmp, rtmp, '\0');
                    fmpz_clear(t);
                    return -1;
                }

                fmpz_get_mpf(s, fmpz_mat_entry(B, kappa, kappa));
                for (j = 0; j < kappa; j++) /* orthogonalization */
                {
                    fmpz_get_mpf(mpf_mat_entry(r, kappa, j),
                                 fmpz_mat_entry(B, kappa, j));
                    for (i = 0; i < j; i++)
                    {
                        mpf_mul(tmp, mpf_mat_entry(r, kappa, i),
                                mpf_mat_entry(mu, j, i));
                        mpf_sub(mpf_mat_entry(r, kappa, j),
                                mpf_mat_entry(r, kappa, j), tmp);
                    }
                    mpf_div(mpf_mat_entry(mu, kappa, j),
                            mpf_mat_entry(r, kappa, j), mpf_mat_entry(r, j,
                                                                      j));
                    mpf_mul(tmp, mpf_mat_entry(mu, kappa, j),
                            mpf_mat_entry(r, kappa, j));
                    mpf_sub(s + j + 1, s + j, tmp);
                }
                mpf_set(mpf_mat_entry(r, kappa, kappa), s + kappa);

                x = _fmpz_vec_init(kappa);
                for (j = kappa - 1; j >= 0; j--)    /* size-reduction */
                {
                    mpf_abs(tmp, mpf_mat_entry(mu, kappa, j));
                    if (mpf_cmp_d(tmp, halfplus) > 0)
                    {
                        test = 1;
                        if (mpf_cmp_d(tmp, onedothalfplus) <= 0)
                        {
                            if (mpf_sgn(mpf_mat_entry(mu, kappa, j)) >= 0)
                            {
                                mpf_set_si(tmp, 1);
                            }
                            else
                            {
                                mpf_set_si(tmp, -1);
                            }
                        }
                        else
                        {
                            mpf_set_d(rtmp, 0.5);
                            if (mpf_sgn(mpf_mat_entry(mu, kappa, j)) < 0)
                            {
                                mpf_neg(tmp, tmp);
                                mpf_sub(tmp, tmp, rtmp);
                                mpf_ceil(tmp, tmp);
                            }
                            else
                            {
                                mpf_add(tmp, tmp, rtmp);
                                mpf_floor(tmp, tmp);
                            }
                        }
                        fmpz_set_mpf(x + j, tmp);
                        for (i = 0; i < j; i++) /* update μ matrix */
                        {
                            mpf_mul(rtmp, tmp, mpf_mat_entry(mu, j, i));
                            mpf_sub(mpf_mat_entry(mu, kappa, i),
                                    mpf_mat_entry(mu, kappa, i), rtmp);
                        }
                    }
                }

                if (test)
                {
                    for (j = 0; j < kappa; j++)
                    {
                        fmpz_pow_ui(t, x + j, 2);
                        fmpz_addmul(fmpz_mat_entry(B, kappa, kappa), t,
                                    fmpz_mat_entry(B, j, j));

                        fmpz_mul(t, x + j, fmpz_mat_entry(B, kappa, j));
                        fmpz_mul_ui(t, t, 2);
                        fmpz_sub(fmpz_mat_entry(B, kappa, kappa),
                                 fmpz_mat_entry(B, kappa, kappa), t);

                        for (i = 0; i < j; i++)
                        {
                            fmpz_mul(t, x + i, x + j);
                            fmpz_mul(t, t, fmpz_mat_entry(B, j, i));
                            fmpz_mul_ui(t, t, 2);
                            fmpz_add(fmpz_mat_entry(B, kappa, kappa),
                                     fmpz_mat_entry(B, kappa, kappa), t);
                        }
                    }
                    for (i = 0; i < d; i++)
                    {
                        if (i < kappa)
                        {
                            for (j = 0; j <= i; j++)
                                fmpz_submul(fmpz_mat_entry(B, kappa, i), x + j,
                                            fmpz_mat_entry(B, i, j));
                            for (j = i + 1; j < kappa; j++)
                                fmpz_submul(fmpz_mat_entry(B, kappa, i), x + j,
                                            fmpz_mat_entry(B, j, i));
                        }
                        else if (i > kappa)
                        {
                            for (j = 0; j < kappa; j++)
                                fmpz_submul(fmpz_mat_entry(B, i, kappa), x + j,
                                            fmpz_mat_entry(B, i, j));
                        }
                    }
                }

                _fmpz_vec_clear(x, kappa);
            } while (test);

            mpf_set_d(rtmp, ctt);
            mpf_mul(tmp, rtmp, mpf_mat_entry(r, kappa - 1, kappa - 1));
            if (mpf_cmp(tmp, s + kappa - 1) <= 0)   /* check LLL condition */
            {
                mpf_mul(tmp, mpf_mat_entry(mu, kappa, kappa - 1),
                        mpf_mat_entry(r, kappa, kappa - 1));
                mpf_sub(mpf_mat_entry(r, kappa, kappa), s + kappa - 1, tmp);
                kappa++;
            }
            else
            {
                int kappa2;
                kappa2 = kappa;
                do
                {
                    kappa--;
                    if (kappa > 0)
                    {
                        mpf_mul(tmp, rtmp,
                                mpf_mat_entry(r, kappa - 1, kappa - 1));
                    }
                } while ((kappa > 0) && (mpf_cmp(tmp, s + kappa - 1) > 0));

                for (j = 0; j < kappa; j++)
                {
                    mpf_set(mpf_mat_entry(mu, kappa, j),
                            mpf_mat_entry(mu, kappa2, j));
                    mpf_set(mpf_mat_entry(r, kappa, j),
                            mpf_mat_entry(r, kappa2, j));
                }
                mpf_set(mpf_mat_entry(r, kappa, kappa), s + kappa);

                for (j = kappa2; j > kappa; j--)
                {
                    for (i = kappa2; i < d; i++)
                        fmpz_swap(fmpz_mat_entry(B, i, j),
                                  fmpz_mat_entry(B, i, j - 1));
                    for (i = 0; i < kappa; i++)
                        fmpz_swap(fmpz_mat_entry(B, j, i),
                                  fmpz_mat_entry(B, j - 1, i));
                }
                for (j = kappa2; j > kappa; j--)
                {
                    for (i = j; i > kappa; i--)
                        fmpz_swap(fmpz_mat_entry(B, j, i),
                                  fmpz_mat_entry(B, j - 1, i - 1));
                }
                for (j = 0; 2 * j < kappa2 - kappa; j++)
                    fmpz_swap(fmpz_mat_entry(B, kappa + j, kappa),
                              fmpz_mat_entry(B, kappa2 - j, kappa));
                kappa++;
            }
        }

        for (i = 0; i < d - 1; i++)
        {
            for (j = i + 1; j < d; j++)
            {
                fmpz_set(fmpz_mat_entry(B, i, j), fmpz_mat_entry(B, j, i));
            }
        }

        mpf_mat_clear(mu);
        mpf_mat_clear(r);
        _mpf_vec_clear(s, d);
        mpf_clears(tmp, rtmp, '\0');
        fmpz_clear(t);
    }
    return 0;
}
