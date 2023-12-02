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

#if defined(FUNC_HEAD) && defined(TYPE)
#ifdef GM
#undef GM
#endif
#define GM ((fl->rt == Z_BASIS) ? A->exactSP : B)

FUNC_HEAD
{
    int newd = 0;
#if TYPE
    int ok = 1;
    fmpz_t rii;
#endif
    if (fl->rt == Z_BASIS && fl->gt == APPROX)
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

        mpf_init_set_d(ctt, (fl->delta + 1) / 2);

        alpha = (int *) flint_malloc(d * sizeof(int));

        mpf_init2(tmp, prec);
        mpf_init2(rtmp, prec);

        mpf_mat_init(mu, d, d, prec);
        mpf_mat_init(r, d, d, prec);
        mpf_mat_init(appB, d, n, prec);
        mpf_mat_init(A->appSP2, d, d, prec);

        if (U != NULL)
        {
            if (U->r != d)
            {
                flint_throw(FLINT_ERROR, "(fmpz_lll_mpf*): Incompatible dimensions of capturing matrix.\n");
            }
        }

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
            _mpf_vec_set_fmpz_vec(appB->rows[i], B->rows[i], n);

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

        zeros = i - 1;          /* all vectors B[i] with i <= zeros are zero vectors */
        kappa = i + 1;
        kappamax = kappa;

        if (zeros < d - 1)
        {
            mpf_set(mpf_mat_entry(r, i, i), mpf_mat_entry(A->appSP2, i, i));
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
                fmpz_lll_check_babai_heuristic(kappa, B, U, mu, r, s, appB,
                                               A, alpha[kappa], zeros,
                                               kappamax, n, tmp, rtmp,
                                               prec, fl);

            if (babai_ok == -1)
            {
                flint_free(alpha);
                mpf_clears(ctt, tmp, rtmp, NULL);
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
                mpf_sub(mpf_mat_entry(r, kappa, kappa), s + kappa - 1, tmp);
                kappa++;
            }
            else
            {

                /* ******************************************* */
                /* Step5: Find the right insertion index kappa */
                /* kappa2 remains the initial kappa            */
                /* ******************************************* */

                kappa2 = kappa;
#if TYPE
                if (kappa == d - 1 && gs_B != NULL)
                {
                    fmpz_init(rii);
                    mpf_mul(tmp, mpf_mat_entry(mu, kappa, kappa - 1),
                            mpf_mat_entry(r, kappa, kappa - 1));
                    mpf_mul_2exp(tmp, tmp, 1);
                    mpf_sub(tmp, s + kappa - 1, tmp);
                    fmpz_set_mpf(rii, tmp); /* using a heuristic lower bound on the final GS norm */
                    if (fmpz_cmp(rii, gs_B) > 0)
                    {
                        d--;
                    }
                    fmpz_clear(rii);
                    if (kappa >= d)
                    {
                        break;
                    }
                }
#endif
                do
                {
                    kappa--;
                    if (kappa > zeros + 1)
                    {
                        mpf_mul(tmp,
                                mpf_mat_entry(r, kappa - 1, kappa - 1), ctt);
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

                if (U != NULL)
                {
                    Btmp = _fmpz_vec_init(U->c);
                    _fmpz_vec_set(Btmp, U->rows[kappa2], U->c);
                    for (i = kappa2; i > kappa; i--)
                        _fmpz_vec_set(U->rows[i], U->rows[i - 1], U->c);
                    _fmpz_vec_set(U->rows[kappa], Btmp, U->c);
                    _fmpz_vec_clear(Btmp, U->c);
                }

                appBtmp = appB->rows[kappa2];
                for (i = kappa2; i > kappa; i--)
                    appB->rows[i] = appB->rows[i - 1];
                appB->rows[kappa] = appBtmp;

                /* *************************** */
                /* Step8: Update appSP: tricky */
                /* *************************** */

                for (i = 0; i <= kappa2; i++)
                    mpf_set(appSPtmp + i, mpf_mat_entry(A->appSP2, kappa2, i));

                for (i = kappa2 + 1; i <= kappamax; i++)
                    mpf_set(appSPtmp + i, mpf_mat_entry(A->appSP2, i, kappa2));

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
                    mpf_set(mpf_mat_entry(A->appSP2, kappa, i), appSPtmp + i);
                mpf_set(mpf_mat_entry(A->appSP2, kappa, kappa),
                        appSPtmp + kappa2);

                for (i = kappa2 + 1; i <= kappamax; i++)
                    mpf_set(mpf_mat_entry(A->appSP2, i, kappa), appSPtmp + i);

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

#if TYPE
        if (gs_B != NULL)
        {
            newd = d;
            fmpz_init(rii);
            for (i = d - 1; (i >= 0) && (ok > 0); i--)
            {
                mpf_div_2exp(tmp, mpf_mat_entry(r, i, i), 1);
                /* rii is the G-S length of ith vector divided by 2 */
                fmpz_set_mpf(rii, tmp);
                if ((ok = fmpz_cmp(rii, gs_B)) > 0)
                {
                    newd--;
                }
            }
            fmpz_clear(rii);
        }
#endif

        flint_free(alpha);
        mpf_clears(ctt, tmp, rtmp, NULL);
        mpf_mat_clear(mu);
        mpf_mat_clear(r);
        mpf_mat_clear(appB);
        mpf_mat_clear(A->appSP2);
        _mpf_vec_clear(s, B->r);
        _mpf_vec_clear(appSPtmp, B->r);
    }
    else
    {
        int kappa, kappa2, d, n, i, j, zeros, kappamax, update_b = 1;
        mpf_mat_t mu, r;
        fmpz_gram_t A;
        mpf *s, *mutmp;
        mpf_t ctt, tmp, rtmp;
        int *alpha;
        fmpz *Btmp;

        n = B->c;
        d = B->r;

        mpf_init_set_d(ctt, (fl->delta + 1) / 2);

        alpha = (int *) flint_malloc(d * sizeof(int));

        mpf_init2(tmp, prec);
        mpf_init2(rtmp, prec);

        mpf_mat_init(mu, d, d, prec);
        mpf_mat_init(r, d, d, prec);
        if (fl->rt == Z_BASIS)
        {
            fmpz_mat_init(A->exactSP, d, d);
        }

        s = _mpf_vec_init(d, prec);

        if (U != NULL)
        {
            if (U->r != d)
            {
                flint_throw(FLINT_ERROR, "(fmpz_lll_mpf*): Incompatible dimensions of capturing matrix.\n");
            }
            else if (U->c == d && n > d && fmpz_mat_is_one(U))
            {
                update_b = 0;
            }
        }

        if (fl->rt == Z_BASIS)
        {
            fmpz_mat_gram(A->exactSP, B);
        }

        /* ********************************* */
        /* Step2: Initializing the main loop */
        /* ********************************* */

        kappamax = 0;
        i = 0;

        do
            ;
        while ((fmpz_cmp_ui(fmpz_mat_entry(GM, i, i), 0) <= 0) && (++i < d));

        zeros = i - 1;          /* all vectors B[i] with i <= zeros are zero vectors */
        kappa = i + 1;
        kappamax = kappa;

        if (zeros < d - 1)
        {
            fmpz_get_mpf(mpf_mat_entry(r, i, i), fmpz_mat_entry(GM, i, i));
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
                fmpz_lll_check_babai_heuristic(kappa, (update_b ? B : NULL), U,
                                               mu, r, s, NULL, A, alpha[kappa],
                                               zeros, kappamax, n, tmp, rtmp,
                                               prec, fl);

            if (babai_ok == -1)
            {
                flint_free(alpha);
                mpf_clears(ctt, tmp, rtmp, NULL);
                mpf_mat_clear(mu);
                mpf_mat_clear(r);
                if (fl->rt == Z_BASIS)
                {
                    fmpz_mat_clear(A->exactSP);
                }
                _mpf_vec_clear(s, d);
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
                mpf_sub(mpf_mat_entry(r, kappa, kappa), s + kappa - 1, tmp);
                kappa++;
            }
            else
            {

                /* ******************************************* */
                /* Step5: Find the right insertion index kappa */
                /* kappa2 remains the initial kappa            */
                /* ******************************************* */

                kappa2 = kappa;
#if TYPE
                if (kappa == d - 1 && gs_B != NULL)
                {
                    fmpz_init(rii);
                    mpf_mul(tmp, mpf_mat_entry(mu, kappa, kappa - 1),
                            mpf_mat_entry(r, kappa, kappa - 1));
                    mpf_mul_2exp(tmp, tmp, 1);
                    mpf_sub(tmp, s + kappa - 1, tmp);
                    fmpz_set_mpf(rii, tmp); /* using a heuristic lower bound on the final GS norm */
                    if (fmpz_cmp(rii, gs_B) > 0)
                    {
                        d--;
                    }
                    fmpz_clear(rii);
                    if (kappa >= d)
                    {
                        break;
                    }
                }
#endif
                do
                {
                    kappa--;
                    if (kappa > zeros + 1)
                    {
                        mpf_mul(tmp,
                                mpf_mat_entry(r, kappa - 1, kappa - 1), ctt);
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

                /* *************** */
                /* Step7: Update B */
                /* *************** */

                if (fl->rt == Z_BASIS && update_b)
                {
                    Btmp = B->rows[kappa2];
                    for (i = kappa2; i > kappa; i--)
                        B->rows[i] = B->rows[i - 1];
                    B->rows[kappa] = Btmp;
                }

                if (U != NULL)
                {
                    Btmp = _fmpz_vec_init(U->c);
                    _fmpz_vec_set(Btmp, U->rows[kappa2], U->c);
                    for (i = kappa2; i > kappa; i--)
                        _fmpz_vec_set(U->rows[i], U->rows[i - 1], U->c);
                    _fmpz_vec_set(U->rows[kappa], Btmp, U->c);
                    _fmpz_vec_clear(Btmp, U->c);
                }

                /* ********************* */
                /* Step8: Update exactSP */
                /* ********************* */

                for (j = kappa2; j > kappa; j--)
                {
                    for (i = kappa2; i < d; i++)
                        fmpz_swap(fmpz_mat_entry(GM, i, j),
                                  fmpz_mat_entry(GM, i, j - 1));
                    for (i = 0; i < kappa; i++)
                        fmpz_swap(fmpz_mat_entry(GM, j, i),
                                  fmpz_mat_entry(GM, j - 1, i));
                }
                for (j = kappa2; j > kappa; j--)
                {
                    for (i = j; i > kappa; i--)
                        fmpz_swap(fmpz_mat_entry(GM, j, i),
                                  fmpz_mat_entry(GM, j - 1, i - 1));
                }
                for (j = 0; 2 * j < kappa2 - kappa; j++)
                    fmpz_swap(fmpz_mat_entry(GM, kappa + j, kappa),
                              fmpz_mat_entry(GM, kappa2 - j, kappa));

                if (mpf_sgn(mpf_mat_entry(r, kappa, kappa)) <= 0.0)
                {
                    zeros++;
                    kappa++;
                    fmpz_get_mpf(mpf_mat_entry(r, kappa, kappa),
                                 fmpz_mat_entry(GM, kappa, kappa));
                }

                kappa++;
            }
        }

        if (fl->rt == GRAM)
        {
            for (i = 0; i < B->r - 1; i++)
            {
                for (j = i + 1; j < B->r; j++)
                {
                    fmpz_set(fmpz_mat_entry(B, i, j), fmpz_mat_entry(B, j, i));
                }
            }
        }
        else if (!update_b)
        {
            fmpz_mat_mul(B, U, B);
        }

#if TYPE
        if (gs_B != NULL)
        {
            newd = d;
            fmpz_init(rii);
            for (i = d - 1; (i >= 0) && (ok > 0); i--)
            {
                mpf_div_2exp(tmp, mpf_mat_entry(r, i, i), 1);
                /* rii is the G-S length of ith vector divided by 2 */
                fmpz_set_mpf(rii, tmp);
                if ((ok = fmpz_cmp(rii, gs_B)) > 0)
                {
                    newd--;
                }
            }
            fmpz_clear(rii);
        }
#endif

        flint_free(alpha);
        mpf_clears(ctt, tmp, rtmp, NULL);
        mpf_mat_clear(mu);
        mpf_mat_clear(r);
        if (fl->rt == Z_BASIS)
        {
            fmpz_mat_clear(A->exactSP);
        }
        _mpf_vec_clear(s, B->r);
    }
    return newd;
}

#undef GM

#endif
