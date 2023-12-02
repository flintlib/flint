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

#if defined(FUNC_HEAD) && defined(CALL_BABAI) && defined(TYPE)
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
    int *expo = NULL;
    if (fl->rt == Z_BASIS && fl->gt == APPROX)
    {
        int kappa, kappa2, d, n, i, j, zeros, kappamax, shift;
        int num_failed_fast = 0;
        int babai_ok = 0;
        int heuristic_fail = 0;
        d_mat_t mu, r, appB;
        fmpz_gram_t A;
        double *s, *mutmp, *appBtmp, *appSPtmp;
        double ctt;
        int *alpha;
        fmpz *Btmp;
        ulong max_exp, iter, max_iter, newvec, newvec_max;

        n = B->c;
        d = B->r;

        ctt = (fl->delta + 1) / 2;

        shift = fmpz_lll_shift(B);

        alpha = (int *) flint_malloc(d * sizeof(int));
        expo = (int *) flint_malloc(d * sizeof(int));

        d_mat_init(mu, d, d);
        d_mat_init(r, d, d);
        d_mat_init(appB, d, n);
        d_mat_init(A->appSP, d, d);

        s = _d_vec_init(d);
        appSPtmp = _d_vec_init(d);

        if (U != NULL)
        {
            if (U->r != d)
            {
                flint_throw(FLINT_ERROR, "(fmpz_lll_d*): Incompatible dimensions of capturing matrix.\n");
            }
        }

        for (i = 0; i < d; i++)
        {
            for (j = 0; j < d; j++)
            {
                d_mat_entry(A->appSP, i, j) = D_NAN;
            }
        }

        /* ************************** */
        /* Step1: Initialization Step */
        /* ************************** */

        max_exp = 0;
        for (i = 0; i < d; i++)
        {
            expo[i] = _fmpz_vec_get_d_vec_2exp(appB->rows[i], B->rows[i], n);
            max_exp = FLINT_MAX(max_exp, expo[i]);
        }
        max_iter =
            (ulong) ((d - 1) +
                     (d - 1) * d * (2 * max_exp +
                                    d_log2(d)) / d_log2(8 / (fl->delta + 7)));

        /* ********************************* */
        /* Step2: Initializing the main loop */
        /* ********************************* */

        kappamax = 0;
        i = 0;

        do
        {
            d_mat_entry(A->appSP, i, i) = _d_vec_norm(appB->rows[i], n);
        } while ((d_mat_entry(A->appSP, i, i) <= 0.0) && (++i < d));    /* Check if this should be D_EPS and not 0.0: done */

        zeros = i - 1;          /* all vectors B[i] with i <= zeros are zero vectors */
        kappa = i + 1;
        kappamax = kappa;

        if (zeros < d - 1)
        {
            d_mat_entry(r, i, i) = d_mat_entry(A->appSP, i, i);
        }

        for (i = zeros + 1; i < d; i++)
            alpha[i] = 0;

        newvec = 0;
        newvec_max = 1;
        iter = 0;
        while (kappa < d)
        {
            int new_kappa;
            double tmp = 0.0;

            if (iter >= max_iter)
            {
                break;
            }
            iter++;

            new_kappa = 0;
            if (kappa > kappamax)
            {
                /* In the first time we hit a new kappa we're going to size-reduce in advance (for knapsack)... */
                kappamax = kappa;
                newvec++;

                if (newvec > newvec_max)
                {
                    newvec_max *= 2;
                    newvec = 0;
                    new_kappa = 1;
                }
            }

            /* ********************************** */
            /* Step3: Call to the Babai algorithm */
            /* ********************************** */
            CALL_BABAI(num_failed_fast, babai_ok, heuristic_fail);

            if (heuristic_fail == -1)
            {
                flint_free(alpha);
                flint_free(expo);
                d_mat_clear(mu);
                d_mat_clear(r);
                d_mat_clear(appB);
                d_mat_clear(A->appSP);
                _d_vec_clear(s);
                _d_vec_clear(appSPtmp);
                /* Need to switch to mpf / arb */
                return -1;
            }

            /* End of the real Babai part... */
            if (new_kappa == 1)
            {
#if TYPE == 2
                /* running ahead to kappa = d, without upsetting LLL... */
                for (kappa2 = d - 1; kappa2 > kappa; kappa2--)
                {
                    babai_ok =
                        fmpz_lll_advance_check_babai(kappa, kappa2, B, U,
                                                     mu, r, s, appB, expo,
                                                     A, alpha[kappa2],
                                                     zeros, kappa + 1, n, fl);
                    if (babai_ok == -1)
                    {
                        heuristic_fail =
                            fmpz_lll_advance_check_babai_heuristic_d(kappa,
                                                                     kappa2,
                                                                     B, U,
                                                                     mu, r,
                                                                     s,
                                                                     appB,
                                                                     expo,
                                                                     A,
                                                                     alpha
                                                                     [kappa2],
                                                                     zeros,
                                                                     kappa
                                                                     + 1,
                                                                     n, fl);
                    }
                }
#endif
            }

            /* ************************************ */
            /* Step4: Success of Lovasz's condition */
            /* ************************************ */

            tmp = d_mat_entry(r, kappa - 1, kappa - 1) * ctt;
            tmp = ldexp(tmp, 2 * (expo[kappa - 1] - expo[kappa]));

            if (tmp <= s[kappa - 1])
            {
                alpha[kappa] = kappa;
                tmp =
                    d_mat_entry(mu, kappa, kappa - 1) * d_mat_entry(r,
                                                                    kappa,
                                                                    kappa - 1);
                d_mat_entry(r, kappa, kappa) = s[kappa - 1] - tmp;
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
                    tmp =
                        2 * d_mat_entry(mu, kappa,
                                        kappa - 1) * d_mat_entry(r, kappa,
                                                                 kappa - 1);
                    fmpz_set_d_2exp(rii, s[kappa - 1] - tmp, 2 * expo[kappa]);    /* using a heuristic lower bound on the final GS norm */
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
                        tmp = d_mat_entry(r, kappa - 1, kappa - 1) * ctt;
                        tmp = ldexp(tmp, 2 * (expo[kappa - 1] - expo[kappa2]));
                    }
                } while ((kappa >= zeros + 2) && (s[kappa - 1] <= tmp));

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

                d_mat_entry(r, kappa, kappa) = s[kappa];

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

                j = expo[kappa2];
                for (i = kappa2; i > kappa; i--)
                    expo[i] = expo[i - 1];
                expo[kappa] = j;

                /* *************************** */
                /* Step8: Update appSP: tricky */
                /* *************************** */

                for (i = 0; i <= kappa2; i++)
                    appSPtmp[i] = d_mat_entry(A->appSP, kappa2, i);

                for (i = kappa2 + 1; i <= kappamax; i++)
                    appSPtmp[i] = d_mat_entry(A->appSP, i, kappa2);

                for (i = kappa2; i > kappa; i--)
                {
                    for (j = 0; j < kappa; j++)
                        d_mat_entry(A->appSP, i, j) =
                            d_mat_entry(A->appSP, i - 1, j);
                    d_mat_entry(A->appSP, i, kappa) = appSPtmp[i - 1];

                    for (j = kappa + 1; j <= i; j++)
                        d_mat_entry(A->appSP, i, j) =
                            d_mat_entry(A->appSP, i - 1, j - 1);

                    for (j = kappa2 + 1; j <= kappamax; j++)
                        d_mat_entry(A->appSP, j, i) =
                            d_mat_entry(A->appSP, j, i - 1);
                }

                for (i = 0; i < kappa; i++)
                    d_mat_entry(A->appSP, kappa, i) = appSPtmp[i];
                d_mat_entry(A->appSP, kappa, kappa) = appSPtmp[kappa2];

                for (i = kappa2 + 1; i <= kappamax; i++)
                    d_mat_entry(A->appSP, i, kappa) = appSPtmp[i];

                if (d_mat_entry(r, kappa, kappa) <= 0.0)
                {
                    zeros++;
                    kappa++;
                    d_mat_entry(A->appSP, kappa, kappa) =
                        _d_vec_norm(appB->rows[kappa], n);
                    d_mat_entry(r, kappa, kappa) =
                        d_mat_entry(A->appSP, kappa, kappa);
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
                /* rii is the G-S length of ith vector divided by 2 */
                fmpz_set_d_2exp(rii, d_mat_entry(r, i, i), 2 * expo[i] - 1);
                if ((ok = fmpz_cmp(rii, gs_B)) > 0)
                {
                    newd--;
                }
            }
            fmpz_clear(rii);
        }
#endif

        flint_free(alpha);
        flint_free(expo);
        d_mat_clear(mu);
        d_mat_clear(r);
        d_mat_clear(appB);
        d_mat_clear(A->appSP);
        _d_vec_clear(s);
        _d_vec_clear(appSPtmp);

        if (kappa < d)
            return -1;
    }
    else
    {
        int kappa, kappa2, d, n, i, j, zeros, kappamax, update_b = 1;
        slong exp;
        int num_failed_fast = 0;
        int babai_ok = 0;
        int heuristic_fail = 0;
        d_mat_t mu, r;
        fmpz_gram_t A;
        double *s, *mutmp;
        double ctt;
        int *alpha;
        fmpz *Btmp;
        ulong max_exp, iter, max_iter;

        n = B->c;
        d = B->r;

        ctt = (fl->delta + 1) / 2;

        alpha = (int *) flint_malloc(d * sizeof(int));
        expo = (int *) flint_malloc(d * sizeof(int));

        d_mat_init(mu, d, d);
        d_mat_init(r, d, d);
        if (fl->rt == Z_BASIS)
        {
            fmpz_mat_init(A->exactSP, d, d);
        }

        s = _d_vec_init(d);

        if (U != NULL)
        {
            if (U->r != d)
            {
                flint_throw(FLINT_ERROR, "(fmpz_lll_d*): Incompatible dimensions of capturing matrix.\n");
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

        /* ************************** */
        /* Step1: Initialization Step */
        /* ************************** */

        max_exp = 0;
        for (i = 0; i < d; i++)
        {
            fmpz_get_d_2exp(&exp, fmpz_mat_entry(GM, i, i));
            max_exp = FLINT_MAX(max_exp, (expo[i] = exp));
        }
        max_iter =
            (ulong) ((d - 1) +
                     (d - 1) * d * max_exp / d_log2(8 / (fl->delta + 7)));

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
            d_mat_entry(r, i, i) =
                fmpz_get_d_2exp(&exp, fmpz_mat_entry(GM, i, i));
            d_mat_entry(r, i, i) = ldexp(d_mat_entry(r, i, i), exp - expo[i]);
        }

        for (i = zeros + 1; i < d; i++)
            alpha[i] = 0;

        iter = 0;
        while (kappa < d)
        {
            double tmp = 0.0;
            if (iter >= max_iter)
            {
                break;
            }
            iter++;

            if (kappa > kappamax)
            {
                kappamax = kappa;
            }

            /* ********************************** */
            /* Step3: Call to the Babai algorithm */
            /* ********************************** */
            if (num_failed_fast < 50)
            {
                babai_ok =
                    fmpz_lll_check_babai(kappa, (update_b ? B : NULL), U, mu,
                                         r, s, NULL, expo, A, alpha[kappa],
                                         zeros, kappamax, n, fl);
            }
            else
            {
                babai_ok = -1;
            }

            if (babai_ok == -1)
            {
                num_failed_fast++;
                heuristic_fail =
                    fmpz_lll_check_babai_heuristic_d(kappa,
                                                     (update_b ? B : NULL), U,
                                                     mu, r, s, NULL, expo, A,
                                                     alpha[kappa], zeros,
                                                     kappamax, n, fl);
            }

            if (heuristic_fail == -1)
            {
                flint_free(alpha);
                flint_free(expo);
                d_mat_clear(mu);
                d_mat_clear(r);
                if (fl->rt == Z_BASIS)
                {
                    fmpz_mat_clear(A->exactSP);
                }
                _d_vec_clear(s);
                /* Need to switch to mpf / arb */
                return -1;
            }

            /* ************************************ */
            /* Step4: Success of Lovasz's condition */
            /* ************************************ */

            tmp = d_mat_entry(r, kappa - 1, kappa - 1) * ctt;
            tmp = ldexp(tmp, (expo[kappa - 1] - expo[kappa]));

            if (tmp <= s[kappa - 1])
            {
                alpha[kappa] = kappa;
                tmp =
                    ldexp(d_mat_entry(mu, kappa, kappa - 1) * d_mat_entry(r,
                                                                          kappa,
                                                                          kappa
                                                                          - 1),
                          (expo[kappa] - expo[kappa - 1]));
                d_mat_entry(r, kappa, kappa) = s[kappa - 1] - tmp;
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
                    tmp =
                        ldexp(2 * d_mat_entry(mu, kappa,
                                              kappa - 1) * d_mat_entry(r,
                                                                       kappa,
                                                                       kappa -
                                                                       1),
                              (expo[kappa] - expo[kappa - 1]));
                    fmpz_set_d(rii, ldexp(s[kappa - 1] - tmp, expo[kappa]));    /* using a heuristic lower bound on the final GS norm */
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
                        tmp = d_mat_entry(r, kappa - 1, kappa - 1) * ctt;
                        tmp = ldexp(tmp, (expo[kappa - 1] - expo[kappa2]));
                    }
                } while ((kappa >= zeros + 2) && (s[kappa - 1] <= tmp));

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

                d_mat_entry(r, kappa, kappa) = s[kappa];

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

                j = expo[kappa2];
                for (i = kappa2; i > kappa; i--)
                    expo[i] = expo[i - 1];
                expo[kappa] = j;

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

                if (d_mat_entry(r, kappa, kappa) <= 0.0)
                {
                    zeros++;
                    kappa++;
                    d_mat_entry(r, kappa, kappa) =
                        fmpz_get_d_2exp(&exp, fmpz_mat_entry
                                        (GM, kappa, kappa));
                    d_mat_entry(r, kappa, kappa) =
                        ldexp(d_mat_entry(r, kappa, kappa), exp - expo[kappa]);
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
                /* rii is the G-S length of ith vector divided by 2 */
                fmpz_set_d(rii, ldexp(d_mat_entry(r, i, i), expo[i] - 1));
                if ((ok = fmpz_cmp(rii, gs_B)) > 0)
                {
                    newd--;
                }
            }
            fmpz_clear(rii);
        }
#endif

        flint_free(alpha);
        flint_free(expo);
        d_mat_clear(mu);
        d_mat_clear(r);
        if (fl->rt == Z_BASIS)
        {
            fmpz_mat_clear(A->exactSP);
        }
        _d_vec_clear(s);

        if (kappa < d)
            return -1;
    }
    return newd;
}

#undef GM

#endif
