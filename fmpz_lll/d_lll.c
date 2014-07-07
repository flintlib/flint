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

#if defined(FUNC_HEAD) && defined(CALL_BABAI) && defined(TYPE)

FUNC_HEAD
{
    int newd = 0, ok = 1;
    fmpz_t rii;
    int *expo = NULL;
    if (fl->rt == Z_BASIS)
    {
        if (fl->gt == APPROX)
        {
            int kappa, kappa2, d, n, i, j, zeros, kappamax, shift, new_kappa;
            int num_failed_fast = 0;
            int babai_ok = 0;
            int heuristic_fail = 0;
            d_mat_t mu, r, appB;
            fmpz_gram_t A;
            double *s, *mutmp, *appBtmp, *appSPtmp;
            double ctt, tmp = 0.0;
            int *alpha;
            fmpz *Btmp;
            ulong max_exp, iter, max_iter, newvec, newvec_max;

            n = B->c;
            d = B->r;

            ctt = (4 * fl->delta + 1) / 5;

            shift = fmpz_lll_shift(B);

            alpha = (int *) flint_malloc(d * sizeof(int));
            expo = (int *) flint_malloc(d * sizeof(int));

            d_mat_init(mu, d, d);
            d_mat_init(r, d, d);
            d_mat_init(appB, d, n);
            d_mat_init(A->appSP, d, d);

            s = _d_vec_init(d);
            appSPtmp = _d_vec_init(d);

            if (fl->store_trans)
            {
                fmpz_mat_one(U);
            }

            for (i = 0; i < d; i++)
            {
                for (j = 0; j < d; j++)
                {
                    d_mat_entry(A->appSP, i, j) = NAN;
                }
            }

            /* ************************** */
            /* Step1: Initialization Step */
            /* ************************** */

            max_exp = 0;
            for (i = 0; i < d; i++)
            {
                expo[i] =
                    _fmpz_vec_get_d_vec_2exp(appB->rows[i], B->rows[i], n);
                max_exp = FLINT_MAX(max_exp, expo[i]);
            }
            max_iter =
                (ulong) (d - 2 * d * (d + 1) * ((max_exp + 3) / log(ctt)));

            /* ********************************* */
            /* Step2: Initializing the main loop */
            /* ********************************* */

            kappamax = 0;
            i = 0;

            do
            {
                d_mat_entry(A->appSP, i, i) = _d_vec_norm(appB->rows[i], n);
            } while ((d_mat_entry(A->appSP, i, i) <= 0.0) && (++i < d));    /* Check if this should be D_EPS and not 0.0: done */

            zeros = i - 1;      /* all vectors B[i] with i <= zeros are zero vectors */
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
                                                         zeros, kappa + 1,
                                                         n, fl);
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
                                                                         n,
                                                                         fl);
                        }
                    }
#endif
                }

#if TYPE
                /* Use the newd stuff here... */
                ok = 1;
                if (kappa == d - 1
                    && ldexp(d_mat_entry(r, kappa, kappa),
                             2 * (expo[kappa] - expo[kappa - 1])) >=
                    d_mat_entry(r, kappa - 1, kappa - 1) / 2.0)
                {
                    fmpz_init(rii);
                    for (i = d - 1; (i >= 0) && (ok > 0); i--)
                    {
                        /* rii is the squared G-S length of ith vector */
                        fmpz_set_d(rii,
                                   ldexp(d_mat_entry(r, i, i), 2 * expo[i]));
                        if ((ok = fmpz_cmp(rii, gs_B)) > 0)
                        {
                            d--;
                        }
                    }
                    fmpz_clear(rii);
                }
                newd = d;
                if (kappa >= d)
                {
                    break;
                }
#endif

                /* ************************************ */
                /* Step4: Success of Lovasz's condition */
                /* ************************************ */

                tmp = d_mat_entry(r, kappa - 1, kappa - 1) * ctt;
                tmp = ldexp(tmp, 2 * (expo[kappa - 1] - expo[kappa]));

                if (tmp <= s[kappa - 1])
                {
                    alpha[kappa] = kappa;
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
                            tmp = d_mat_entry(r, kappa - 1, kappa - 1) * ctt;
                            tmp =
                                ldexp(tmp,
                                      2 * (expo[kappa - 1] - expo[kappa2]));
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

                    Btmp = _fmpz_vec_init(d);
                    if (fl->store_trans)
                    {
                        _fmpz_vec_set(Btmp, U->rows[kappa2], d);
                        for (i = kappa2; i > kappa; i--)
                            _fmpz_vec_set(U->rows[i], U->rows[i - 1], d);
                        _fmpz_vec_set(U->rows[kappa], Btmp, d);
                    }
                    _fmpz_vec_clear(Btmp, d);

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
            int kappa, kappa2, d, n, i, j, zeros, kappamax, shift, new_kappa;
            slong exp;
            int num_failed_fast = 0;
            int babai_ok = 0;
            int heuristic_fail = 0;
            d_mat_t mu, r, appB;
            fmpz_gram_t A;
            double *s, *mutmp, *appBtmp;
            fmpz *exactSPtmp;
            double ctt;
            int *alpha;
            fmpz *Btmp;
            ulong max_exp, iter, max_iter, newvec, newvec_max;

            n = B->c;
            d = B->r;

            ctt = (4 * fl->delta + 1) / 5;

            shift = fmpz_lll_shift(B);

            alpha = (int *) flint_malloc(d * sizeof(int));
            expo = (int *) flint_malloc(d * sizeof(int));

            d_mat_init(mu, d, d);
            d_mat_init(r, d, d);
            d_mat_init(appB, d, n);
            fmpz_mat_init(A->exactSP, d, d);

            s = _d_vec_init(d);
            exactSPtmp = _fmpz_vec_init(d);

            if (fl->store_trans)
            {
                fmpz_mat_one(U);
            }

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

            max_exp = 0;
            for (i = 0; i < d; i++)
            {
                expo[i] =
                    _fmpz_vec_get_d_vec_2exp(appB->rows[i], B->rows[i], n);
                max_exp = FLINT_MAX(max_exp, expo[i]);
            }
            max_iter =
                (ulong) (d - 2 * d * (d + 1) * ((max_exp + 3) / log(ctt)));

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
                d_mat_entry(r, i, i) =
                    fmpz_get_d_2exp(&exp, fmpz_mat_entry(A->exactSP, i, i));
                d_mat_entry(r, i, i) =
                    ldexp(d_mat_entry(r, i, i), exp - 2 * expo[i]);
            }

            for (i = zeros + 1; i < d; i++)
                alpha[i] = 0;

            newvec = 0;
            newvec_max = 1;
            iter = 0;
            while (kappa < d)
            {
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
                if (num_failed_fast < 50)
                {
                    babai_ok =
                        fmpz_lll_check_babai(kappa, B, U, mu, r, s, appB, expo,
                                             A, alpha[kappa], zeros, kappamax,
                                             FLINT_MIN(kappamax + 1 + shift,
                                                       n), fl);
                }
                else
                {
                    babai_ok = -1;
                }

                if (babai_ok == -1)
                {
                    num_failed_fast++;
                    heuristic_fail =
                        fmpz_lll_check_babai_heuristic_d(kappa, B, U, mu, r, s,
                                                         appB, expo, A,
                                                         alpha[kappa], zeros,
                                                         kappamax,
                                                         FLINT_MIN(kappamax +
                                                                   1 + shift,
                                                                   n), fl);
                }

                if (heuristic_fail == -1)
                {
                    flint_free(alpha);
                    flint_free(expo);
                    d_mat_clear(mu);
                    d_mat_clear(r);
                    d_mat_clear(appB);
                    fmpz_mat_clear(A->exactSP);
                    _d_vec_clear(s);
                    _fmpz_vec_clear(exactSPtmp, d);
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
                                                         zeros, kappa + 1,
                                                         FLINT_MIN(kappamax +
                                                                   1 + shift,
                                                                   n), fl);
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
                                                                         FLINT_MIN
                                                                         (kappamax
                                                                          + 1 +
                                                                          shift,
                                                                          n),
                                                                         fl);
                        }
                    }
#endif
                }

#if TYPE
                /* Use the newd stuff here... */
                ok = 1;
                if (kappa == d - 1
                    && ldexp(d_mat_entry(r, kappa, kappa),
                             2 * (expo[kappa] - expo[kappa - 1])) >=
                    d_mat_entry(r, kappa - 1, kappa - 1) / 2.0)
                {
                    fmpz_init(rii);
                    for (i = d - 1; (i >= 0) && (ok > 0); i--)
                    {
                        /* rii is the squared G-S length of ith vector */
                        fmpz_set_d(rii,
                                   ldexp(d_mat_entry(r, i, i), 2 * expo[i]));
                        if ((ok = fmpz_cmp(rii, gs_B)) > 0)
                        {
                            d--;
                        }
                    }
                    fmpz_clear(rii);
                }
                newd = d;
                if (kappa >= d)
                {
                    break;
                }
#endif

                /* ************************************ */
                /* Step4: Success of Lovasz's condition */
                /* ************************************ */

                tmp = d_mat_entry(r, kappa - 1, kappa - 1) * ctt;
                tmp = ldexp(tmp, 2 * (expo[kappa - 1] - expo[kappa]));

                if (tmp <= s[kappa - 1])
                {
                    alpha[kappa] = kappa;
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
                            tmp = d_mat_entry(r, kappa - 1, kappa - 1) * ctt;
                            tmp =
                                ldexp(tmp,
                                      2 * (expo[kappa - 1] - expo[kappa2]));
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

                    Btmp = _fmpz_vec_init(d);
                    if (fl->store_trans)
                    {
                        _fmpz_vec_set(Btmp, U->rows[kappa2], d);
                        for (i = kappa2; i > kappa; i--)
                            _fmpz_vec_set(U->rows[i], U->rows[i - 1], d);
                        _fmpz_vec_set(U->rows[kappa], Btmp, d);
                    }
                    _fmpz_vec_clear(Btmp, d);

                    appBtmp = appB->rows[kappa2];
                    for (i = kappa2; i > kappa; i--)
                        appB->rows[i] = appB->rows[i - 1];
                    appB->rows[kappa] = appBtmp;

                    j = expo[kappa2];
                    for (i = kappa2; i > kappa; i--)
                        expo[i] = expo[i - 1];
                    expo[kappa] = j;

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

                    if (d_mat_entry(r, kappa, kappa) <= 0.0)
                    {
                        zeros++;
                        kappa++;
                        _fmpz_vec_dot(fmpz_mat_entry(A->exactSP, kappa, kappa),
                                      B->rows[kappa], B->rows[kappa], n);
                        d_mat_entry(r, kappa, kappa) =
                            fmpz_get_d_2exp(&exp, fmpz_mat_entry
                                            (A->exactSP, kappa, kappa));
                        d_mat_entry(r, kappa, kappa) =
                            ldexp(d_mat_entry(r, kappa, kappa),
                                  exp - 2 * expo[kappa]);
                    }

                    kappa++;
                }
            }

            flint_free(alpha);
            flint_free(expo);
            d_mat_clear(mu);
            d_mat_clear(r);
            d_mat_clear(appB);
            fmpz_mat_clear(A->exactSP);
            _d_vec_clear(s);
            _fmpz_vec_clear(exactSPtmp, d);

            if (kappa < d)
                return -1;
        }
    }
    else
    {
        int kappa, d, i, j, test, max_expo = INT_MAX;
        d_mat_t mu, r;
        double *s;
        double ctt, halfplus, onedothalfplus;
        fmpz *x;
        fmpz_t t, dmax;
        ulong loops;

        d = B->r;

        d_mat_init(mu, d, d);
        d_mat_init(r, d, d);

        s = _d_vec_init(d);

        ctt = (4 * fl->delta + 1) / 5;
        halfplus = (4 * fl->eta + 0.5) / 5;
        onedothalfplus = 1.0 + halfplus;

        fmpz_init(t);
        fmpz_init(dmax);
        fmpz_set_d(dmax, DBL_MAX);

        if (fl->store_trans)
        {
            fmpz_mat_one(U);
        }

        i = 0;
        do
            ;
        while ((fmpz_cmp_ui(fmpz_mat_entry(B, i, i), 0) <= 0) && (++i < d));

        kappa = i + 1;

        if (fmpz_cmpabs(fmpz_mat_entry(B, i, i), dmax) > 0)
        {
            d_mat_clear(mu);
            d_mat_clear(r);
            _d_vec_clear(s);
            fmpz_clear(t);
            fmpz_clear(dmax);
            return -1;
        }
        d_mat_entry(r, i, i) = fmpz_get_d(fmpz_mat_entry(B, i, i));

        while (kappa < d)
        {
            loops = 0;
            do
            {
                test = 0;

                for (j = 0; j < kappa; j++) /* orthogonalization */
                {
                    if (fmpz_cmpabs(fmpz_mat_entry(B, kappa, j), dmax) > 0)
                    {
                        d_mat_clear(mu);
                        d_mat_clear(r);
                        _d_vec_clear(s);
                        fmpz_clear(t);
                        fmpz_clear(dmax);
                        return -1;
                    }
                    d_mat_entry(r, kappa, j) =
                        fmpz_get_d(fmpz_mat_entry(B, kappa, j));
                    for (i = 0; i < j; i++)
                    {
                        d_mat_entry(r, kappa, j) -=
                            d_mat_entry(r, kappa, i) * d_mat_entry(mu, j, i);
                    }
                    d_mat_entry(mu, kappa, j) =
                        d_mat_entry(r, kappa, j) / d_mat_entry(r, j, j);
                }

                if (loops >= 20)
                {
                    int new_max_expo = INT_MIN;
                    for (j = 0; j < kappa; j++)
                    {
                        int expo2;
                        frexp(d_mat_entry(mu, kappa, j), &expo2);
                        new_max_expo = FLINT_MAX(new_max_expo, expo2);
                    }
                    if (new_max_expo > max_expo - SIZE_RED_FAILURE_THRESH)
                    {
                        d_mat_clear(mu);
                        d_mat_clear(r);
                        _d_vec_clear(s);
                        fmpz_clear(t);
                        fmpz_clear(dmax);
                        return -1;
                    }
                    max_expo = new_max_expo;
                }

                x = _fmpz_vec_init(kappa);
                for (j = kappa - 1; j >= 0; j--)    /* size-reduction */
                {
                    double tmp = fabs(d_mat_entry(mu, kappa, j));
                    if (tmp > halfplus)
                    {
                        test = 1;
                        if (tmp <= onedothalfplus)
                        {
                            if (d_mat_entry(mu, kappa, j) >= 0)
                            {
                                tmp = 1;
                            }
                            else
                            {
                                tmp = -1;
                            }
                        }
                        else
                        {
                            tmp = d_mat_entry(mu, kappa, j);
                            if (tmp < 0)
                                tmp = ceil(tmp - 0.5);
                            else
                                tmp = floor(tmp + 0.5);
                        }
                        fmpz_set_d(x + j, tmp);
                        for (i = 0; i < j; i++) /* update μ matrix */
                        {
                            d_mat_entry(mu, kappa, i) -=
                                tmp * d_mat_entry(mu, j, i);
                        }
                        if (fl->store_trans)
                        {
                            _fmpz_vec_scalar_submul_fmpz(U->rows[kappa],
                                                         U->rows[j], d, x + j);
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
                loops++;
            } while (test);

            if (fmpz_cmpabs(fmpz_mat_entry(B, kappa, kappa), dmax) > 0)
            {
                d_mat_clear(mu);
                d_mat_clear(r);
                _d_vec_clear(s);
                fmpz_clear(t);
                fmpz_clear(dmax);
                return -1;
            }
            s[0] = fmpz_get_d(fmpz_mat_entry(B, kappa, kappa));
            for (j = 0; j < kappa; j++)
            {
                s[j + 1] =
                    s[j] - d_mat_entry(mu, kappa, j) * d_mat_entry(r,
                                                                   kappa, j);
            }
            d_mat_entry(r, kappa, kappa) = s[kappa];

#if TYPE
            /* Use the newd stuff here... */
            ok = 1;
            if (kappa == d - 1
                && d_mat_entry(r, kappa, kappa) >= d_mat_entry(r, kappa - 1,
                                                               kappa -
                                                               1) / 2.0)
            {
                fmpz_init(rii);
                for (i = d - 1; (i >= 0) && (ok > 0); i--)
                {
                    /* rii is the squared G-S length of ith vector */
                    fmpz_set_d(rii, d_mat_entry(r, i, i));
                    if ((ok = fmpz_cmp(rii, gs_B)) > 0)
                    {
                        d--;
                    }
                }
                fmpz_clear(rii);
            }
            newd = d;
            if (kappa >= d)
            {
                break;
            }
#else
            if (!ok)
            {
                fmpz_init(rii);
                fmpz_clear(rii);
            }
#endif

            if (ctt * d_mat_entry(r, kappa - 1, kappa - 1) <= s[kappa - 1]) /* check LLL condition */
            {
                kappa++;
            }
            else
            {
                int kappa2;
                kappa2 = kappa;
                do
                {
                    kappa--;
                } while ((kappa > 0)
                         && (ctt * d_mat_entry(r, kappa - 1, kappa - 1) >
                             s[kappa - 1]));

                for (j = 0; j < kappa; j++)
                {
                    d_mat_entry(mu, kappa, j) = d_mat_entry(mu, kappa2, j);
                    d_mat_entry(r, kappa, j) = d_mat_entry(r, kappa2, j);
                }
                d_mat_entry(r, kappa, kappa) = s[kappa];

                x = _fmpz_vec_init(d);
                if (fl->store_trans)
                {
                    _fmpz_vec_set(x, U->rows[kappa2], d);
                    for (i = kappa2; i > kappa; i--)
                        _fmpz_vec_set(U->rows[i], U->rows[i - 1], d);
                    _fmpz_vec_set(U->rows[kappa], x, d);
                }
                _fmpz_vec_clear(x, d);

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

        for (i = 0; i < B->r - 1; i++)
        {
            for (j = i + 1; j < B->r; j++)
            {
                fmpz_set(fmpz_mat_entry(B, i, j), fmpz_mat_entry(B, j, i));
            }
        }

        d_mat_clear(mu);
        d_mat_clear(r);
        _d_vec_clear(s);
        fmpz_clear(t);
        fmpz_clear(dmax);
    }
    return newd;
}

#endif
