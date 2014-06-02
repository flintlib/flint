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

    Copyright (C) 2005, 2006 Damien Stehlé
    Copyright (C) 2009, 2010 William Hart
    Copyright (C) 2009, 2010 Andy Novocin
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2014 Abhinav Baid

******************************************************************************/

#include "fmpz_lll.h"

int
fmpz_lll_d(fmpz_mat_t B, const fmpz_lll_t fl)
{
    if (fl->is_gram == 0)
    {
        if (fl->gram_type == 0)
        {
            int kappa, kappa2, d, n, i, j, zeros, kappamax, shift;
            int num_failed_fast = 0;
            int babai_ok = 0;
            int heuristic_fail = 0;
            d_mat_t mu, r, appB;
            fmpz_gram_t A;
            double *s, *mutmp, *appBtmp, *appSPtmp;
            double ctt;
            int *expo, *alpha;
            fmpz *Btmp;

            n = B->c;
            d = B->r;

            ctt = (4 * fl->delta + 1) / 5;

            shift = fmpz_lll_shift(B);

            alpha = (int *) malloc(d * sizeof(int));
            expo = (int *) malloc(d * sizeof(int));

            d_mat_init(mu, d, d);
            d_mat_init(r, d, d);
            d_mat_init(appB, d, n);
            d_mat_init(A->appSP, d, d);

            s = _d_vec_init(d);
            appSPtmp = _d_vec_init(d);

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

            for (i = 0; i < d; i++)
                expo[i] =
                    _fmpz_vec_get_d_vec_2exp(appB->rows[i], B->rows[i], n);

            /* ********************************* */
            /* Step2: Initializing the main loop */
            /* ********************************* */

            kappamax = 0;
            i = 0;

            do
                d_mat_entry(A->appSP, i, i) = _d_vec_norm(appB->rows[i], n);
            while ((d_mat_entry(A->appSP, i, i) <= 0.0) && (++i < d));  /* Check if this should be D_EPS and not 0.0: done */

            zeros = i - 1;      /* all vectors B[i] with i <= zeros are zero vectors */
            kappa = i + 1;
            kappamax = kappa;

            if (zeros < d - 1)
            {
                d_mat_entry(r, i, i) = d_mat_entry(A->appSP, i, i);
            }

            for (i = zeros + 1; i < d; i++)
                alpha[i] = 0;

            while (kappa < d)
            {
                double tmp = 0.0;

                if (kappa > kappamax)
                    kappamax = kappa;

                /* ********************************** */
                /* Step3: Call to the Babai algorithm */
                /* ********************************** */
                if (num_failed_fast < 50)
                {
                    babai_ok =
                        fmpz_lll_check_babai(kappa, B, mu, r, s, appB, expo, A,
                                             alpha[kappa], zeros, kappamax,
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
                        fmpz_lll_check_babai_heuristic_d(kappa, B, mu, r, s,
                                                         appB, expo, A,
                                                         alpha[kappa], zeros,
                                                         kappamax,
                                                         FLINT_MIN(kappamax +
                                                                   1 + shift,
                                                                   n), fl);
                }

                if (heuristic_fail == -1)
                {
                    free(alpha);
                    free(expo);
                    d_mat_clear(mu);
                    d_mat_clear(r);
                    d_mat_clear(appB);
                    d_mat_clear(A->appSP);
                    _d_vec_clear(s);
                    _d_vec_clear(appSPtmp);
                    /* Need to switch to mpf / arb */
                    return -1;
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
                                                                        kappa -
                                                                        1);
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

            free(alpha);
            free(expo);
            d_mat_clear(mu);
            d_mat_clear(r);
            d_mat_clear(appB);
            d_mat_clear(A->appSP);
            _d_vec_clear(s);
            _d_vec_clear(appSPtmp);
        }
        else
        {
            int kappa, kappa2, d, n, i, j, zeros, kappamax, shift;
            int num_failed_fast = 0;
            int babai_ok = 0;
            int heuristic_fail = 0;
            d_mat_t mu, r, appB;
            fmpz_gram_t A;
            double *s, *mutmp, *appBtmp;
            fmpz *exactSPtmp;
            double ctt;
            int *expo, *alpha;
            fmpz *Btmp;

            n = B->c;
            d = B->r;

            ctt = (4 * fl->delta + 1) / 5;

            shift = fmpz_lll_shift(B);

            alpha = (int *) malloc(d * sizeof(int));
            expo = (int *) malloc(d * sizeof(int));

            d_mat_init(mu, d, d);
            d_mat_init(r, d, d);
            d_mat_init(appB, d, n);
            fmpz_mat_init(A->exactSP, d, d);

            s = _d_vec_init(d);
            exactSPtmp = _fmpz_vec_init(d);

            fmpz_mat_gram(A->exactSP, B);

            /* ************************** */
            /* Step1: Initialization Step */
            /* ************************** */

            for (i = 0; i < d; i++)
                expo[i] =
                    _fmpz_vec_get_d_vec_2exp(appB->rows[i], B->rows[i], n);

            /* ********************************* */
            /* Step2: Initializing the main loop */
            /* ********************************* */

            kappamax = 0;
            i = 0;

            do
                ;
            while ((fmpz_cmp_ui(fmpz_mat_entry(A->exactSP, i, i), 0) <= 0)
                   && (++i < d));


            zeros = i - 1;      /* all vectors B[i] with i <= zeros are zero vectors */
            kappa = i + 1;
            kappamax = kappa;

            if (zeros < d - 1)
            {
                d_mat_entry(r, i, i) =
                    fmpz_get_d(fmpz_mat_entry(A->exactSP, i, i)) / ldexp(1,
                                                                         2 *
                                                                         expo
                                                                         [i]);
            }

            for (i = zeros + 1; i < d; i++)
                alpha[i] = 0;

            while (kappa < d)
            {
                double tmp = 0.0;

                if (kappa > kappamax)
                    kappamax = kappa;

                /* ********************************** */
                /* Step3: Call to the Babai algorithm */
                /* ********************************** */
                if (num_failed_fast < 50)
                {
                    babai_ok =
                        fmpz_lll_check_babai(kappa, B, mu, r, s, appB, expo, A,
                                             alpha[kappa], zeros, kappamax,
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
                        fmpz_lll_check_babai_heuristic_d(kappa, B, mu, r, s,
                                                         appB, expo, A,
                                                         alpha[kappa], zeros,
                                                         kappamax,
                                                         FLINT_MIN(kappamax +
                                                                   1 + shift,
                                                                   n), fl);
                }

                if (heuristic_fail == -1)
                {
                    free(alpha);
                    free(expo);
                    d_mat_clear(mu);
                    d_mat_clear(r);
                    d_mat_clear(appB);
                    fmpz_mat_clear(A->exactSP);
                    _d_vec_clear(s);
                    _fmpz_vec_clear(exactSPtmp, d);
                    /* Need to switch to mpf / arb */
                    return -1;
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
                                                                        kappa -
                                                                        1);
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
                            fmpz_get_d(fmpz_mat_entry
                                       (A->exactSP, kappa, kappa)) / ldexp(1,
                                                                           2 *
                                                                           expo
                                                                           [kappa]);
                    }

                    kappa++;
                }
            }

            free(alpha);
            free(expo);
            d_mat_clear(mu);
            d_mat_clear(r);
            d_mat_clear(appB);
            fmpz_mat_clear(A->exactSP);
            _d_vec_clear(s);
            _fmpz_vec_clear(exactSPtmp, d);
        }
    }
    return 0;
}
