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

#if defined(FUNC_HEAD) && defined(LIMIT) && defined(COMPUTE) && defined(TYPE)

FUNC_HEAD
{
    if (fl->rt == Z_BASIS)
    {
        if (fl->gt == APPROX)
        {
            int i, j, k, test, aa, exponent, max_expo = INT_MAX;
            slong xx;
            double tmp, rtmp, halfplus, onedothalfplus;
            ulong loops;

            aa = (a > zeros) ? a : zeros + 1;

            halfplus = (4 * fl->eta + 0.5) / 5;
            onedothalfplus = 1.0 + halfplus;

            loops = 0;

            do
            {
                test = 0;

                /* ************************************** */
                /* Step2: compute the GSO for stage kappa */
                /* ************************************** */

                for (j = aa; j < LIMIT; j++)
                {
                    if (d_is_nan(d_mat_entry(A->appSP, kappa, j)))
                    {
                        COMPUTE(A->appSP, kappa, j, n);
                    }

                    if (j > zeros + 2)
                    {
                        tmp =
                            d_mat_entry(mu, j, zeros + 1) * d_mat_entry(r,
                                                                        kappa,
                                                                        zeros +
                                                                        1);
                        rtmp = d_mat_entry(A->appSP, kappa, j) - tmp;

                        for (k = zeros + 2; k < j - 1; k++)
                        {
                            tmp =
                                d_mat_entry(mu, j, k) * d_mat_entry(r, kappa,
                                                                    k);
                            rtmp = rtmp - tmp;
                        }

                        tmp =
                            d_mat_entry(mu, j, j - 1) * d_mat_entry(r, kappa,
                                                                    j - 1);
                        d_mat_entry(r, kappa, j) = rtmp - tmp;
                    }
                    else if (j == zeros + 2)
                    {
                        tmp =
                            d_mat_entry(mu, j, zeros + 1) * d_mat_entry(r,
                                                                        kappa,
                                                                        zeros +
                                                                        1);
                        d_mat_entry(r, kappa, j) =
                            d_mat_entry(A->appSP, kappa, j) - tmp;
                    }
                    else
                        d_mat_entry(r, kappa, j) =
                            d_mat_entry(A->appSP, kappa, j);

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
                        new_max_expo =
                            FLINT_MAX(new_max_expo,
                                      expo[kappa] - expo[j] + expo2);
                    }
                    if (new_max_expo > max_expo - SIZE_RED_FAILURE_THRESH)
                    {
                        return -1;
                    }
                    max_expo = new_max_expo;
                }

                /* **************************** */
                /* Step3--5: compute the X_j's  */
                /* **************************** */

                for (j = LIMIT - 1; j > zeros; j--)
                {
                    /* test of the relaxed size-reduction condition */
                    tmp = fabs(d_mat_entry(mu, kappa, j));
                    tmp = ldexp(tmp, expo[kappa] - expo[j]);

                    if (tmp > halfplus)
                    {
                        test = 1;
                        exponent = expo[j] - expo[kappa];

                        /* we consider separately the cases X = +-1 */
                        if (tmp <= onedothalfplus)
                        {
                            if (d_mat_entry(mu, kappa, j) >= 0) /* in this case, X is 1 */
                            {
                                for (k = zeros + 1; k < j; k++)
                                {
                                    tmp =
                                        ldexp(d_mat_entry(mu, j, k), exponent);
                                    d_mat_entry(mu, kappa, k) =
                                        d_mat_entry(mu, kappa, k) - tmp;
                                }
                                _fmpz_vec_sub(B->rows[kappa], B->rows[kappa],
                                              B->rows[j], n);
                                if (fl->store_trans)
                                {
                                    _fmpz_vec_sub(U->rows[kappa],
                                                  U->rows[kappa], U->rows[j],
                                                  B->r);
                                }
                            }
                            else    /* otherwise X is -1 */
                            {
                                for (k = zeros + 1; k < j; k++)
                                {
                                    tmp =
                                        ldexp(d_mat_entry(mu, j, k), exponent);
                                    d_mat_entry(mu, kappa, k) =
                                        d_mat_entry(mu, kappa, k) + tmp;
                                }
                                _fmpz_vec_add(B->rows[kappa], B->rows[kappa],
                                              B->rows[j], n);
                                if (fl->store_trans)
                                {
                                    _fmpz_vec_add(U->rows[kappa],
                                                  U->rows[kappa], U->rows[j],
                                                  B->r);
                                }
                            }
                        }
                        else    /* we must have |X| >= 2 */
                        {
                            tmp = ldexp(d_mat_entry(mu, kappa, j), -exponent);
                            if ((tmp < (double) MAX_LONG)
                                && (tmp > (double) -MAX_LONG))
                            {
                                if (tmp < 0)
                                    tmp = ceil(tmp - 0.5);
                                else
                                    tmp = floor(tmp + 0.5);

                                for (k = zeros + 1; k < j; k++)
                                {
                                    rtmp = tmp * d_mat_entry(mu, j, k);
                                    rtmp = ldexp(rtmp, exponent);
                                    d_mat_entry(mu, kappa, k) =
                                        d_mat_entry(mu, kappa, k) - rtmp;
                                }

                                xx = (slong) tmp;
                                _fmpz_vec_scalar_submul_si(B->rows[kappa],
                                                           B->rows[j], n, xx);
                                if (fl->store_trans)
                                {
                                    _fmpz_vec_scalar_submul_si(U->rows[kappa],
                                                               U->rows[j],
                                                               B->r, xx);
                                }
                            }
                            else
                            {
                                tmp =
                                    frexp(d_mat_entry(mu, kappa, j),
                                          &exponent);

                                tmp = tmp * MAX_LONG;
                                xx = (slong) tmp;
                                exponent += expo[kappa] - expo[j] - CPU_SIZE_1;

                                /* This case is extremely rare: never happened for me. Check this: done */
                                if (exponent <= 0)
                                {
                                    /* flint_printf("rare case kappa = %d, j = %d ******\n",
                                       kappa, j); */
                                    xx = xx << -exponent;
                                    exponent = 0;

                                    _fmpz_vec_scalar_submul_si(B->rows[kappa],
                                                               B->rows[j], n,
                                                               xx);
                                    if (fl->store_trans)
                                    {
                                        _fmpz_vec_scalar_submul_si(U->rows
                                                                   [kappa],
                                                                   U->rows[j],
                                                                   B->r, xx);
                                    }

                                    for (k = zeros + 1; k < j; k++)
                                    {
                                        rtmp =
                                            ((double) xx) * d_mat_entry(mu, j,
                                                                        k);
                                        rtmp =
                                            ldexp(rtmp, expo[j] - expo[kappa]);
                                        d_mat_entry(mu, kappa, k) =
                                            d_mat_entry(mu, kappa, k) - rtmp;
                                    }
                                }
                                else
                                {
                                    _fmpz_vec_scalar_submul_si_2exp(B->rows
                                                                    [kappa],
                                                                    B->rows[j],
                                                                    n, xx,
                                                                    exponent);
                                    if (fl->store_trans)
                                    {
                                        _fmpz_vec_scalar_submul_si_2exp(U->rows
                                                                        [kappa],
                                                                        U->rows
                                                                        [j],
                                                                        B->r,
                                                                        xx,
                                                                        exponent);
                                    }

                                    for (k = zeros + 1; k < j; k++)
                                    {
                                        rtmp =
                                            ((double) xx) * d_mat_entry(mu, j,
                                                                        k);
                                        rtmp =
                                            ldexp(rtmp,
                                                  exponent + expo[j] -
                                                  expo[kappa]);
                                        d_mat_entry(mu, kappa, k) =
                                            d_mat_entry(mu, kappa, k) - rtmp;
                                    }
                                }
                            }
                        }
                    }
                }

                if (test)       /* Anything happened? */
                {
                    expo[kappa] =
                        _fmpz_vec_get_d_vec_2exp(appB->rows[kappa],
                                                 B->rows[kappa], n);
                    aa = zeros + 1;

                    for (i = zeros + 1; i <= LIMIT; i++)
                        d_mat_entry(A->appSP, kappa, i) = NAN;

                    for (i = LIMIT + 1; i <= kappamax; i++)
                        d_mat_entry(A->appSP, i, kappa) = NAN;
                }
                else
                {
#if TYPE == 2
                    for (i = zeros + 1; i <= LIMIT; i++)
                        d_mat_entry(A->appSP, kappa, i) = NAN;
#endif
                }
                loops++;
            } while (test);

#if TYPE == 1
            if (d_is_nan(d_mat_entry(A->appSP, kappa, kappa)))
            {
                COMPUTE(A->appSP, kappa, kappa, n);
            }

            s[zeros + 1] = d_mat_entry(A->appSP, kappa, kappa);

            for (k = zeros + 1; k < kappa; k++)
            {
                tmp = d_mat_entry(mu, kappa, k) * d_mat_entry(r, kappa, k);
                s[k + 1] = s[k] - tmp;
            }
            d_mat_entry(r, kappa, kappa) = s[kappa];
#endif
        }
        else
        {
            int i, j, k, test, aa, exponent, max_expo = INT_MAX;
            slong exp;
            slong xx;
            double *appBtmp;
            double tmp, rtmp, halfplus, onedothalfplus;
            ulong loops;

            aa = (a > zeros) ? a : zeros + 1;

            appBtmp = _d_vec_init(n);

            halfplus = (4 * fl->eta + 0.5) / 5;
            onedothalfplus = 1.0 + halfplus;

            loops = 0;

            do
            {
                test = 0;

                /* ************************************** */
                /* Step2: compute the GSO for stage kappa */
                /* ************************************** */

                for (j = aa; j < LIMIT; j++)
                {
                    if (fmpz_cmp_si
                        (fmpz_mat_entry(A->exactSP, kappa, j), WORD_MIN) == 0)
                    {
                        _fmpz_vec_dot(fmpz_mat_entry(A->exactSP, kappa, j),
                                      B->rows[kappa], B->rows[j], n);
                    }

                    if (j > zeros + 2)
                    {
                        tmp =
                            d_mat_entry(mu, j, zeros + 1) * d_mat_entry(r,
                                                                        kappa,
                                                                        zeros +
                                                                        1);
                        rtmp =
                            fmpz_get_d_2exp(&exp,
                                            fmpz_mat_entry(A->exactSP, kappa,
                                                           j));
                        rtmp =
                            ldexp(rtmp, exp - (expo[kappa] + expo[j])) - tmp;

                        for (k = zeros + 2; k < j - 1; k++)
                        {
                            tmp =
                                d_mat_entry(mu, j, k) * d_mat_entry(r, kappa,
                                                                    k);
                            rtmp = rtmp - tmp;
                        }

                        tmp =
                            d_mat_entry(mu, j, j - 1) * d_mat_entry(r, kappa,
                                                                    j - 1);
                        d_mat_entry(r, kappa, j) = rtmp - tmp;
                    }
                    else if (j == zeros + 2)
                    {
                        tmp =
                            d_mat_entry(mu, j, zeros + 1) * d_mat_entry(r,
                                                                        kappa,
                                                                        zeros +
                                                                        1);
                        d_mat_entry(r, kappa, j) =
                            fmpz_get_d_2exp(&exp,
                                            fmpz_mat_entry(A->exactSP, kappa,
                                                           j));
                        d_mat_entry(r, kappa, j) =
                            ldexp(d_mat_entry(r, kappa, j),
                                  exp - (expo[kappa] + expo[j])) - tmp;
                    }
                    else
                    {
                        d_mat_entry(r, kappa, j) =
                            fmpz_get_d_2exp(&exp, fmpz_mat_entry
                                            (A->exactSP, kappa, j));
                        d_mat_entry(r, kappa, j) =
                            ldexp(d_mat_entry(r, kappa, j),
                                  exp - (expo[kappa] + expo[j]));
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
                        new_max_expo =
                            FLINT_MAX(new_max_expo,
                                      expo[kappa] - expo[j] + expo2);
                    }
                    if (new_max_expo > max_expo - SIZE_RED_FAILURE_THRESH)
                    {
                        _d_vec_clear(appBtmp);
                        return -1;
                    }
                    max_expo = new_max_expo;
                }

                /* **************************** */
                /* Step3--5: compute the X_j's  */
                /* **************************** */

                for (j = LIMIT - 1; j > zeros; j--)
                {
                    /* test of the relaxed size-reduction condition */
                    tmp = fabs(d_mat_entry(mu, kappa, j));
                    tmp = ldexp(tmp, expo[kappa] - expo[j]);

                    if (tmp > halfplus)
                    {
                        test = 1;
                        exponent = expo[j] - expo[kappa];

                        /* we consider separately the cases X = +-1 */
                        if (tmp <= onedothalfplus)
                        {
                            if (d_mat_entry(mu, kappa, j) >= 0) /* in this case, X is 1 */
                            {
                                for (k = zeros + 1; k < j; k++)
                                {
                                    tmp =
                                        ldexp(d_mat_entry(mu, j, k), exponent);
                                    d_mat_entry(mu, kappa, k) =
                                        d_mat_entry(mu, kappa, k) - tmp;
                                }
                                _fmpz_vec_sub(B->rows[kappa], B->rows[kappa],
                                              B->rows[j], n);
                                if (fl->store_trans)
                                {
                                    _fmpz_vec_sub(U->rows[kappa],
                                                  U->rows[kappa], U->rows[j],
                                                  B->r);
                                }
                            }
                            else    /* otherwise X is -1 */
                            {
                                for (k = zeros + 1; k < j; k++)
                                {
                                    tmp =
                                        ldexp(d_mat_entry(mu, j, k), exponent);
                                    d_mat_entry(mu, kappa, k) =
                                        d_mat_entry(mu, kappa, k) + tmp;
                                }
                                _fmpz_vec_add(B->rows[kappa], B->rows[kappa],
                                              B->rows[j], n);
                                if (fl->store_trans)
                                {
                                    _fmpz_vec_add(U->rows[kappa],
                                                  U->rows[kappa], U->rows[j],
                                                  B->r);
                                }
                            }
                        }
                        else    /* we must have |X| >= 2 */
                        {
                            tmp = ldexp(d_mat_entry(mu, kappa, j), -exponent);
                            if ((tmp < (double) MAX_LONG)
                                && (tmp > (double) -MAX_LONG))
                            {
                                if (tmp < 0)
                                    tmp = ceil(tmp - 0.5);
                                else
                                    tmp = floor(tmp + 0.5);

                                for (k = zeros + 1; k < j; k++)
                                {
                                    rtmp = tmp * d_mat_entry(mu, j, k);
                                    rtmp = ldexp(rtmp, exponent);
                                    d_mat_entry(mu, kappa, k) =
                                        d_mat_entry(mu, kappa, k) - rtmp;
                                }

                                xx = (slong) tmp;
                                _fmpz_vec_scalar_submul_si(B->rows[kappa],
                                                           B->rows[j], n, xx);
                                if (fl->store_trans)
                                {
                                    _fmpz_vec_scalar_submul_si(U->rows[kappa],
                                                               U->rows[j],
                                                               B->r, xx);
                                }
                            }
                            else
                            {
                                tmp =
                                    frexp(d_mat_entry(mu, kappa, j),
                                          &exponent);

                                tmp = tmp * MAX_LONG;
                                xx = (slong) tmp;
                                exponent += expo[kappa] - expo[j] - CPU_SIZE_1;

                                /* This case is extremely rare: never happened for me. Check this: done */
                                if (exponent <= 0)
                                {
                                    /* flint_printf("rare case kappa = %d, j = %d ******\n",
                                       kappa, j); */
                                    xx = xx << -exponent;
                                    exponent = 0;

                                    _fmpz_vec_scalar_submul_si(B->rows[kappa],
                                                               B->rows[j], n,
                                                               xx);
                                    if (fl->store_trans)
                                    {
                                        _fmpz_vec_scalar_submul_si(U->rows
                                                                   [kappa],
                                                                   U->rows[j],
                                                                   B->r, xx);
                                    }

                                    for (k = zeros + 1; k < j; k++)
                                    {
                                        rtmp =
                                            ((double) xx) * d_mat_entry(mu, j,
                                                                        k);
                                        rtmp =
                                            ldexp(rtmp, expo[j] - expo[kappa]);
                                        d_mat_entry(mu, kappa, k) =
                                            d_mat_entry(mu, kappa, k) - rtmp;
                                    }
                                }
                                else
                                {
                                    _fmpz_vec_scalar_submul_si_2exp(B->rows
                                                                    [kappa],
                                                                    B->rows[j],
                                                                    n, xx,
                                                                    exponent);
                                    if (fl->store_trans)
                                    {
                                        _fmpz_vec_scalar_submul_si_2exp(U->rows
                                                                        [kappa],
                                                                        U->rows
                                                                        [j],
                                                                        B->r,
                                                                        xx,
                                                                        exponent);
                                    }

                                    for (k = zeros + 1; k < j; k++)
                                    {
                                        rtmp =
                                            ((double) xx) * d_mat_entry(mu, j,
                                                                        k);
                                        rtmp =
                                            ldexp(rtmp,
                                                  exponent + expo[j] -
                                                  expo[kappa]);
                                        d_mat_entry(mu, kappa, k) =
                                            d_mat_entry(mu, kappa, k) - rtmp;
                                    }
                                }
                            }
                        }
                    }
                }

                if (test)       /* Anything happened? */
                {
                    expo[kappa] =
                        _fmpz_vec_get_d_vec_2exp(appBtmp, B->rows[kappa], n);
                    aa = zeros + 1;

                    for (i = zeros + 1; i <= LIMIT; i++)
                    {
                        fmpz_set_si(fmpz_mat_entry(A->exactSP, kappa, i),
                                    WORD_MIN);
                    }

                    for (i = LIMIT + 1; i <= kappamax; i++)
                    {
                        fmpz_set_si(fmpz_mat_entry(A->exactSP, i, kappa),
                                    WORD_MIN);
                    }
                }
                else
                {
#if TYPE == 2
                    for (i = zeros + 1; i <= LIMIT; i++)
                    {
                        fmpz_set_si(fmpz_mat_entry(A->exactSP, kappa, i),
                                    WORD_MIN);
                    }
#endif
                }
                loops++;
            } while (test);

#if TYPE == 1
            if (fmpz_cmp_si(fmpz_mat_entry(A->exactSP, kappa, kappa), WORD_MIN)
                == 0)
            {
                _fmpz_vec_dot(fmpz_mat_entry(A->exactSP, kappa, kappa),
                              B->rows[kappa], B->rows[kappa], n);
            }

            s[zeros + 1] =
                fmpz_get_d_2exp(&exp,
                                fmpz_mat_entry(A->exactSP, kappa, kappa));
            s[zeros + 1] = ldexp(s[zeros + 1], exp - 2 * expo[kappa]);

            for (k = zeros + 1; k < kappa; k++)
            {
                tmp = d_mat_entry(mu, kappa, k) * d_mat_entry(r, kappa, k);
                s[k + 1] = s[k] - tmp;
            }
            d_mat_entry(r, kappa, kappa) = s[kappa];
#endif
            _d_vec_clear(appBtmp);
        }
    }
    return 0;
}

#endif
