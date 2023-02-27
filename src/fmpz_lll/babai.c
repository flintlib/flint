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

#if defined(FUNC_HEAD) && defined(LIMIT) && defined(COMPUTE) && defined(TYPE)
#ifdef GM
#undef GM
#endif
#define GM ((fl->rt == Z_BASIS) ? A->exactSP : B)

FUNC_HEAD
{
    if (fl->rt == Z_BASIS && fl->gt == APPROX)
    {
        int i, j, k, test, aa, exponent, max_expo = INT_MAX;
        slong xx;
        double tmp, rtmp, halfplus, onedothalfplus;
        ulong loops;

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
                                                                    zeros + 1);
                    rtmp = d_mat_entry(A->appSP, kappa, j) - tmp;

                    for (k = zeros + 2; k < j - 1; k++)
                    {
                        tmp = d_mat_entry(mu, j, k) * d_mat_entry(r, kappa, k);
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
                                                                    zeros + 1);
                    d_mat_entry(r, kappa, j) =
                        d_mat_entry(A->appSP, kappa, j) - tmp;
                }
                else
                    d_mat_entry(r, kappa, j) = d_mat_entry(A->appSP, kappa, j);

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
                        FLINT_MAX(new_max_expo, expo[kappa] - expo[j] + expo2);
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
                                tmp = ldexp(d_mat_entry(mu, j, k), exponent);
                                d_mat_entry(mu, kappa, k) =
                                    d_mat_entry(mu, kappa, k) - tmp;
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
                                tmp = ldexp(d_mat_entry(mu, j, k), exponent);
                                d_mat_entry(mu, kappa, k) =
                                    d_mat_entry(mu, kappa, k) + tmp;
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
                        tmp = ldexp(d_mat_entry(mu, kappa, j), -exponent);
                        if ((tmp < (double) FMPZ_LLL_MAX_LONG)
                            && (tmp > (double) -FMPZ_LLL_MAX_LONG))
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
                            if (U != NULL)
                            {
                                _fmpz_vec_scalar_submul_si(U->rows[kappa],
                                                           U->rows[j],
                                                           U->c, xx);
                            }
                        }
                        else
                        {
                            tmp = frexp(d_mat_entry(mu, kappa, j), &exponent);

                            tmp = tmp * FMPZ_LLL_MAX_LONG;
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
                                                           B->rows[j], n, xx);
                                if (U != NULL)
                                {
                                    _fmpz_vec_scalar_submul_si(U->rows
                                                               [kappa],
                                                               U->rows[j],
                                                               U->c, xx);
                                }

                                for (k = zeros + 1; k < j; k++)
                                {
                                    rtmp =
                                        ((double) xx) * d_mat_entry(mu, j, k);
                                    rtmp = ldexp(rtmp, expo[j] - expo[kappa]);
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
                                if (U != NULL)
                                {
                                    _fmpz_vec_scalar_submul_si_2exp(U->rows
                                                                    [kappa],
                                                                    U->rows
                                                                    [j],
                                                                    U->c,
                                                                    xx,
                                                                    exponent);
                                }

                                for (k = zeros + 1; k < j; k++)
                                {
                                    rtmp =
                                        ((double) xx) * d_mat_entry(mu, j, k);
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

            if (test)           /* Anything happened? */
            {
                expo[kappa] =
                    _fmpz_vec_get_d_vec_2exp(appB->rows[kappa],
                                             B->rows[kappa], n);
                aa = zeros + 1;

                for (i = zeros + 1; i <= LIMIT; i++)
                    d_mat_entry(A->appSP, kappa, i) = D_NAN;

                for (i = LIMIT + 1; i <= kappamax; i++)
                    d_mat_entry(A->appSP, i, kappa) = D_NAN;
            }
            else
            {
#if TYPE == 2
                for (i = zeros + 1; i <= LIMIT; i++)
                    d_mat_entry(A->appSP, kappa, i) = D_NAN;
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

        for (k = zeros + 1; k < kappa - 1; k++)
        {
            tmp = d_mat_entry(mu, kappa, k) * d_mat_entry(r, kappa, k);
            s[k + 1] = s[k] - tmp;
        }
#endif
    }
    else
    {
        int i, j, k, test, aa, exponent, max_expo = INT_MAX;
        slong exp;
        slong xx;
        double tmp, rtmp, halfplus, onedothalfplus;
        fmpz_t t;
        ulong loops;

        aa = (a > zeros) ? a : zeros + 1;

        fmpz_init(t);

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
                    tmp =
                        ldexp(d_mat_entry(mu, j, zeros + 1) * d_mat_entry(r,
                                                                          kappa,
                                                                          zeros
                                                                          + 1),
                              (expo[j] - expo[zeros + 1]));
                    rtmp = fmpz_get_d_2exp(&exp, fmpz_mat_entry(GM, kappa, j));
                    rtmp = ldexp(rtmp, (exp - expo[kappa])) - tmp;

                    for (k = zeros + 2; k < j - 1; k++)
                    {
                        tmp =
                            ldexp(d_mat_entry(mu, j, k) *
                                  d_mat_entry(r, kappa, k),
                                  (expo[j] - expo[k]));
                        rtmp = rtmp - tmp;
                    }

                    tmp =
                        ldexp(d_mat_entry(mu, j, j - 1) * d_mat_entry(r, kappa,
                                                                      j - 1),
                              (expo[j] - expo[j - 1]));
                    d_mat_entry(r, kappa, j) = rtmp - tmp;
                }
                else if (j == zeros + 2)
                {
                    tmp =
                        ldexp(d_mat_entry(mu, j, zeros + 1) * d_mat_entry(r,
                                                                          kappa,
                                                                          zeros
                                                                          + 1),
                              (expo[j] - expo[zeros + 1]));
                    d_mat_entry(r, kappa, j) =
                        fmpz_get_d_2exp(&exp, fmpz_mat_entry(GM, kappa, j));
                    d_mat_entry(r, kappa, j) =
                        ldexp(d_mat_entry(r, kappa, j),
                              (exp - expo[kappa])) - tmp;
                }
                else
                {
                    d_mat_entry(r, kappa, j) =
                        fmpz_get_d_2exp(&exp, fmpz_mat_entry(GM, kappa, j));
                    d_mat_entry(r, kappa, j) =
                        ldexp(d_mat_entry(r, kappa, j), (exp - expo[kappa]));
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
                        FLINT_MAX(new_max_expo, expo[kappa] - expo[j] + expo2);
                }
                if (new_max_expo > max_expo - SIZE_RED_FAILURE_THRESH)
                {
                    fmpz_clear(t);
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
                            fmpz_set_ui(x + j - zeros - 1, 1);
                            for (k = zeros + 1; k < j; k++)
                            {
                                tmp = ldexp(d_mat_entry(mu, j, k), exponent);
                                d_mat_entry(mu, kappa, k) =
                                    d_mat_entry(mu, kappa, k) - tmp;
                            }
                            if (fl->rt == Z_BASIS && B != NULL)
                            {
                                _fmpz_vec_sub(B->rows[kappa],
                                              B->rows[kappa], B->rows[j], n);
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
                            fmpz_set_si(x + j - zeros - 1, -WORD(1));
                            for (k = zeros + 1; k < j; k++)
                            {
                                tmp = ldexp(d_mat_entry(mu, j, k), exponent);
                                d_mat_entry(mu, kappa, k) =
                                    d_mat_entry(mu, kappa, k) + tmp;
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
                        tmp = ldexp(d_mat_entry(mu, kappa, j), -exponent);
                        if ((tmp < (double) FMPZ_LLL_MAX_LONG)
                            && (tmp > (double) -FMPZ_LLL_MAX_LONG))
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
                            fmpz_set_si(x + j - zeros - 1, xx);
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
                            tmp = frexp(d_mat_entry(mu, kappa, j), &exponent);

                            tmp = tmp * FMPZ_LLL_MAX_LONG;
                            xx = (slong) tmp;
                            exponent += expo[kappa] - expo[j] - CPU_SIZE_1;

                            /* This case is extremely rare: never happened for me. Check this: done */
                            if (exponent <= 0)
                            {
                                /* flint_printf("rare case kappa = %d, j = %d ******\n",
                                   kappa, j); */
                                xx = xx << -exponent;
                                exponent = 0;

                                fmpz_set_si(x + j - zeros - 1, xx);
                                if (fl->rt == Z_BASIS && B != NULL)
                                {
                                    _fmpz_vec_scalar_submul_si(B->rows
                                                               [kappa],
                                                               B->rows[j],
                                                               n, xx);
                                }
                                if (U != NULL)
                                {
                                    _fmpz_vec_scalar_submul_si(U->rows
                                                               [kappa],
                                                               U->rows[j],
                                                               U->c, xx);
                                }

                                for (k = zeros + 1; k < j; k++)
                                {
                                    rtmp =
                                        ((double) xx) * d_mat_entry(mu, j, k);
                                    rtmp = ldexp(rtmp, expo[j] - expo[kappa]);
                                    d_mat_entry(mu, kappa, k) =
                                        d_mat_entry(mu, kappa, k) - rtmp;
                                }
                            }
                            else
                            {
                                fmpz_set_si(x + j - zeros - 1, xx);
                                fmpz_mul_2exp(x + j - zeros - 1, x + j - zeros - 1, exponent);
                                if (fl->rt == Z_BASIS && B != NULL)
                                {
                                    _fmpz_vec_scalar_submul_si_2exp(B->rows
                                                                    [kappa],
                                                                    B->rows
                                                                    [j], n,
                                                                    xx,
                                                                    exponent);
                                }
                                if (U != NULL)
                                {
                                    _fmpz_vec_scalar_submul_si_2exp(U->rows
                                                                    [kappa],
                                                                    U->rows
                                                                    [j],
                                                                    U->c,
                                                                    xx,
                                                                    exponent);
                                }

                                for (k = zeros + 1; k < j; k++)
                                {
                                    rtmp =
                                        ((double) xx) * d_mat_entry(mu, j, k);
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

            if (test)           /* Anything happened? */
            {
                aa = zeros + 1;

                for (j = zeros + 1; j < kappa; j++)
                {
                    fmpz_pow_ui(t, x + j - zeros - 1, 2);
                    fmpz_addmul(fmpz_mat_entry(GM, kappa, kappa),
                                t, fmpz_mat_entry(GM, j, j));

                    fmpz_mul(t, x + j - zeros - 1, fmpz_mat_entry(GM, kappa, j));
                    fmpz_mul_2exp(t, t, 1);
                    fmpz_sub(fmpz_mat_entry(GM, kappa, kappa),
                             fmpz_mat_entry(GM, kappa, kappa), t);

                    for (i = zeros + 1; i < j; i++)
                    {
                        fmpz_mul(t, x + i - zeros - 1, x + j - zeros - 1);
                        fmpz_mul(t, t, fmpz_mat_entry(GM, j, i));
                        fmpz_mul_2exp(t, t, 1);
                        fmpz_add(fmpz_mat_entry(GM, kappa, kappa),
                                 fmpz_mat_entry(GM, kappa, kappa), t);
                    }
                }

                fmpz_get_d_2exp(&exp, fmpz_mat_entry(GM, kappa, kappa));
                expo[kappa] = exp;

                for (i = zeros + 1; i < kappa; i++)
                {
                    for (j = zeros + 1; j <= i; j++)
                        fmpz_submul(fmpz_mat_entry(GM, kappa, i),
                                    x + j - zeros - 1, fmpz_mat_entry(GM, i, j));
                    for (j = i + 1; j < kappa; j++)
                        fmpz_submul(fmpz_mat_entry(GM, kappa, i),
                                    x + j - zeros - 1, fmpz_mat_entry(GM, j, i));
                }

                for (i = kappa + 1; i < GM->r; i++)
                {
                    for (j = zeros + 1; j < kappa; j++)
                        fmpz_submul(fmpz_mat_entry(GM, i, kappa),
                                    x + j - zeros - 1, fmpz_mat_entry(GM, i, j));
                }
            }

            _fmpz_vec_clear(x, kappa - 1 - zeros);
            loops++;
        } while (test);

        s[zeros + 1] = fmpz_get_d_2exp(&exp, fmpz_mat_entry(GM, kappa, kappa));
        s[zeros + 1] = ldexp(s[zeros + 1], exp - expo[kappa]);

        for (k = zeros + 1; k < kappa - 1; k++)
        {
            tmp =
                ldexp(d_mat_entry(mu, kappa, k) * d_mat_entry(r, kappa, k),
                      (expo[kappa] - expo[k]));
            s[k + 1] = s[k] - tmp;
        }

        fmpz_clear(t);
    }
    return 0;
}

#undef GM

#else

void osxdummy18464823876() /* OSX doesn't like empty files */
{
}

#endif
