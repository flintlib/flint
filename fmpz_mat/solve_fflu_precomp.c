/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_mat.h"

#define XX(ii,jj) fmpz_mat_entry(X,(ii),(jj))
#define BB(ii,jj) fmpz_mat_entry(B,(ii),(jj))
#define LU(ii,jj) fmpz_mat_entry(FFLU,(ii),(jj))

void
fmpz_mat_set_perm(fmpz_mat_t X, const slong * perm, const fmpz_mat_t B)
{
    if (X == B)
    {
        /* Not implemented */
        flint_abort();
    }
    else
    {
        slong i, j;

        if (perm == NULL)
            flint_abort();

        for (i = 0; i < fmpz_mat_nrows(B); i++)
            for (j = 0; j < fmpz_mat_ncols(B); j++)
                fmpz_set(fmpz_mat_entry(X, i, j),
                         fmpz_mat_entry(B, perm[i], j));
    }
}

void
fmpz_mat_solve_fflu_precomp(fmpz_mat_t X,
                    const slong * perm,
                    const fmpz_mat_t FFLU, const fmpz_mat_t B)
{
    fmpz_t T;
    slong i, j, k, m, n, norm = 0;
    ulong p1h, p1m, p1l, p2h, p2m, p2l, uden = 0, dinv = 0, quo;
    ulong FLINT_SET_BUT_UNUSED(rem);
    slong fbits = fmpz_mat_max_bits(FFLU);
    slong bbits = fmpz_mat_max_bits(B);
    int small = FLINT_ABS(fbits) <= FLINT_BITS - 2
             && FLINT_ABS(bbits) <= FLINT_BITS - 2;
    int sgn, dsgn = 0, den1 = 0, work_to_do;

    n = X->r;
    m = X->c;

    fmpz_init(T);
    fmpz_mat_set_perm(X, perm, B);

    for (k = 0; k < m; k++)
    {
        /* Fraction-free forward substitution */
        for (i = 0; i < n - 1; i++)
        {
            if (small)
            {
                den1 = 0;

                if (i > 0)
                {
                    uden = FLINT_ABS((slong)(*LU(i-1, i-1)));
                    dsgn = 0 > (slong)(*LU(i-1, i-1));
                    count_leading_zeros(norm, uden);
                    invert_limb(dinv, uden << norm);
                    den1 = fmpz_is_one(LU(i-1, i-1));
                }

                work_to_do = (!den1 || !fmpz_is_zero(XX(i, k)) ||
                      !fmpz_is_one(LU(i, i)));

                if (work_to_do)
                {
                    for (j = i + 1; j < n; j++)
                    {
                        smul_ppmm(p1h, p1l, *XX(j, k), *LU(i, i));
                        smul_ppmm(p2h, p2l, *LU(j, i), *XX(i, k));
                        sub_ddmmss(p1h, p1l, p1h, p1l, p2h, p2l);

                        sgn = 0 > (slong) p1h;

                        if (sgn) /* take absolute value */
                           sub_ddmmss(p1h, p1l, UWORD(0), UWORD(0), p1h, p1l);

                        if (i > 0 && !den1)
                        {
                            if (p1h >= uden)
                            {
                                fmpz_set_uiui(XX(j, k), p1h, p1l);

                                if (sgn)
                                    fmpz_neg(XX(j, k), XX(j, k));

                                fmpz_divexact(XX(j, k), XX(j, k), LU(i-1, i-1));

                                small = 0;
                            } else
                            {
                                udiv_qrnnd_preinv(quo, rem,
                                  (p1h << norm) + r_shift(p1l, (FLINT_BITS - norm)),
                                      p1l << norm, uden << norm, dinv);

                                if (sgn ^ dsgn)
                                    fmpz_neg_ui(XX(j, k), quo);
                                else
                                    fmpz_set_ui(XX(j, k), quo);

                                if (quo > COEFF_MAX)
                                    small = 0;
                            }
                        } else
                        {
                            if (p1h > 0)
                            {
                                fmpz_set_uiui(XX(j, k), p1h, p1l);

                                small = 0;
                            } else
                            {
                                fmpz_set_ui(XX(j, k), p1l);

                                if (p1l > COEFF_MAX)
                                    small = 0;
                            }

                            if (sgn)
                                fmpz_neg(XX(j, k), XX(j, k));
                        }
                    }
                }
            } else
            {
                for (j = i + 1; j < n; j++)
                {
                    fmpz_mul(XX(j, k), XX(j, k), LU(i, i));
                    fmpz_mul(T, LU(j, i), XX(i, k));
                    fmpz_sub(XX(j, k), XX(j, k), T);
                    if (i > 0)
                        fmpz_divexact(XX(j, k), XX(j, k), LU(i-1, i-1));
                }
            }
        }

        /* Fraction-free back substitution */
        for (i = n - 2; i >= 0; i--)
        {
            if (small)
            {
                smul_ppmm(p1m, p1l, *XX(i, k), *LU(n-1, n-1));
                p1h = FLINT_SIGN_EXT(p1m);

                uden = FLINT_ABS((slong)(*LU(i, i)));
                dsgn = 0 > (slong)(*LU(i, i));
                count_leading_zeros(norm, uden);
                invert_limb(dinv, uden << norm);

                for (j = i + 1; j < n; j++)
                {
                    if (!fmpz_is_zero(LU(i, j)))
                    {
                        smul_ppmm(p2m, p2l, *XX(j, k), *LU(i, j));
                        p2h = FLINT_SIGN_EXT(p2m);
                        sub_dddmmmsss(p1h, p1m, p1l, p1h, p1m, p1l, p2h, p2m, p2l);
                    }
                }

                sgn = 0 > (slong) p1h;

                if (sgn) /* take absolute value */
                   sub_dddmmmsss(p1h, p1m, p1l, UWORD(0), UWORD(0), UWORD(0), p1h, p1m, p1l);

                if (p1h != 0 || p1m >= uden)
                {
                    fmpz_set_signed_uiuiui(XX(i, k), p1h, p1m, p1l);

                    fmpz_divexact(XX(i, k), XX(i, k), LU(i, i));

                    if (sgn)
                        fmpz_neg(XX(i, k), XX(i, k));

                    small = 0;
                } else
                {
                    udiv_qrnnd_preinv(quo, rem,
                      (p1m << norm) + r_shift(p1l, (FLINT_BITS - norm)),
                          p1l << norm, uden << norm, dinv);

                    if (sgn ^ dsgn)
                        fmpz_neg_ui(XX(i, k), quo);
                    else
                        fmpz_set_ui(XX(i, k), quo);

                    if (quo > COEFF_MAX)
                        small = 0;
                }
            } else
            {
                fmpz_mul(XX(i, k), XX(i, k), LU(n-1, n-1));
                for (j = i + 1; j < n; j++)
                {
                    fmpz_mul(T, XX(j, k), LU(i, j));
                    fmpz_sub(XX(i, k), XX(i, k), T);
                }
                fmpz_divexact(XX(i, k), XX(i, k), LU(i, i));
            }
        }
    }

    fmpz_clear(T);
}
