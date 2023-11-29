/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2020-2021 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mat.h"

#define XX(ii,jj) fmpz_mat_entry(X,(ii),(jj))
#define XXx(ii,jj) fmpz_mat_entry(Xx,(ii),(jj))
#define BB(ii,jj) fmpz_mat_entry(B,(ii),(jj))
#define LU(ii,jj) fmpz_mat_entry(FFLU,(ii),(jj))

void
fmpz_mat_set_perm(fmpz_mat_t X, const slong * perm, const fmpz_mat_t B)
{
    if (X == B)
    {
        /* Not implemented */
        flint_throw(FLINT_ERROR, "(%s): Not implemented\n", __func__);
    }
    else
    {
        slong i, j;

        if (perm == NULL)
            flint_throw(FLINT_ERROR, "(%s): perm == NULL\n", __func__);

        for (i = 0; i < fmpz_mat_nrows(B); i++)
            for (j = 0; j < fmpz_mat_ncols(B); j++)
                fmpz_set(fmpz_mat_entry(X, i, j),
                         fmpz_mat_entry(B, perm[i], j));
    }
}

int
fmpz_mat_solve_fflu_precomp(fmpz_mat_t X,
                    const slong * perm,
                    const fmpz_mat_t FFLU, const fmpz_mat_t B)
{
    fmpz_t T;
    slong c, i, j, k, l, m, n, rnk, norm = 0;
    ulong p1h, p1m, p1l, p2h, p2m, p2l, uden = 0, dinv = 0, quo;
    ulong FLINT_SET_BUT_UNUSED(rem);
    slong fbits = fmpz_mat_max_bits(FFLU);
    slong bbits = fmpz_mat_max_bits(B);
    int small = FLINT_ABS(fbits) <= SMALL_FMPZ_BITCOUNT_MAX
             && FLINT_ABS(bbits) <= SMALL_FMPZ_BITCOUNT_MAX;
    int sgn, dsgn = 0, den1 = 0, work_to_do, flag = 1;
    fmpz_mat_t Xx;
    fmpz * diag;
    slong * piv;

    n = B->r;
    m = B->c;
    c = FFLU->c;

    fmpz_init(T);
    fmpz_mat_init(Xx, B->r, B->c);

    fmpz_mat_set_perm(Xx, perm, B);

    diag = _fmpz_vec_init(n);
    piv = (slong *) flint_malloc(n*sizeof(slong));

    rnk = 0;
    for (i = 0, l = 0; i < n; i++, l++)
    {
        while (l < c && fmpz_is_zero(LU(i, l)))
            l++;

        piv[i] = l;

        if (l < c)
        {
            fmpz_set(diag + i, LU(i, l));
            rnk += 1;
        }
    }

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
                    uden = FLINT_ABS((slong)(diag[i - 1]));
                    dsgn = 0 > (slong)(diag[i - 1]);
                    if (uden != 0) /* see #1029 */
                    {
                       norm = flint_clz(uden);
                       dinv = n_preinvert_limb_prenorm(uden << norm);
                    } else
                       dinv = 0;
                    den1 = fmpz_is_one(diag + i - 1);
                }

                work_to_do = (!den1 || !fmpz_is_zero(XXx(i, k)) ||
                      !fmpz_is_one(diag + i));

                if (work_to_do)
                {
                    for (j = i + 1; j < n; j++)
                    {
                        if (i < c && piv[i] < c)
                        {
                            smul_ppmm(p1h, p1l, *XXx(j, k), diag[i]);
                            smul_ppmm(p2h, p2l, *LU(j, piv[i]), *XXx(i, k));
                            sub_ddmmss(p1h, p1l, p1h, p1l, p2h, p2l);

                            sgn = 0 > (slong) p1h;

                            if (sgn) /* take absolute value */
                                sub_ddmmss(p1h, p1l, UWORD(0), UWORD(0), p1h, p1l);

                            if (i > 0 && !den1 && i < c && piv[i - 1] < c)
                            {
                                if (p1h >= uden)
                                {
                                    fmpz_set_uiui(XXx(j, k), p1h, p1l);

                                    if (sgn)
                                        fmpz_neg(XXx(j, k), XXx(j, k));

                                    if (i < c && piv[i - 1] < c)
                                    {
                                        flag = fmpz_divides(XXx(j, k), XXx(j, k), diag + i - 1);
                                        if (!flag)
                                            goto cleanup;
                                    }

                                    small = 0;
                                } else
                                {
                                    udiv_qrnnd_preinv(quo, rem,
                                      (p1h << norm) + r_shift(p1l, (FLINT_BITS - norm)),
                                          p1l << norm, uden << norm, dinv);

                                    flag = rem == 0;
                                    if (!flag)
                                        goto cleanup;

                                    if (sgn ^ dsgn)
                                        fmpz_neg_ui(XXx(j, k), quo);
                                    else
                                        fmpz_set_ui(XXx(j, k), quo);

                                    if (quo > COEFF_MAX)
                                        small = 0;
                                }
                            } else
                            {
                                if (p1h > 0)
                                {
                                    fmpz_set_uiui(XXx(j, k), p1h, p1l);

                                    small = 0;
                                } else
                                {
                                    fmpz_set_ui(XXx(j, k), p1l);

                                    if (p1l > COEFF_MAX)
                                        small = 0;
                                }

                                if (sgn)
                                    fmpz_neg(XXx(j, k), XXx(j, k));
                            }
                        } else if (i > 0)
                        {
                            if (i < c && piv[i - 1] < c)
                            {
                                flag = fmpz_divides(XXx(j, k), XXx(j, k), diag + i - 1);
                                if (!flag)
                                   goto cleanup;
                            }
                        }
                    }
                }
            } else
            {
                for (j = i + 1; j < n; j++)
                {
                    if (i < c && piv[i] < c)
                    {
                        fmpz_mul(XXx(j, k), XXx(j, k), diag + i);
                        fmpz_mul(T, LU(j, piv[i]), XXx(i, k));
                        fmpz_sub(XXx(j, k), XXx(j, k), T);
                    }
                    if (i > 0)
                    {
                        if (i < c && piv[i - 1] < c)
                        {
                            flag = fmpz_divides(XXx(j, k), XXx(j, k), diag + i - 1);
                            if (!flag)
                               goto cleanup;
                        }
                    }
                }
            }
        }

        l = rnk - 1;
        /* Fraction-free back substitution */
        for (i = c - 1; i >= 0; i--)
        {
            if (l > -1 && i == piv[l])
            {
                if (small)
                {
                    if (rnk != 0)
                       smul_ppmm(p1m, p1l, *XXx(l, k), diag[rnk - 1]);
                        else
                    {
                       p1l = *XXx(l, k);
                       p1m = FLINT_SIGN_EXT(p1l);
                    }
                    p1h = FLINT_SIGN_EXT(p1m);

                    uden = FLINT_ABS((slong)(diag[l]));
                    dsgn = 0 > (slong)(diag[l]);
                    if (uden != 0) /* see #1029 */
                    {
                       norm = flint_clz(uden);
                       dinv = n_preinvert_limb_prenorm(uden << norm);
                    } else
                       dinv = 0;

                    for (j = piv[l] + 1; j < c; j++)
                    {
                        if (!fmpz_is_zero(LU(l, j)))
                        {
                            smul_ppmm(p2m, p2l, *XX(j, k), *LU(l, j));
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

                        flag = fmpz_divides(XX(i, k), XX(i, k), diag + l);
                        if (!flag)
                            goto cleanup;

                        if (sgn)
                            fmpz_neg(XX(i, k), XX(i, k));

                        small = 0;
                    } else
                    {
                        udiv_qrnnd_preinv(quo, rem,
                          (p1m << norm) + r_shift(p1l, (FLINT_BITS - norm)),
                              p1l << norm, uden << norm, dinv);

                        flag = rem == 0;
                        if (!flag)
                            goto cleanup;

                        if (sgn ^ dsgn)
                            fmpz_neg_ui(XX(i, k), quo);
                        else
                            fmpz_set_ui(XX(i, k), quo);

                        if (quo > COEFF_MAX)
                            small = 0;
                    }
                } else
                {
                    if (rnk != 0)
                        fmpz_mul(XX(i, k), XXx(l, k), diag + rnk - 1);
                    else
                        fmpz_set(XX(i, k), XXx(l, k));
                    for (j = piv[l] + 1; j < c; j++)
                    {
                        fmpz_mul(T, XX(j, k), LU(l, j));
                        fmpz_sub(XX(i, k), XX(i, k), T);
                    }
                    flag = fmpz_divides(XX(i, k), XX(i, k), diag + l);
                        if (!flag)
                        goto cleanup;
                }
                l--;
            } else
                 fmpz_zero(XX(i, k));
        }
    }

cleanup:

    _fmpz_vec_clear(diag, n);
    flint_free(piv);
    fmpz_mat_clear(Xx);
    fmpz_clear(T);

    return flag;
}
