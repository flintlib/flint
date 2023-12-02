/*
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpf_mat.h"
#include "gmpcompat.h"

#ifdef __GNUC__
# define ceil __builtin_ceil
#else
# include <math.h>
#endif

void
mpf_mat_gso(mpf_mat_t B, const mpf_mat_t A)
{
    slong i, j, k;
    int flag;
    mpf_t t, s, tmp, eps;
    flint_bitcnt_t exp;

    if (B->r != A->r || B->c != A->c)
    {
        flint_throw(FLINT_ERROR, "Exception (mpf_mat_gso). Incompatible dimensions.\n");
    }

    if (B == A)
    {
        mpf_mat_t T;
        mpf_mat_init(T, A->r, A->c, B->prec);
        mpf_mat_gso(T, A);
        mpf_mat_swap_entrywise(B, T);
        mpf_mat_clear(T);
        return;
    }

    if (A->r == 0)
    {
        return;
    }

    mpf_init2(t, B->prec);
    mpf_init2(s, B->prec);
    mpf_init2(tmp, B->prec);
    mpf_init2(eps, B->prec);
    exp = ceil((A->prec) / 64.0) * 64;
    flint_mpf_set_ui(eps, 1);
    mpf_div_2exp(eps, eps, exp);

    for (k = 0; k < A->c; k++)
    {
        for (j = 0; j < A->r; j++)
        {
            mpf_set(mpf_mat_entry(B, j, k), mpf_mat_entry(A, j, k));
        }
        flag = 1;
        while (flag)
        {
            flint_mpf_set_ui(t, 0);
            for (i = 0; i < k; i++)
            {
                flint_mpf_set_ui(s, 0);
                for (j = 0; j < A->r; j++)
                {
                    mpf_mul(tmp, mpf_mat_entry(B, j, i),
                            mpf_mat_entry(B, j, k));
                    mpf_add(s, s, tmp);
                }
                mpf_mul(tmp, s, s);
                mpf_add(t, t, tmp);
                for (j = 0; j < A->r; j++)
                {
                    mpf_mul(tmp, s, mpf_mat_entry(B, j, i));
                    mpf_sub(mpf_mat_entry(B, j, k), mpf_mat_entry(B, j, k),
                            tmp);
                }
            }
            flint_mpf_set_ui(s, 0);
            for (j = 0; j < A->r; j++)
            {
                mpf_mul(tmp, mpf_mat_entry(B, j, k), mpf_mat_entry(B, j, k));
                mpf_add(s, s, tmp);
            }
            mpf_add(t, t, s);
            flag = 0;
            if (mpf_cmp(s, t) < 0)
            {
                if (mpf_cmp(s, eps) < 0)
                    flint_mpf_set_ui(s, 0);
                else
                    flag = 1;
            }
        }
        mpf_sqrt(s, s);
        if (flint_mpf_cmp_ui(s, 0) != 0)
            mpf_ui_div(s, 1, s);
        for (j = 0; j < A->r; j++)
        {
            mpf_mul(mpf_mat_entry(B, j, k), mpf_mat_entry(B, j, k), s);
        }
    }
    mpf_clears(t, s, tmp, eps, NULL);
}
