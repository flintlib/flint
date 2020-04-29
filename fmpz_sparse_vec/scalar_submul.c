/*
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by th e Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz_sparse_vec.h"

void fmpz_sparse_vec_scalar_submul_fmpz(fmpz_sparse_vec_t w, const fmpz_sparse_vec_t u, const fmpz_sparse_vec_t v, const fmpz_t c)
{
    slong unnz = u->nnz, vnnz = v->nnz, wnnz, k;
    fmpz_t tmp, tmp2;
    fmpz_sparse_entry_struct *ue, *ve, *we;

/*     flint_printf("c = "), fmpz_print(c); flint_printf("\n"); */
    /* Check for simpler operations first */
    if (fmpz_is_zero(c) || vnnz == 0) fmpz_sparse_vec_set(w, u, 0);
    else if (fmpz_is_one(c)) fmpz_sparse_vec_sub(w, u, v);
    else if (fmpz_equal_si(c, WORD(-1))) fmpz_sparse_vec_add(w, u, v);
    else if (unnz == 0) 
    {
        fmpz_sparse_vec_scalar_mul_fmpz(w, v, c);
        fmpz_sparse_vec_neg(w, w);
    }
    else
    {
        fmpz_init(tmp);
        fmpz_init_set(tmp2, c);
        wnnz = _fmpz_sparse_vec_union_nnz (u, v);
        _fmpz_sparse_vec_resize(w, wnnz);
        ue = u->entries + unnz, ve = v->entries + vnnz, we = w->entries + wnnz;
        /* flint_printf("c = "), fmpz_print(c); flint_printf(": %wd\n", c[0]); */
        while ((k = _fmpz_sparse_vector_merge_descend (&we, &ue, &ve, u, v)) >= 0)
        {
/*             flint_printf("\t\t\tc = "), fmpz_print(c); flint_printf("\n");
            flint_printf("\t\t\ttmp = %wd, tmp2 = %wd\n", tmp[0], tmp2[0]); */
            switch(k)
            {
                case 0: fmpz_set(we->val, ue->val); break;
                case 1: fmpz_mul(we->val, ve->val, tmp2); fmpz_neg(we->val, we->val); break;
                default: 
/*                 flint_printf("\t\tSubtracting: ");
                fmpz_print(we->val); flint_printf(" <- ");
                fmpz_print(ue->val); flint_printf("-"); 
                fmpz_print(ve->val); flint_printf("*"); fmpz_print(c); flint_printf("\n");
                flint_printf("%wd <- %wd - %wd * %wd\n", we->val[0], ue->val[0], ve->val[0], c[0]);
                flint_printf("c1: "); fmpz_print(c); flint_printf("\n"); */
                fmpz_mul(tmp, ve->val, tmp2);
/*                 flint_printf("%wd = %wd * %wd\n", tmp[0], ve->val[0], c[0]);
                flint_printf("c2: "); fmpz_print(c); flint_printf("\n"); */
                fmpz_sub(we->val, ue->val, tmp);
/*                 flint_printf("%wd = %wd - %wd, %wd\n", we->val[0], ue->val[0], tmp[0], c[0]);
                flint_printf("c3: "); fmpz_print(c); flint_printf("\n"); */
                if (fmpz_is_zero(we->val)) we++;
/*                 flint_printf("c4: "); fmpz_print(c); flint_printf("\n"); */
            }
/*             flint_printf("\t\tAfter subtraction: %wd, %wd\n", k, we->ind);
            flint_printf("\t\t\tr: "); fmpz_sparse_vec_print_pretty(w, 0, 0); flint_printf("\n");      */
        }
        _fmpz_sparse_vector_shift_left (w, we - w->entries);
        fmpz_clear(tmp);
        fmpz_clear(tmp2);
    }
}
