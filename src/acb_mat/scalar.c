/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb_mat.h"

void acb_mat_scalar_mul_2exp_si(acb_mat_t B, const acb_mat_t A, slong c)
{
    slong i, j;

    for (i = 0; i < acb_mat_nrows(A); i++)
        for (j = 0; j < acb_mat_ncols(A); j++)
            acb_mul_2exp_si(acb_mat_entry(B, i, j), acb_mat_entry(A, i, j), c);
}

#define MAT_SCALAR_OP(func_name, scalar_type, scalar_op) \
void func_name(acb_mat_t B, const acb_mat_t A, scalar_type c, slong prec) \
{ \
    slong i, j; \
    for (i = 0; i < acb_mat_nrows(A); i++) \
        for (j = 0; j < acb_mat_ncols(A); j++) \
            scalar_op(acb_mat_entry(B, i, j), acb_mat_entry(A, i, j), c, prec); \
}

MAT_SCALAR_OP(acb_mat_scalar_addmul_si, slong, acb_addmul_si)
MAT_SCALAR_OP(   acb_mat_scalar_mul_si, slong,    acb_mul_si)
MAT_SCALAR_OP(   acb_mat_scalar_div_si, slong,    acb_div_si)

MAT_SCALAR_OP(acb_mat_scalar_addmul_fmpz, const fmpz_t, acb_addmul_fmpz)
MAT_SCALAR_OP(   acb_mat_scalar_mul_fmpz, const fmpz_t,    acb_mul_fmpz)
MAT_SCALAR_OP(   acb_mat_scalar_div_fmpz, const fmpz_t,    acb_div_fmpz)

MAT_SCALAR_OP(acb_mat_scalar_addmul_arb, const arb_t, acb_addmul_arb)
MAT_SCALAR_OP(   acb_mat_scalar_mul_arb, const arb_t,    acb_mul_arb)
MAT_SCALAR_OP(   acb_mat_scalar_div_arb, const arb_t,    acb_div_arb)

MAT_SCALAR_OP(acb_mat_scalar_addmul_acb, const acb_t, acb_addmul)
MAT_SCALAR_OP(   acb_mat_scalar_mul_acb, const acb_t,    acb_mul)
MAT_SCALAR_OP(   acb_mat_scalar_div_acb, const acb_t,    acb_div)

#undef MAT_SCALAR_OP
