/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"
#include "arb_mat.h"

void arb_mat_scalar_mul_2exp_si(arb_mat_t B, const arb_mat_t A, slong c)
{
    slong i, j;

    for (i = 0; i < arb_mat_nrows(A); i++)
        for (j = 0; j < arb_mat_ncols(A); j++)
            arb_mul_2exp_si(arb_mat_entry(B, i, j), arb_mat_entry(A, i, j), c);
}

#define MAT_SCALAR_OP(func_name, scalar_type, scalar_op) \
void func_name(arb_mat_t B, const arb_mat_t A, scalar_type c, slong prec) \
{ \
    slong i, j; \
    for (i = 0; i < arb_mat_nrows(A); i++) \
        for (j = 0; j < arb_mat_ncols(A); j++) \
            scalar_op(arb_mat_entry(B, i, j), arb_mat_entry(A, i, j), c, prec); \
}

MAT_SCALAR_OP(arb_mat_scalar_addmul_si, slong, arb_addmul_si)
MAT_SCALAR_OP(   arb_mat_scalar_mul_si, slong,    arb_mul_si)
MAT_SCALAR_OP(   arb_mat_scalar_div_si, slong,    arb_div_si)

MAT_SCALAR_OP(arb_mat_scalar_addmul_fmpz, const fmpz_t, arb_addmul_fmpz)
MAT_SCALAR_OP(   arb_mat_scalar_mul_fmpz, const fmpz_t,    arb_mul_fmpz)
MAT_SCALAR_OP(   arb_mat_scalar_div_fmpz, const fmpz_t,    arb_div_fmpz)

MAT_SCALAR_OP(arb_mat_scalar_addmul_arb, const arb_t, arb_addmul)
MAT_SCALAR_OP(   arb_mat_scalar_mul_arb, const arb_t,    arb_mul)
MAT_SCALAR_OP(   arb_mat_scalar_div_arb, const arb_t,    arb_div)

#undef MAT_SCALAR_OP
