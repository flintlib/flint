/*
    Copyright (C) 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mat.h"
#include "fmpq_mat.h"
#include "acb.h"
#include "acb_mat.h"

/* helper macros *************************************************************/

#define mat_entry(mat, i, j) ((mat)->entries + (i) * (mat)->stride + (j))
#define mat_nrows(mat) ((mat)->r)
#define mat_ncols(mat) ((mat)->c)

/* comparisons ***************************************************************/

/* Checks if matrix fulfills a criteria */
#define IS_OP(func_name, T, OP)                 \
int func_name(const T am)                       \
{                                               \
    slong ix, jx;                               \
                                                \
    for (ix = 0; ix < mat_nrows(am); ix++)      \
        for (jx = 0; jx < mat_ncols(am); jx++)  \
            if (!OP(mat_entry(am, ix, jx)))     \
                return 0;                       \
                                                \
    return 1;                                   \
}

/* Compares two matrices of possibly different types */
#define COMPARISON_OP(func_name, S, T, OP)      \
int func_name(const S am, const T bm)           \
{                                               \
    slong ix, jx;                               \
                                                \
    if (mat_nrows(am) != mat_nrows(bm) || mat_ncols(am) != mat_ncols(bm)) \
        return 0;                               \
                                                \
    for (ix = 0; ix < mat_nrows(am); ix++)      \
        for (jx = 0; jx < mat_ncols(am); jx++)  \
            if (!OP(mat_entry(am, ix, jx), mat_entry(bm, ix, jx))) \
                return 0;                       \
                                                \
    return 1;                                   \
}

#define NOT_COMPARISON_OP(func_name, S, T, OP)  \
int func_name(const S am, const T bm)           \
{                                               \
    slong ix, jx;                               \
                                                \
    if (mat_nrows(am) != mat_nrows(bm) || mat_ncols(am) != mat_nrows(bm)) \
        return 1;                               \
                                                \
    for (ix = 0; ix < mat_nrows(am); ix++)      \
        for (jx = 0; jx < mat_ncols(am); jx++)  \
            if (OP(mat_entry(am, ix, jx), mat_entry(bm, ix, jx))) \
                return 1;                       \
                                                \
    return 0;                                   \
}

/* setters and getters *******************************************************/

/* Sets matrix to one single thing */
#define SET_OP_1(func_name, T, OP)              \
void func_name(T rm)                            \
{                                               \
    slong ix, jx;                               \
                                                \
    for (ix = 0; ix < mat_nrows(rm); ix++)      \
        for (jx = 0; jx < mat_ncols(rm); jx++)  \
            OP(mat_entry(rm, ix, jx));          \
}

/* Set-swap operation of same type type */
#define SET_OP(func_name, S, T, OP)             \
void func_name(S am, T bm)                      \
{                                               \
    slong ix, jx;                               \
                                                \
    for (ix = 0; ix < mat_nrows(bm); ix++)      \
        for (jx = 0; jx < mat_ncols(bm); jx++)  \
            OP(mat_entry(am, ix, jx), mat_entry(bm, ix, jx)); \
}

/* Set operation with precision */
#define SET_PREC_OP(func_name, S, T, OP)        \
void func_name(S rm, const T am, slong prec)    \
{                                               \
    slong ix, jx;                               \
                                                \
    for (ix = 0; ix < mat_nrows(rm); ix++)      \
        for (jx = 0; jx < mat_ncols(rm); jx++)  \
            OP(mat_entry(rm, ix, jx), mat_entry(am, ix, jx), prec); \
}

/* linear arithmetic *********************************************************/

#define AORS_NATIVE_OP(func_name, T, OP)        \
void func_name(T rm, const T am, const T bm, slong prec) \
{                                               \
    slong ix, jx;                               \
                                                \
    for (ix = 0; ix < mat_nrows(am); ix++)      \
        for (jx = 0; jx < mat_ncols(am); jx++)  \
            OP(mat_entry(rm, ix, jx),           \
               mat_entry(am, ix, jx),           \
               mat_entry(bm, ix, jx), prec);    \
}

#define SCALAR_OP(func_name, S, T, OP)   \
void func_name(S rm, const S am, T bs, slong prec) \
{                                               \
    slong ix, jx;                               \
    for (ix = 0; ix < mat_nrows(am); ix++)      \
        for (jx = 0; jx < mat_ncols(am); jx++)  \
            OP(mat_entry(rm, ix, jx), mat_entry(am, ix, jx), bs, prec); \
}

/* define functions **********************************************************/

IS_OP(acb_mat_is_exact,  acb_mat_t, acb_is_exact)
IS_OP(acb_mat_is_finite, acb_mat_t, acb_is_finite)
IS_OP(acb_mat_is_zero,   acb_mat_t, acb_is_zero)
IS_OP(acb_mat_is_real,   acb_mat_t, acb_is_real)

COMPARISON_OP(acb_mat_eq,                acb_mat_t, acb_mat_t,  acb_eq)
COMPARISON_OP(acb_mat_equal,             acb_mat_t, acb_mat_t,  acb_equal)
COMPARISON_OP(acb_mat_overlaps,          acb_mat_t, acb_mat_t,  acb_overlaps)
COMPARISON_OP(acb_mat_contains,          acb_mat_t, acb_mat_t,  acb_contains)
COMPARISON_OP(acb_mat_contains_fmpz_mat, acb_mat_t, fmpz_mat_t, acb_contains_fmpz)
COMPARISON_OP(acb_mat_contains_fmpq_mat, acb_mat_t, fmpq_mat_t, acb_contains_fmpq)

NOT_COMPARISON_OP(acb_mat_ne, acb_mat_t, acb_mat_t, acb_ne)

SET_OP_1(acb_mat_indeterminate, acb_mat_t, acb_indeterminate)
SET_OP_1(acb_mat_zero,          acb_mat_t, acb_zero)
SET_OP_1(acb_mat_ones,          acb_mat_t, acb_one)

SET_OP(acb_mat_swap_entrywise, acb_mat_t, acb_mat_t,        acb_swap)
SET_OP(acb_mat_get_mid,        acb_mat_t, const acb_mat_t,  acb_get_mid)
SET_OP(acb_mat_neg,            acb_mat_t, const acb_mat_t,  acb_neg)
SET_OP(acb_mat_set_fmpz_mat,   acb_mat_t, const fmpz_mat_t, acb_set_fmpz)
SET_OP(acb_mat_set_arb_mat,    acb_mat_t, const arb_mat_t,  acb_set_arb)
SET_OP(acb_mat_conjugate,      acb_mat_t, const acb_mat_t,  acb_conj)

SET_PREC_OP(acb_mat_set_round_fmpz_mat, acb_mat_t, fmpz_mat_t, acb_set_round_fmpz)
SET_PREC_OP(acb_mat_set_fmpq_mat, acb_mat_t, fmpq_mat_t, acb_set_fmpq)
SET_PREC_OP(acb_mat_set_round_arb_mat,  acb_mat_t, arb_mat_t,  acb_set_round_arb)

AORS_NATIVE_OP(acb_mat_add, acb_mat_t, acb_add)
AORS_NATIVE_OP(acb_mat_sub, acb_mat_t, acb_sub)

SCALAR_OP(acb_mat_scalar_mul_si,   acb_mat_t, slong,        acb_mul_si)
SCALAR_OP(acb_mat_scalar_mul_fmpz, acb_mat_t, const fmpz_t, acb_mul_fmpz)
SCALAR_OP(acb_mat_scalar_mul_arb,  acb_mat_t, const arb_t,  acb_mul_arb)
SCALAR_OP(acb_mat_scalar_mul_acb,  acb_mat_t, const acb_t,  acb_mul)

SCALAR_OP(acb_mat_scalar_addmul_si,   acb_mat_t, slong,        acb_addmul_si)
SCALAR_OP(acb_mat_scalar_addmul_fmpz, acb_mat_t, const fmpz_t, acb_addmul_fmpz)
SCALAR_OP(acb_mat_scalar_addmul_arb,  acb_mat_t, const arb_t,  acb_addmul_arb)
SCALAR_OP(acb_mat_scalar_addmul_acb,  acb_mat_t, const acb_t,  acb_addmul)

SCALAR_OP(acb_mat_scalar_div_si,   acb_mat_t, slong,        acb_div_si)
SCALAR_OP(acb_mat_scalar_div_fmpz, acb_mat_t, const fmpz_t, acb_div_fmpz)
SCALAR_OP(acb_mat_scalar_div_arb,  acb_mat_t, const arb_t,  acb_div_arb)
SCALAR_OP(acb_mat_scalar_div_acb,  acb_mat_t, const acb_t,  acb_div)
