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
#include "arb.h"
#include "arb_mat.h"

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

IS_OP(arb_mat_is_exact,  arb_mat_t, arb_is_exact)
IS_OP(arb_mat_is_finite, arb_mat_t, arb_is_finite)
IS_OP(arb_mat_is_zero,   arb_mat_t, arb_is_zero)

COMPARISON_OP(arb_mat_eq,                arb_mat_t, arb_mat_t,  arb_eq)
COMPARISON_OP(arb_mat_equal,             arb_mat_t, arb_mat_t,  arb_equal)
COMPARISON_OP(arb_mat_overlaps,          arb_mat_t, arb_mat_t,  arb_overlaps)
COMPARISON_OP(arb_mat_contains,          arb_mat_t, arb_mat_t,  arb_contains)
COMPARISON_OP(arb_mat_contains_fmpz_mat, arb_mat_t, fmpz_mat_t, arb_contains_fmpz)
COMPARISON_OP(arb_mat_contains_fmpq_mat, arb_mat_t, fmpq_mat_t, arb_contains_fmpq)

NOT_COMPARISON_OP(arb_mat_ne, arb_mat_t, arb_mat_t, arb_ne)

SET_OP_1(arb_mat_indeterminate, arb_mat_t, arb_indeterminate)
SET_OP_1(arb_mat_zero,          arb_mat_t, arb_zero)
SET_OP_1(arb_mat_ones,          arb_mat_t, arb_one)

SET_OP(arb_mat_swap_entrywise, arb_mat_t, arb_mat_t,        arb_swap)
SET_OP(arb_mat_get_mid,        arb_mat_t, const arb_mat_t,  arb_get_mid_arb)
SET_OP(arb_mat_neg,            arb_mat_t, const arb_mat_t,  arb_neg)
SET_OP(arb_mat_set_fmpz_mat,   arb_mat_t, const fmpz_mat_t, arb_set_fmpz)

SET_PREC_OP(arb_mat_set_round_fmpz_mat, arb_mat_t, fmpz_mat_t, arb_set_round_fmpz)
SET_PREC_OP(arb_mat_set_fmpq_mat, arb_mat_t, fmpq_mat_t, arb_set_fmpq)

AORS_NATIVE_OP(arb_mat_add, arb_mat_t, arb_add)
AORS_NATIVE_OP(arb_mat_sub, arb_mat_t, arb_sub)

SCALAR_OP(arb_mat_scalar_mul_si,   arb_mat_t, slong,        arb_mul_si)
SCALAR_OP(arb_mat_scalar_mul_fmpz, arb_mat_t, const fmpz_t, arb_mul_fmpz)
SCALAR_OP(arb_mat_scalar_mul_arb,  arb_mat_t, const arb_t,  arb_mul)

SCALAR_OP(arb_mat_scalar_addmul_si,   arb_mat_t, slong,        arb_addmul_si)
SCALAR_OP(arb_mat_scalar_addmul_fmpz, arb_mat_t, const fmpz_t, arb_addmul_fmpz)
SCALAR_OP(arb_mat_scalar_addmul_arb,  arb_mat_t, const arb_t,  arb_addmul)

SCALAR_OP(arb_mat_scalar_div_si,   arb_mat_t, slong,        arb_div_si)
SCALAR_OP(arb_mat_scalar_div_fmpz, arb_mat_t, const fmpz_t, arb_div_fmpz)
SCALAR_OP(arb_mat_scalar_div_arb,  arb_mat_t, const arb_t,  arb_div)
