/*
    Copyright (C) 2012 Fredrik Johansson
    Copyright (C) 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"

arb_ptr _arb_vec_entry_ptr(arb_ptr vec, slong i)
{
    return vec + i;
}

void _arb_vec_scalar_mul_fmpz(arb_ptr res, arb_srcptr vec, slong len, const fmpz_t c, slong prec)
{
    slong i;
    arf_t t;
    arf_init(t);
    arf_set_fmpz(t, c);
    for (i = 0; i < len; i++)
        arb_mul_arf(res + i, vec + i, t, prec);
    arf_clear(t);
}

void _arb_vec_scalar_mul_2exp_si(arb_ptr res, arb_srcptr src, slong len, slong c)
{
    slong i;
    for (i = 0; i < len; i++)
        arb_mul_2exp_si(res + i, src + i, c);
}

slong _arb_vec_bits(arb_srcptr x, slong len)
{
    slong i, b, c;

    b = 0;
    for (i = 0; i < len; i++)
    {
        c = arb_bits(x + i);
        b = FLINT_MAX(b, c);
    }

    return b;
}

slong _arb_vec_allocated_bytes(arb_srcptr vec, slong len)
{
    slong i, size;

    size = len * sizeof(arb_struct);

    for (i = 0; i < len; i++)
        size += arb_allocated_bytes(vec + i);

    return size;
}

double _arb_vec_estimate_allocated_bytes(slong len, slong prec)
{
    double size;

    size = len * (double) sizeof(arb_struct);

    if (prec > ARF_NOPTR_LIMBS * FLINT_BITS)
        size += len * (double) ((prec + FLINT_BITS - 1) / FLINT_BITS) * sizeof(ulong);

    return size;
}

#define IS_OP(func_name, T, OP) \
int func_name(T ap, slong len)  \
{                               \
    slong ix;                   \
                                \
    for (ix = 0; ix < len; ix++)\
        if (!OP(ap + ix))       \
            return 0;           \
                                \
    return 1;                   \
}

IS_OP(_arb_vec_is_zero,   arb_srcptr, arb_is_zero)
IS_OP(_arb_vec_is_finite, arb_srcptr, arb_is_finite)

#define SET_OP_1(func_name, T, OP) \
void func_name(T ap, slong len) \
{                               \
    slong ix;                   \
                                \
    for (ix = 0; ix < len; ix++)\
        OP(ap + ix);            \
}

SET_OP_1(_arb_vec_zero,          arb_ptr, arb_zero)
SET_OP_1(_arb_vec_indeterminate, arb_ptr, arb_indeterminate)

#define SET_OP(func_name, S, T, OP) \
void func_name(S ap, T bp, slong len) \
{                               \
    slong ix;                   \
                                \
    for (ix = 0; ix < len; ix++)\
        OP(ap + ix, bp + ix);   \
}

SET_OP(_arb_vec_set,               arb_ptr, arb_srcptr, arb_set)
SET_OP(_arb_vec_swap,              arb_ptr, arb_ptr,    arb_swap)
SET_OP(_arb_vec_neg,               arb_ptr, arb_srcptr, arb_neg)
SET_OP(_arb_vec_trim,              arb_ptr, arb_srcptr, arb_trim)
SET_OP(_arb_vec_add_error_arf_vec, arb_ptr, arf_srcptr, arb_add_error_arf)
SET_OP(_arb_vec_add_error_mag_vec, arb_ptr, mag_srcptr, arb_add_error_mag)

#define SET_PREC_OP(func_name, S, T, OP) \
void func_name(S ap, T bp, slong len, slong prec) \
{                               \
    slong ix;                   \
                                \
    for (ix = 0; ix < len; ix++)\
        OP(ap + ix, bp + ix, prec); \
}

SET_PREC_OP(_arb_vec_set_round, arb_ptr, arb_srcptr, arb_set_round)

#define AORS_OP(func_name, S, T, OP) \
void func_name(S rp, T ap, T bp, slong len, slong prec) \
{                               \
    slong ix;                   \
                                \
    for (ix = 0; ix < len; ix++)\
        OP(rp + ix, ap + ix, bp + ix, prec); \
}

AORS_OP(_arb_vec_add, arb_ptr, arb_srcptr, arb_add)
AORS_OP(_arb_vec_sub, arb_ptr, arb_srcptr, arb_sub)

#define SCALAR_OP(func_name, S, T, U, OP) \
void func_name(S rp, T ap, slong len, U b, slong prec) \
{                               \
    slong ix;                   \
                                \
    for (ix = 0; ix < len; ix++)\
        OP(rp + ix, ap + ix, b, prec); \
}

SCALAR_OP(_arb_vec_scalar_mul,      arb_ptr, arb_srcptr, const arb_t,  arb_mul)
SCALAR_OP(_arb_vec_scalar_div,      arb_ptr, arb_srcptr, const arb_t,  arb_div)
SCALAR_OP(_arb_vec_scalar_addmul,   arb_ptr, arb_srcptr, const arb_t, arb_addmul)

#define COMPARISON_OP(func_name, S, T, OP) \
int func_name(S ap, T bp, slong len)\
{                                   \
    slong ix;                       \
                                    \
    for (ix = 0; ix < len; ix++)    \
        if (!OP(ap + ix, bp + ix))  \
            return 0;               \
                                    \
    return 1;                       \
}

COMPARISON_OP(_arb_vec_equal,               arb_srcptr, arb_srcptr, arb_equal)
COMPARISON_OP(_arb_vec_overlaps,            arb_srcptr, arb_srcptr, arb_overlaps)
COMPARISON_OP(_arb_vec_contains,            arb_srcptr, arb_srcptr, arb_contains)
COMPARISON_OP(_arb_vec_get_unique_fmpz_vec, fmpz *,     arb_srcptr, arb_get_unique_fmpz)
