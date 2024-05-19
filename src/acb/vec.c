/*
    Copyright (C) 2012 Fredrik Johansson
    Copyright (C) 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"

acb_ptr _acb_vec_entry_ptr(acb_ptr vec, slong i)
{
    return vec + i;
}

void _acb_vec_get_real(arb_ptr re, acb_srcptr vec, slong len)
{
    slong i;

    for (i = 0; i < len; i++)
        arb_set(re + i, acb_realref(vec + i));
}

void _acb_vec_get_imag(arb_ptr im, acb_srcptr vec, slong len)
{
    slong i;

    for (i = 0; i < len; i++)
        arb_set(im + i, acb_imagref(vec + i));
}

void _acb_vec_set_real_imag(acb_ptr vec, arb_srcptr re, arb_srcptr im, slong len)
{
    slong i;

    for (i = 0; i < len; i++)
        acb_set_arb_arb(vec + i, re + i, im + i);
}

slong _acb_vec_bits(acb_srcptr vec, slong len)
{
    return _arb_vec_bits((arb_srcptr) vec, 2 * len);
}

slong _acb_vec_allocated_bytes(acb_srcptr vec, slong len)
{
    return _arb_vec_allocated_bytes((arb_srcptr) vec, 2 * len);
}

double _acb_vec_estimate_allocated_bytes(slong len, slong prec)
{
    return 2 * _arb_vec_estimate_allocated_bytes(len, prec);
}

void _acb_vec_scalar_mul_2exp_si(acb_ptr res, acb_srcptr vec, slong len, slong c)
{
    slong i;
    for (i = 0; i < len; i++)
        acb_mul_2exp_si(res + i, vec + i, c);
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

IS_OP(_acb_vec_is_zero,   acb_srcptr, acb_is_zero)
IS_OP(_acb_vec_is_real,   acb_srcptr, acb_is_real)
IS_OP(_acb_vec_is_finite, acb_srcptr, acb_is_finite)

#define SET_OP_1(func_name, T, OP) \
void func_name(T ap, slong len) \
{                               \
    slong ix;                   \
                                \
    for (ix = 0; ix < len; ix++)\
        OP(ap + ix);            \
}

SET_OP_1(_acb_vec_zero,          acb_ptr, acb_zero)
SET_OP_1(_acb_vec_indeterminate, acb_ptr, acb_indeterminate)

#define SET_OP(func_name, S, T, OP) \
void func_name(S ap, T bp, slong len) \
{                               \
    slong ix;                   \
                                \
    for (ix = 0; ix < len; ix++)\
        OP(ap + ix, bp + ix);   \
}

SET_OP(_acb_vec_set,               acb_ptr, acb_srcptr, acb_set)
SET_OP(_acb_vec_swap,              acb_ptr, acb_ptr,    acb_swap)
SET_OP(_acb_vec_neg,               acb_ptr, acb_srcptr, acb_neg)
SET_OP(_acb_vec_trim,              acb_ptr, acb_srcptr, acb_trim)
SET_OP(_acb_vec_add_error_arf_vec, acb_ptr, arf_srcptr, acb_add_error_arf)
SET_OP(_acb_vec_add_error_mag_vec, acb_ptr, mag_srcptr, acb_add_error_mag)
SET_OP(_acb_vec_scalar_mul_onei,   acb_ptr, acb_srcptr, acb_mul_onei)

#define SET_PREC_OP(func_name, S, T, OP) \
void func_name(S ap, T bp, slong len, slong prec) \
{                               \
    slong ix;                   \
                                \
    for (ix = 0; ix < len; ix++)\
        OP(ap + ix, bp + ix, prec); \
}

SET_PREC_OP(_acb_vec_set_round, acb_ptr, acb_srcptr, acb_set_round)
SET_PREC_OP(_acb_vec_sqr,       acb_ptr, acb_srcptr, acb_sqr)

#define AORS_OP(func_name, S, T, OP) \
void func_name(S rp, T ap, T bp, slong len, slong prec) \
{                               \
    slong ix;                   \
                                \
    for (ix = 0; ix < len; ix++)\
        OP(rp + ix, ap + ix, bp + ix, prec); \
}

AORS_OP(_acb_vec_add, acb_ptr, acb_srcptr, acb_add)
AORS_OP(_acb_vec_sub, acb_ptr, acb_srcptr, acb_sub)

#define SCALAR_OP(func_name, S, T, U, OP) \
void func_name(S rp, T ap, slong len, U b, slong prec) \
{                               \
    slong ix;                   \
                                \
    for (ix = 0; ix < len; ix++)\
        OP(rp + ix, ap + ix, b, prec); \
}

SCALAR_OP(_acb_vec_scalar_mul_ui,   acb_ptr, acb_srcptr, ulong,        acb_mul_ui)
SCALAR_OP(_acb_vec_scalar_mul_fmpz, acb_ptr, acb_srcptr, const fmpz_t, acb_mul_fmpz)
SCALAR_OP(_acb_vec_scalar_mul_arb,  acb_ptr, acb_srcptr, const arb_t,  acb_mul_arb)
SCALAR_OP(_acb_vec_scalar_mul,      acb_ptr, acb_srcptr, const acb_t,  acb_mul)

SCALAR_OP(_acb_vec_scalar_div_ui,   acb_ptr, acb_srcptr, ulong,        acb_div_ui)
SCALAR_OP(_acb_vec_scalar_div_fmpz, acb_ptr, acb_srcptr, const fmpz_t, acb_div_fmpz)
SCALAR_OP(_acb_vec_scalar_div_arb,  acb_ptr, acb_srcptr, const arb_t,  acb_div_arb)
SCALAR_OP(_acb_vec_scalar_div,      acb_ptr, acb_srcptr, const acb_t,  acb_div)

SCALAR_OP(_acb_vec_scalar_addmul,   acb_ptr, acb_srcptr, const acb_t, acb_addmul)
SCALAR_OP(_acb_vec_scalar_submul,   acb_ptr, acb_srcptr, const acb_t, acb_submul)

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

COMPARISON_OP(_acb_vec_equal,               acb_srcptr, acb_srcptr, acb_equal)
COMPARISON_OP(_acb_vec_overlaps,            acb_srcptr, acb_srcptr, acb_overlaps)
COMPARISON_OP(_acb_vec_contains,            acb_srcptr, acb_srcptr, acb_contains)
COMPARISON_OP(_acb_vec_get_unique_fmpz_vec, fmpz *,     acb_srcptr, acb_get_unique_fmpz)
