/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef CA_VEC_H
#define CA_VEC_H

#ifdef CA_VEC_INLINES_C
#define CA_VEC_INLINE
#else
#define CA_VEC_INLINE static inline
#endif

#include "ca.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Vector object */

typedef struct
{
    ca_ptr entries;
    slong alloc;
    slong length;
}
ca_vec_struct;

typedef ca_vec_struct ca_vec_t[1];

#define ca_vec_entry(vec, i) ((vec)->entries + (i))

CA_VEC_INLINE ca_ptr
ca_vec_entry_ptr(ca_vec_t vec, slong i)
{
    return ca_vec_entry(vec, i);
}

/* Memory management */

ca_ptr _ca_vec_init(slong len, ca_ctx_t ctx);
void ca_vec_init(ca_vec_t vec, slong len, ca_ctx_t ctx);

void _ca_vec_clear(ca_ptr v, slong len, ca_ctx_t ctx);
void ca_vec_clear(ca_vec_t vec, ca_ctx_t ctx);

void _ca_vec_swap(ca_ptr vec1, ca_ptr vec2, slong len, ca_ctx_t ctx);

CA_VEC_INLINE void
ca_vec_swap(ca_vec_t vec1, ca_vec_t vec2, ca_ctx_t ctx)
{
    ca_vec_struct t = *vec1;
    *vec1 = *vec2;
    *vec2 = t;
}

/* Length */

CA_VEC_INLINE
slong ca_vec_length(const ca_vec_t vec, ca_ctx_t ctx)
{
    return vec->length;
}

void _ca_vec_fit_length(ca_vec_t vec, slong len, ca_ctx_t ctx);
void ca_vec_set_length(ca_vec_t res, slong len, ca_ctx_t ctx);

/* Assignment */

void _ca_vec_set(ca_ptr res, ca_srcptr src, slong len, ca_ctx_t ctx);
void ca_vec_set(ca_vec_t res, const ca_vec_t src, ca_ctx_t ctx);

/* Special vectors */

void _ca_vec_zero(ca_ptr res, slong len, ca_ctx_t ctx);
void ca_vec_zero(ca_vec_t res, slong len, ca_ctx_t ctx);

CA_VEC_INLINE
void _ca_vec_unknown(ca_ptr vec, slong len, ca_ctx_t ctx)
{
    slong i;
    for (i = 0; i < len; i++)
        ca_unknown(vec + i, ctx);
}

CA_VEC_INLINE
void _ca_vec_undefined(ca_ptr vec, slong len, ca_ctx_t ctx)
{
    slong i;
    for (i = 0; i < len; i++)
        ca_undefined(vec + i, ctx);
}

/* Input and output */

void ca_vec_print(const ca_vec_t vec, ca_ctx_t ctx);
void ca_vec_printn(const ca_vec_t vec, slong digits, ca_ctx_t ctx);

/* List operations */

CA_VEC_INLINE void
ca_vec_append(ca_vec_t vec, const ca_t f, ca_ctx_t ctx)
{
    _ca_vec_fit_length(vec, vec->length + 1, ctx);
    ca_set(vec->entries + vec->length, f, ctx);
    vec->length++;
}

/* Arithmetic */

void _ca_vec_neg(ca_ptr res, ca_srcptr src, slong len, ca_ctx_t ctx);
void ca_vec_neg(ca_vec_t res, const ca_vec_t src, ca_ctx_t ctx);

void _ca_vec_add(ca_ptr res, ca_srcptr vec1, ca_srcptr vec2, slong len, ca_ctx_t ctx);
void _ca_vec_sub(ca_ptr res, ca_srcptr vec1, ca_srcptr vec2, slong len, ca_ctx_t ctx);
void _ca_vec_scalar_mul_ca(ca_ptr res, ca_srcptr src, slong len, const ca_t c, ca_ctx_t ctx);
void _ca_vec_scalar_div_ca(ca_ptr res, ca_srcptr src, slong len, const ca_t c, ca_ctx_t ctx);
void _ca_vec_scalar_addmul_ca(ca_ptr res, ca_srcptr vec, slong len, const ca_t c, ca_ctx_t ctx);
void _ca_vec_scalar_submul_ca(ca_ptr res, ca_srcptr vec, slong len, const ca_t c, ca_ctx_t ctx);

/* Comparisons and predicates */

truth_t _ca_vec_check_is_zero(ca_srcptr vec, slong len, ca_ctx_t ctx);

/* Internal representation */

CA_VEC_INLINE int
_ca_vec_is_fmpq_vec(ca_srcptr vec, slong len, ca_ctx_t ctx)
{
    slong i;
    for (i = 0; i < len; i++)
        if (!CA_IS_QQ(vec + i, ctx))
                return 0;
    return 1;
}

CA_VEC_INLINE int
_ca_vec_fmpq_vec_is_fmpz_vec(ca_srcptr vec, slong len, ca_ctx_t ctx)
{
    slong i;
    for (i = 0; i < len; i++)
        if (!fmpz_is_one(CA_FMPQ_DENREF(vec + i)))
                return 0;
    return 1;
}

CA_VEC_INLINE void
_ca_vec_fmpq_vec_get_fmpz_vec_den(fmpz * c, fmpz_t den, ca_srcptr vec, slong len, ca_ctx_t ctx)
{
    slong i;

    fmpz_one(den);

    if (_ca_vec_fmpq_vec_is_fmpz_vec(vec, len, ctx))
    {
        for (i = 0; i < len; i++)
            fmpz_set(c + i, CA_FMPQ_NUMREF(vec + i));
    }
    else
    {
        for (i = 0; i < len; i++)
            fmpz_lcm(den, den, CA_FMPQ_DENREF(vec + i));

        for (i = 0; i < len; i++)
        {
            fmpz_divexact(c + i, den, CA_FMPQ_DENREF(vec + i));
            fmpz_mul(c + i, c + i, CA_FMPQ_NUMREF(vec + i));
        }
    }
}

CA_VEC_INLINE void
_ca_vec_set_fmpz_vec_div_fmpz(ca_ptr res, const fmpz * v, const fmpz_t den, slong len, ca_ctx_t ctx)
{
    slong i;

    if (fmpz_is_one(den))
    {
        for (i = 0; i < len; i++)
            ca_set_fmpz(res + i, v + i, ctx);
    }
    else
    {
        for (i = 0; i < len; i++)
        {
            ca_set_fmpz(res + i, v + i, ctx);       /* todo: optimize */
            ca_div_fmpz(res + i, res + i, den, ctx);
        }
    }
}

#ifdef __cplusplus
}
#endif

#endif
