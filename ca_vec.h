/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#ifndef CA_VEC_H
#define CA_VEC_H

#ifdef CA_VEC_INLINES_C
#define CA_VEC_INLINE
#else
#define CA_VEC_INLINE static __inline__
#endif

#include <stdio.h>
#include "flint/flint.h"
#include "flint/fmpz_poly.h"
#include "flint/fmpq_poly.h"
#include "arb_poly.h"
#include "acb_poly.h"
#include "antic/nf.h"
#include "antic/nf_elem.h"
#include "ca.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Vector object */

typedef struct
{
    ca_ptr coeffs;
    slong length;
    slong alloc;
}
ca_vec_struct;

typedef ca_vec_struct ca_vec_t[1];

#define ca_vec_entry(vec, i) ((vec)->coeffs + (i))

/* Memory management */

ca_ptr _ca_vec_init(slong len, ca_ctx_t ctx);
void ca_vec_init(ca_vec_t vec, slong len, ca_ctx_t ctx);

void _ca_vec_clear(ca_ptr v, slong len, ca_ctx_t ctx);
void ca_vec_clear(ca_vec_t vec, ca_ctx_t ctx);

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

void ca_vec_set_length(ca_vec_t res, slong len, ca_ctx_t ctx);

/* Assignment */

void _ca_vec_set(ca_ptr res, ca_srcptr src, slong len, ca_ctx_t ctx);
void ca_vec_set(ca_vec_t res, const ca_vec_t src, ca_ctx_t ctx);

/* Special vectors */

void _ca_vec_zero(ca_ptr res, slong len, ca_ctx_t ctx);
void ca_vec_zero(ca_vec_t res, slong len, ca_ctx_t ctx);

/* Arithmetic */

void _ca_vec_neg(ca_ptr res, ca_srcptr src, slong len, ca_ctx_t ctx);
void ca_vec_neg(ca_vec_t res, const ca_vec_t src, ca_ctx_t ctx);

void _ca_vec_add(ca_ptr res, ca_srcptr vec1, ca_srcptr vec2, slong len, ca_ctx_t ctx);
void _ca_vec_sub(ca_ptr res, ca_srcptr vec1, ca_srcptr vec2, slong len, ca_ctx_t ctx);
void _ca_vec_scalar_mul_ca(ca_ptr res, ca_srcptr src, slong len, const ca_t c, ca_ctx_t ctx);
void _ca_vec_scalar_addmul_ca(ca_ptr res, ca_srcptr vec, slong len, const ca_t c, ca_ctx_t ctx);


#ifdef __cplusplus
}
#endif

#endif
