/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#ifndef CA_POLY_H
#define CA_POLY_H

#ifdef CA_POLY_INLINES_C
#define CA_POLY_INLINE
#else
#define CA_POLY_INLINE static __inline__
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

/* Polynomial object */

typedef struct
{
    ca_ptr coeffs;
    slong length;
    slong alloc;
}
ca_poly_struct;

typedef ca_poly_struct ca_poly_t[1];


/* Memory management */

void ca_poly_init(ca_poly_t poly, ca_ctx_t ctx);
void ca_poly_init2(ca_poly_t poly, slong len, ca_ctx_t ctx);
void ca_poly_clear(ca_poly_t poly, ca_ctx_t ctx);
void ca_poly_fit_length(ca_poly_t poly, slong len, ca_ctx_t ctx);
void _ca_poly_set_length(ca_poly_t poly, slong len, ca_ctx_t ctx);
void _ca_poly_normalise(ca_poly_t poly, ca_ctx_t ctx);

CA_POLY_INLINE void
ca_poly_swap(ca_poly_t poly1, ca_poly_t poly2, ca_ctx_t ctx)
{
    ca_poly_struct t = *poly1;
    *poly1 = *poly2;
    *poly2 = t;
}

/* Random generation */

void ca_poly_randtest(ca_poly_t poly, flint_rand_t state, slong len, slong depth, slong bits, ca_ctx_t ctx);
void ca_poly_randtest_rational(ca_poly_t poly, flint_rand_t state, slong len, slong bits, ca_ctx_t ctx);

/* Basic operations */

/* todo: document */
void _ca_poly_reverse(ca_ptr res, ca_srcptr poly, slong len, slong n, ca_ctx_t ctx);

/* Input and output */

void ca_poly_print(const ca_poly_t poly, ca_ctx_t ctx);
void ca_poly_printn(const ca_poly_t poly, slong digits, ca_ctx_t ctx);

/* Comparisons */

/* todo: document */
truth_t _ca_poly_check_equal(ca_srcptr poly1, slong len1, ca_srcptr poly2, slong len2, ca_ctx_t ctx);
truth_t ca_poly_check_equal(const ca_poly_t poly1, const ca_poly_t poly2, ca_ctx_t ctx);


/* Arithmetic */

/* todo: document */
CA_POLY_INLINE void
ca_poly_mul_ca(ca_poly_t res, const ca_poly_t poly, const ca_t c, ca_ctx_t ctx)
{
    ca_poly_fit_length(res, poly->length, ctx);
    ca_vec_scalar_mul_ca(res->coeffs, poly->coeffs, poly->length, c, ctx);
    _ca_poly_set_length(res, poly->length, ctx);
    _ca_poly_normalise(res, ctx);
}

#ifdef __cplusplus
}
#endif

#endif
