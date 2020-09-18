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

/* Input and output */

void ca_poly_print(const ca_poly_t poly, ca_ctx_t ctx);
void ca_poly_printn(const ca_poly_t poly, slong digits, ca_ctx_t ctx);


#ifdef __cplusplus
}
#endif

#endif
