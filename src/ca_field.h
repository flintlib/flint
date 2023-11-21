/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef CA_FIELD_H
#define CA_FIELD_H

#ifdef CA_FIELD_INLINES_C
#define CA_FIELD_INLINE
#else
#define CA_FIELD_INLINE static inline
#endif

#include "ca.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Types *********************************************************************/

/* note: types and macros are defined in ca.h since they are needed there */

#define CA_FIELD_HASH_C UWORD(100003)

void ca_field_init_qq(ca_field_t K, ca_ctx_t ctx);
void ca_field_init_nf(ca_field_t K, const qqbar_t x, ca_ctx_t ctx);
void ca_field_init_const(ca_field_t K, calcium_func_code func, ca_ctx_t ctx);
void ca_field_init_fx(ca_field_t K, calcium_func_code func, const ca_t x, ca_ctx_t ctx);
void ca_field_init_fxy(ca_field_t K, calcium_func_code func, const ca_t x, const ca_t y, ca_ctx_t ctx);
void ca_field_init_multi(ca_field_t K, slong len, ca_ctx_t ctx);
void ca_field_clear(ca_field_t K, ca_ctx_t ctx);

void ca_field_set_ext(ca_field_t K, slong i, ca_ext_srcptr x, ca_ctx_t ctx);

void ca_field_print(const ca_field_t K, ca_ctx_t ctx);

int ca_field_cmp(const ca_field_t K1, const ca_field_t K2, ca_ctx_t ctx);

void ca_field_build_ideal(ca_field_t K, ca_ctx_t ctx);
void ca_field_build_ideal_erf(ca_field_t K, ca_ctx_t ctx);
void ca_field_build_ideal_gamma(ca_field_t K, ca_ctx_t ctx);

void ca_field_cache_init(ca_field_cache_t cache, ca_ctx_t ctx);
void ca_field_cache_clear(ca_field_cache_t cache, ca_ctx_t ctx);
ca_field_ptr ca_field_cache_insert_ext(ca_field_cache_t cache, ca_ext_struct ** x, slong length, ca_ctx_t ctx);

#ifdef __cplusplus
}
#endif

#endif

