/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef CA_EXT_H
#define CA_EXT_H

#ifdef CA_EXT_INLINES_C
#define CA_EXT_INLINE
#else
#define CA_EXT_INLINE static inline
#endif

#include "ca.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Types *********************************************************************/

/* note: types and macros are defined in ca.h since they are needed there */

void ca_ext_init_qqbar(ca_ext_t res, const qqbar_t x, ca_ctx_t ctx);
void ca_ext_init_const(ca_ext_t res, calcium_func_code func, ca_ctx_t ctx);
void ca_ext_init_fx(ca_ext_t res, calcium_func_code func, const ca_t x, ca_ctx_t ctx);
void ca_ext_init_fxy(ca_ext_t res, calcium_func_code func, const ca_t x, const ca_t y, ca_ctx_t ctx);
void ca_ext_init_fxn(ca_ext_t res, calcium_func_code func, ca_srcptr x, slong nargs, ca_ctx_t ctx);

/* todo: this could avoid rehashing, ... */
CA_EXT_INLINE void
ca_ext_init_set(ca_ext_t res, const ca_ext_t x, ca_ctx_t ctx)
{
    if (CA_EXT_HEAD(x) == CA_QQBar)
    {
        ca_ext_init_qqbar(res, CA_EXT_QQBAR(x), ctx);
    }
    else
    {
        ca_ext_init_fxn(res, CA_EXT_HEAD(x), CA_EXT_FUNC_ARGS(x), CA_EXT_FUNC_NARGS(x), ctx);
    }
}

void ca_ext_clear(ca_ext_t res, ca_ctx_t ctx);

CA_EXT_INLINE slong ca_ext_nargs(const ca_ext_t x, ca_ctx_t ctx)
{
    if (CA_EXT_HEAD(x) == CA_QQBar)
        return 0;
    else
        return CA_EXT_FUNC_NARGS(x);
}

CA_EXT_INLINE void ca_ext_get_arg(ca_t res, const ca_ext_t x, slong i, ca_ctx_t ctx)
{
    if (CA_EXT_HEAD(x) == CA_QQBar || i < 0 || i >= CA_EXT_FUNC_NARGS(x))
    {
        flint_throw(FLINT_ERROR, "ca_ext_get_arg: index out of range\n");
    }
    else
    {
        ca_set(res, CA_EXT_FUNC_ARGS(x) + i, ctx);
    }
}

CA_EXT_INLINE ulong ca_ext_hash(const ca_ext_t x, ca_ctx_t ctx)
{
    return CA_EXT_HASH(x);
}

int ca_ext_equal_repr(const ca_ext_t x, const ca_ext_t y, ca_ctx_t ctx);
int ca_ext_cmp_repr(const ca_ext_t x, const ca_ext_t y, ca_ctx_t ctx);

void ca_ext_print(const ca_ext_t x, ca_ctx_t ctx);

void ca_ext_get_acb_raw(acb_t res, ca_ext_t x, slong prec, ca_ctx_t ctx);


void ca_ext_cache_init(ca_ext_cache_t cache, ca_ctx_t ctx);
void ca_ext_cache_clear(ca_ext_cache_t cache, ca_ctx_t ctx);
ca_ext_ptr ca_ext_cache_insert(ca_ext_cache_t cache, const ca_ext_t x, ca_ctx_t ctx);

#ifdef __cplusplus
}
#endif

#endif

