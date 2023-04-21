/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fft_small.h"

FLINT_TLS_PREFIX mpn_ctx_t default_mpn_ctx;
FLINT_TLS_PREFIX int default_mpn_ctx_initialized = 0;

#define DEFAULT_PRIME UWORD(0x0003f00000000001)

void
mpn_ctx_cleanup(void)
{
    if (default_mpn_ctx_initialized)
    {
        default_mpn_ctx_initialized = 0;
        mpn_ctx_clear(default_mpn_ctx);
    }
}

mpn_ctx_struct * get_default_mpn_ctx(void)
{
    if (!default_mpn_ctx_initialized)
    {
        mpn_ctx_init(default_mpn_ctx, DEFAULT_PRIME);
        flint_register_cleanup_function(mpn_ctx_cleanup);
        default_mpn_ctx_initialized = 1;
    }

    return default_mpn_ctx;
}

void
mpn_mul_default_mpn_ctx(mp_ptr r1, mp_srcptr i1, mp_size_t n1, mp_srcptr i2, mp_size_t n2)
{
    mpn_ctx_mpn_mul(get_default_mpn_ctx(), r1, i1, n1, i2, n2);
}

#if defined(__AVX2__)

void
_nmod_poly_mul_mid_default_mpn_ctx(mp_ptr res, slong zl, slong zh, mp_srcptr a, slong an, mp_srcptr b, slong bn, nmod_t mod)
{
    _nmod_poly_mul_mid_mpn_ctx(res, zl, zh, a, an, b, bn, mod, get_default_mpn_ctx());
}

#else

void _nmod_poly_divrem_mpn_ctx(
    ulong* q,
    ulong* r,
    const ulong* a, ulong an,
    const ulong* b, ulong bn,
    nmod_t mod,
    mpn_ctx_t R)
{
    flint_abort();
}

void
_nmod_poly_mul_mid_default_mpn_ctx(mp_ptr res, slong zl, slong zh, mp_srcptr a, slong an, mp_srcptr b, slong bn, nmod_t mod)
{
    flint_abort();
}

#endif

