/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fft_small.h"

FLINT_TLS_PREFIX mpn_ctx_t default_mpn_ctx;
FLINT_TLS_PREFIX int default_mpn_ctx_initialized = 0;

#define DEFAULT_PRIME UWORD(0x0003f00000000001)

static void
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
mpn_mul_default_mpn_ctx(nn_ptr r1, nn_srcptr i1, slong n1, nn_srcptr i2, slong n2)
{
    mpn_ctx_mpn_mul(get_default_mpn_ctx(), r1, i1, n1, i2, n2);
}

void
_nmod_poly_mul_mid_default_mpn_ctx(nn_ptr res, slong zl, slong zh, nn_srcptr a, slong an, nn_srcptr b, slong bn, nmod_t mod)
{
    _nmod_poly_mul_mid_mpn_ctx(res, zl, zh, a, an, b, bn, mod, get_default_mpn_ctx());
}

int
_fmpz_poly_mul_mid_default_mpn_ctx(fmpz * res, slong zl, slong zh, const fmpz * a, slong an, const fmpz * b, slong bn)
{
    return _fmpz_poly_mul_mid_mpn_ctx(res, zl, zh, a, an, b, bn, get_default_mpn_ctx());
}
