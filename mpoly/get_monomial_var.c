/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "mpoly.h"

/* this file does not need to change with new orderings */

ulong mpoly_get_monomial_var_exp_ui_sp(const ulong * poly_exps,
                        slong var, flint_bitcnt_t bits, const mpoly_ctx_t mctx)
{
    slong offset, shift;

    mpoly_gen_offset_shift_sp(&offset, &shift, var, bits, mctx);

    return (poly_exps[offset] >> shift) & ((-UWORD(1)) >> (FLINT_BITS - bits));    
}

ulong mpoly_get_monomial_var_exp_ui_mp(const ulong * poly_exps,
                        slong var, flint_bitcnt_t bits, const mpoly_ctx_t mctx)
{
    slong offset;
    slong j;
    ulong wpf = bits/FLINT_BITS;
    ulong r, check;

    offset = mpoly_gen_offset_mp(var, bits, mctx);

    r = poly_exps[offset + 0];
    check = 0;
    for (j = 1; j < wpf; j++)
        check |= poly_exps[offset + j];

    if (check != 0)
        flint_throw(FLINT_ERROR, "Exponent does not fit a ulong.");

    return r;
}

slong mpoly_get_monomial_var_exp_si_mp(const ulong * poly_exps,
                        slong var, flint_bitcnt_t bits, const mpoly_ctx_t mctx)
{
    slong offset;
    slong j;
    ulong wpf = bits/FLINT_BITS;
    ulong r, check;

    offset = mpoly_gen_offset_mp(var, bits, mctx);

    r = poly_exps[offset + 0];
    check = FLINT_SIGN_EXT(r);
    for (j = 1; j < wpf; j++)
        check |= poly_exps[offset + j];

    if (check != 0)
        flint_throw(FLINT_ERROR, "Exponent does not fit an slong.");

    return (slong) r;
}
