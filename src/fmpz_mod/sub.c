/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod.h"

void _fmpz_mod_sub1(fmpz_t a, const fmpz_t b, const fmpz_t c,
                                                     const fmpz_mod_ctx_t ctx)
{
    mp_limb_t a0, b0, c0;

    FLINT_ASSERT(fmpz_mod_is_canonical(b, ctx));
    FLINT_ASSERT(fmpz_mod_is_canonical(c, ctx));

    b0 = fmpz_get_ui(b);
    c0 = fmpz_get_ui(c);
    a0 = nmod_sub(b0, c0, ctx->mod);
    fmpz_set_ui(a, a0);

    FLINT_ASSERT(fmpz_mod_is_canonical(a, ctx));
}

void _fmpz_mod_sub2s(fmpz_t a, const fmpz_t b, const fmpz_t c,
                                                     const fmpz_mod_ctx_t ctx)
{
    mp_limb_t a0, b0, c0;

    FLINT_ASSERT(fmpz_mod_is_canonical(b, ctx));
    FLINT_ASSERT(fmpz_mod_is_canonical(c, ctx));

    b0 = fmpz_get_ui(b);
    c0 = fmpz_get_ui(c);
    a0 = b0 - c0;
    fmpz_set_ui(a, a0);

    FLINT_ASSERT(fmpz_mod_is_canonical(a, ctx));
}

void _fmpz_mod_sub2(fmpz_t a, const fmpz_t b, const fmpz_t c,
                                                     const fmpz_mod_ctx_t ctx)
{
    mp_limb_t a2, a1, a0, b1, b0, c1, c0;

    FLINT_ASSERT(fmpz_mod_is_canonical(b, ctx));
    FLINT_ASSERT(fmpz_mod_is_canonical(c, ctx));

    fmpz_get_uiui(&b1, &b0, b);
    fmpz_get_uiui(&c1, &c0, c);
    sub_dddmmmsss(a2, a1, a0, 0, b1, b0, 0, c1, c0);
    if (a2 != 0)
    {
        add_ssaaaa(a1, a0, a1, a0, ctx->n_limbs[1], ctx->n_limbs[0]);
    }
    fmpz_set_uiui(a, a1, a0);

    FLINT_ASSERT(fmpz_mod_is_canonical(a, ctx));
}

void _fmpz_mod_subN(fmpz_t a, const fmpz_t b, const fmpz_t c,
                                                     const fmpz_mod_ctx_t ctx)
{
    FLINT_ASSERT(fmpz_mod_is_canonical(b, ctx));
    FLINT_ASSERT(fmpz_mod_is_canonical(c, ctx));
    fmpz_sub(a, b, c);
    if (fmpz_sgn(a) < 0)
    {
        fmpz_add(a, a, ctx->n);
    }
    FLINT_ASSERT(fmpz_mod_is_canonical(a, ctx));
}
