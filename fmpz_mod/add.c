/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod.h"

void _fmpz_mod_add1(fmpz_t a, const fmpz_t b, const fmpz_t c,
                                                     const fmpz_mod_ctx_t ctx)
{
    ulong a0, b0, c0;

    FLINT_ASSERT(fmpz_mod_is_canonical(b, ctx));
    FLINT_ASSERT(fmpz_mod_is_canonical(c, ctx));

    b0 = fmpz_get_ui(b);
    c0 = fmpz_get_ui(c);
    a0 = nmod_add(b0, c0, ctx->mod);
    fmpz_set_ui(a, a0);

    FLINT_ASSERT(fmpz_mod_is_canonical(a, ctx));
}

void _fmpz_mod_add2s(fmpz_t a, const fmpz_t b, const fmpz_t c,
                                                     const fmpz_mod_ctx_t ctx)
{
    mp_limb_t a0, b0, c0;

    FLINT_ASSERT(fmpz_mod_is_canonical(b, ctx));
    FLINT_ASSERT(fmpz_mod_is_canonical(c, ctx));

    b0 = fmpz_get_ui(b);
    c0 = fmpz_get_ui(c);
    a0 = b0 + c0;
    fmpz_set_ui(a, a0);

    FLINT_ASSERT(fmpz_mod_is_canonical(a, ctx));
}

void _fmpz_mod_add2(fmpz_t a, const fmpz_t b, const fmpz_t c,
                                                     const fmpz_mod_ctx_t ctx)
{
    mp_limb_t t2, t1, t0, a2, a1, a0, b1, b0, c1, c0;

    FLINT_ASSERT(fmpz_mod_is_canonical(b, ctx));
    FLINT_ASSERT(fmpz_mod_is_canonical(c, ctx));

    fmpz_get_uiui(&b1, &b0, b);
    fmpz_get_uiui(&c1, &c0, c);
    add_sssaaaaaa(a2, a1, a0, 0, b1, b0, 0, c1, c0);
    /* at most one subtraction of n */
    sub_dddmmmsss(t2, t1, t0, a2, a1, a0, 0, ctx->n_limbs[1], ctx->n_limbs[0]);
    if ((slong)(t2) >= 0)
    {
        a1 = t1;
        a0 = t0;
    }
    fmpz_set_uiui(a, a1, a0);

    FLINT_ASSERT(fmpz_mod_is_canonical(a, ctx));
}

void _fmpz_mod_addN(fmpz_t a, const fmpz_t b, const fmpz_t c,
                                                     const fmpz_mod_ctx_t ctx)
{
    FLINT_ASSERT(fmpz_mod_is_canonical(b, ctx));
    FLINT_ASSERT(fmpz_mod_is_canonical(c, ctx));
    fmpz_add(a, b, c);
    if (fmpz_cmp(a, ctx->n) >= 0)
    {
        fmpz_sub(a, a, ctx->n);
    }
    FLINT_ASSERT(fmpz_mod_is_canonical(a, ctx));
}
