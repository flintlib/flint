/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod.h"
#include "fmpz.h"
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
                                                     const fmpz_mod_ctx_t FLINT_UNUSED(ctx))
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

void fmpz_mod_sub_fmpz(fmpz_t a, const fmpz_t b, const fmpz_t c,
                                                      const fmpz_mod_ctx_t ctx)
{
    FLINT_ASSERT(fmpz_mod_is_canonical(b, ctx));

    fmpz_sub(a, b, c);
    fmpz_mod(a, a, ctx->n);

    FLINT_ASSERT(fmpz_mod_is_canonical(a, ctx));
}

void fmpz_mod_sub_ui(fmpz_t a, const fmpz_t b, ulong c,
                                                      const fmpz_mod_ctx_t ctx)
{
    FLINT_ASSERT(fmpz_mod_is_canonical(b, ctx));

    fmpz_sub_ui(a, b, c);
    fmpz_mod(a, a, ctx->n);

    FLINT_ASSERT(fmpz_mod_is_canonical(a, ctx));
}

void fmpz_mod_sub_si(fmpz_t a, const fmpz_t b, slong c,
                                                      const fmpz_mod_ctx_t ctx)
{
    FLINT_ASSERT(fmpz_mod_is_canonical(b, ctx));

    fmpz_sub_si(a, b, c);
    fmpz_mod(a, a, ctx->n);

    FLINT_ASSERT(fmpz_mod_is_canonical(a, ctx));
}

void fmpz_mod_fmpz_sub(fmpz_t a, const fmpz_t b, const fmpz_t c,
                                                      const fmpz_mod_ctx_t ctx)
{
    FLINT_ASSERT(fmpz_mod_is_canonical(c, ctx));

    fmpz_sub(a, b, c);
    fmpz_mod(a, a, ctx->n);

    FLINT_ASSERT(fmpz_mod_is_canonical(a, ctx));
}

void fmpz_mod_ui_sub(fmpz_t a, ulong b, const fmpz_t c,
                                                      const fmpz_mod_ctx_t ctx)
{
    FLINT_ASSERT(fmpz_mod_is_canonical(c, ctx));

    fmpz_sub_ui(a, c, b);
    fmpz_neg(a, a);
    fmpz_mod(a, a, ctx->n);

    FLINT_ASSERT(fmpz_mod_is_canonical(a, ctx));
}

void fmpz_mod_si_sub(fmpz_t a, slong b, const fmpz_t c,
                                                      const fmpz_mod_ctx_t ctx)
{
    FLINT_ASSERT(fmpz_mod_is_canonical(c, ctx));

    fmpz_sub_si(a, c, b);
    fmpz_neg(a, a);
    fmpz_mod(a, a, ctx->n);

    FLINT_ASSERT(fmpz_mod_is_canonical(a, ctx));
}
