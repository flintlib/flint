/*
    Copyright (C) 2026 Brian Heckel

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    Batched Montgomery 4x4 matrix power A^exp over Z/nZ (n odd), AVX2.

    Batched Montgomery only pays off when many multiplies happen between a single
    pair of to/from-Montgomery conversions, so this converts A once, runs the
    whole square-and-multiply chain in Montgomery form, and converts back once.
    The limb count s = ceil((bits+2)/31) is a runtime parameter (see
    mpn_mod_batched_mont.h), so the work scales with the modulus size.

    Falls back to the generic gr_mat_pow_ui for anything this path does not
    handle: non-4x4 matrices, even moduli, moduli too large for the batched
    representation, or a non-AVX2 build.
*/

// TODO: Not sure if flint has a specific function for memset or
// an alternative
#include <string.h>
#include "mpn_mod_batched_mont.h"

#ifdef __AVX2__

/* The square-and-multiply buffers are rotated as plain (non-const) pointers;
   cast them to the const-qualified bm_mat the kernel reads.  A direct implicit
   conversion of a pointer-to-array to a const-qualified one is a -Wpedantic
   warning in C before C23. */
#define BM_CONST(m) ((const ulong (*)[BM_DIM][BM_MAXLIMB]) (m))

int
mpn_mod_mat_pow_ui_batched_mont(gr_mat_t res, const gr_mat_t A, ulong exp, gr_ctx_t ctx)
{
    bm_ctx_t bc;
    bm_mat buf0, buf1, buf2;
    ulong (*base)[BM_DIM][BM_MAXLIMB];
    ulong (*acc)[BM_DIM][BM_MAXLIMB];
    ulong (*scr)[BM_DIM][BM_MAXLIMB];
    ulong (*sw)[BM_DIM][BM_MAXLIMB];
    slong bits = MPN_MOD_CTX_MODULUS_BITS(ctx);
    slong s = (bits + 2 + BM_BITS - 1) / BM_BITS;
    slong bit;

    // check if we don't match conditions for batched montgomery
    if (A->r != 4 || A->c != 4 || res->r != 4 || res->c != 4
        || !(MPN_MOD_CTX_MODULUS(ctx)[0] & 1) || s > BM_MAXLIMB)
        return gr_mat_pow_ui(res, A, exp, ctx);

    if (exp == 0)
        return gr_mat_one(res, ctx);

    memset(buf0, 0, sizeof(bm_mat));
    memset(buf1, 0, sizeof(bm_mat));
    memset(buf2, 0, sizeof(bm_mat));

    bm_ctx_init(&bc, ctx);
    bm_load(buf0, A, ctx, &bc);
    memcpy(buf1, buf0, sizeof(bm_mat));
    base = buf0; acc = buf1; scr = buf2;

    // left-to-right square-and-multiply; the top set bit is already in acc
    for (bit = (slong) FLINT_BIT_COUNT(exp) - 2; bit >= 0; bit--)
    {
        bm_matmul(scr, BM_CONST(acc), BM_CONST(acc), &bc);       /* square */
        sw = acc; acc = scr; scr = sw;
        if ((exp >> bit) & 1)
        {
            bm_matmul(scr, BM_CONST(acc), BM_CONST(base), &bc);  /* multiply by base */
            sw = acc; acc = scr; scr = sw;
        }
    }

    bm_store(res, BM_CONST(acc), ctx, &bc);
    return GR_SUCCESS;
}

#else

int
mpn_mod_mat_pow_ui_batched_mont(gr_mat_t res, const gr_mat_t A, ulong exp, gr_ctx_t ctx)
{
    return gr_mat_pow_ui(res, A, exp, ctx);
}

#endif
