/*
    Copyright (C) 2026 Brian Heckel

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    One-shot batched Montgomery 4x4 matrix multiply over Z/nZ (n odd, AVX2),
    using lazy reduction: for each output entry the 4 raw products are summed
    into a wide accumulator and reduced once. Summing 4 products reaches
    4n^2, so single-subtract REDC is valid only for n < 2^246;
    the dispatch guard enforces that.

    NOTE: unlike a chain, a one-shot matmul cannot amortize the to/from-
    Montgomery conversion (A, B in and C out), which adds ~12 CIOS multiplies on
    top of the matmul. So this is expected to be competitive only when the
    conversion is amortized (matrix powers / long chains); it is provided here
    mainly to measure the one-shot case. Falls back to gr_mat_mul_classical for
    non-4x4, even n, n >= 2^246, or a non-AVX2 build.

    See mpn_mod_batched_mont.h for the batched kernel and its representation.
*/

#include "mpn_mod_batched_mont.h"

#ifdef __AVX2__

/* cast a (non-const) bm_mat to the const-qualified pointer the kernel reads;
   the implicit conversion is a -Wpedantic warning in C before C23 */
#define BM_CONST(m) ((const ulong (*)[BM_DIM][BM_MAXLIMB]) (m))

int
mpn_mod_mat_mul_batched_mont(gr_mat_t C, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx)
{
    slong bits = MPN_MOD_CTX_MODULUS_BITS(ctx);
    slong s = (bits + 2 + BM_BITS - 1) / BM_BITS;

    // check for either lazy montgomery reduction or regular matrix multiplication
    if (A->r == 4 && A->c == 4 && B->r == 4 && B->c == 4 && C->r == 4 && C->c == 4
        && (MPN_MOD_CTX_MODULUS(ctx)[0] & 1) && s <= BM_MAXLIMB)
    {
        bm_ctx_t bc;
        bm_mat mA, mB, mC;
        bm_ctx_init(&bc, ctx);
        bm_load(mA, A, ctx, &bc);
        bm_load(mB, B, ctx, &bc);
        bm_matmul(mC, BM_CONST(mA), BM_CONST(mB), &bc);
        bm_store(C, BM_CONST(mC), ctx, &bc);
        return GR_SUCCESS;
    }

    return gr_mat_mul_classical(C, A, B, ctx);
}

#else

int
mpn_mod_mat_mul_batched_mont(gr_mat_t C, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx)
{
    return gr_mat_mul_classical(C, A, B, ctx);
}

#endif
