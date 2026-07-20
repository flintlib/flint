/*
    Copyright (C) 2026 FLINT contributors

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpoly.h"

/*
    A recurring idiom throughout the mpoly modules: an in-progress
    computation discovers (via an overflowed guard bit) that the current
    packing width is not wide enough, and needs to widen it by at least
    one bit, reallocate storage sized for that wider packing, and repack
    the monomials computed so far into it. This function does exactly
    that: it increases *bits_ptr (rounded up to whatever width mctx
    actually supports via mpoly_fix_bits), recomputes *N_ptr for the new
    width, and returns a freshly flint_malloc'd array holding `len`
    monomials repacked from `old_exps` (which was packed at `old_bits`),
    with room for `alloc` monomials in total.

    The caller owns the returned array. `old_exps` is left untouched --
    it is very often part of an input polynomial the caller does not
    own, so freeing it (when it is safe to do so) remains the caller's
    responsibility, exactly as before this helper existed.
*/
ulong * mpoly_monomials_repack_wider(flint_bitcnt_t * bits_ptr, slong * N_ptr,
        const ulong * old_exps, flint_bitcnt_t old_bits, slong len,
        slong alloc, const mpoly_ctx_t mctx)
{
    ulong * new_exps;

    *bits_ptr = mpoly_fix_bits(*bits_ptr + 1, mctx);
    *N_ptr = mpoly_words_per_exp(*bits_ptr, mctx);

    new_exps = FLINT_ARRAY_ALLOC(*N_ptr * alloc, ulong);
    mpoly_repack_monomials(new_exps, *bits_ptr, old_exps, old_bits, len, mctx);

    return new_exps;
}

/*
    As mpoly_monomials_repack_wider, and additionally grows *cmpmask_ptr
    (via flint_realloc) to match the new bits/N and refreshes its
    contents. This covers the common combined idiom where a cmpmask
    array tracking the same packing width needs to stay in sync with
    the repack; *cmpmask_ptr must already be a flint_malloc/flint_realloc
    (not TMP_ALLOC) allocation, since it is grown with flint_realloc here.
*/
ulong * mpoly_monomials_repack_wider_cmpmask(flint_bitcnt_t * bits_ptr,
        slong * N_ptr, ulong ** cmpmask_ptr, const ulong * old_exps,
        flint_bitcnt_t old_bits, slong len, slong alloc,
        const mpoly_ctx_t mctx)
{
    ulong * new_exps;

    *bits_ptr = mpoly_fix_bits(*bits_ptr + 1, mctx);
    *N_ptr = mpoly_words_per_exp(*bits_ptr, mctx);

    *cmpmask_ptr = (ulong *) flint_realloc(*cmpmask_ptr, *N_ptr * sizeof(ulong));
    mpoly_get_cmpmask(*cmpmask_ptr, *N_ptr, *bits_ptr, mctx);

    new_exps = FLINT_ARRAY_ALLOC(*N_ptr * alloc, ulong);
    mpoly_repack_monomials(new_exps, *bits_ptr, old_exps, old_bits, len, mctx);

    return new_exps;
}
