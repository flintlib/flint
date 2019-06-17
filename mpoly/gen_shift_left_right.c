/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "mpoly.h"

/*
    shift_right (resp. left) divides (resp. multiplies) A in place by x^amount,
    where x is the variable of index var.
    No reallocation on Abits is performed: The shift right version asserts
    that the division was exact, and the shift left version tries to assert
    that the result fits in the same number of bits.
*/

void _mpoly_gen_shift_right(ulong * Aexp, flint_bitcnt_t Abits, slong Alength,
                               slong var, ulong amount, const mpoly_ctx_t mctx)
{
    slong i, N;
    ulong *one;
#if WANT_ASSERT
    ulong mask;
#endif
    TMP_INIT;

    FLINT_ASSERT(Abits <= FLINT_BITS);

#if WANT_ASSERT
    mask = 0;
    for (i = 0; i < FLINT_BITS/Abits; i++)
        mask = (mask << Abits) + (UWORD(1) << (Abits - 1));
#endif

    TMP_START;

    N = mpoly_words_per_exp(Abits, mctx);
    one = (ulong *) TMP_ALLOC(N*sizeof(ulong));

    mpoly_gen_monomial_sp(one, var, Abits, mctx);
    for (i = 0; i < Alength; i++)
    {
        mpoly_monomial_msub(Aexp + N*i, Aexp + N*i, amount, one, N);
        FLINT_ASSERT(!mpoly_monomial_overflows(Aexp + N*i, N, mask));
    }

    TMP_END;
}

void _mpoly_gen_shift_left(ulong * Aexp, flint_bitcnt_t Abits, slong Alength,
                               slong var, ulong amount, const mpoly_ctx_t mctx)
{
    slong i, N;
    ulong *one;
#if WANT_ASSERT
    ulong mask;
#endif
    TMP_INIT;

    FLINT_ASSERT(Abits <= FLINT_BITS);

#if WANT_ASSERT
    mask = 0;
    for (i = 0; i < FLINT_BITS/Abits; i++)
        mask = (mask << Abits) + (UWORD(1) << (Abits - 1));
#endif

    TMP_START;

    N = mpoly_words_per_exp(Abits, mctx);
    one = (ulong *) TMP_ALLOC(N*sizeof(ulong));

    mpoly_gen_monomial_sp(one, var, Abits, mctx);

    for (i = 0; i < Alength; i++)
    {
        mpoly_monomial_madd(Aexp + N*i, Aexp + N*i, amount, one, N);
        FLINT_ASSERT(!mpoly_monomial_overflows(Aexp + N*i, N, mask));
    }

    TMP_END;
}
