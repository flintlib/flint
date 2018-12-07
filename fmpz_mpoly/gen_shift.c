/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"

/*
    shift_right (resp. left) divides (resp. multiplies) A in place by x^amount,
    where x is the variable of index var.
    No reallocation on A->bits is performed: The shift right version asserts
    that the division was exact, and the shift left version tries to assert
    that the result fits in the same number of bits.
*/

void _fmpz_mpoly_gen_shift_right(fmpz_mpoly_t A, slong var, ulong amount,
                                                    const fmpz_mpoly_ctx_t ctx)
{
    slong i, N;
    slong off, sh; /* not actually used */
    ulong *one;
#if WANT_ASSERT
    ulong mask;
#endif
    TMP_INIT;

    FLINT_ASSERT(A->bits <= FLINT_BITS);

#if WANT_ASSERT
    mask = 0;
    for (i = 0; i < FLINT_BITS/A->bits; i++)
        mask = (mask << A->bits) + (UWORD(1) << (A->bits - 1));
#endif

    TMP_START;

    N = mpoly_words_per_exp(A->bits, ctx->minfo);
    one = (ulong *) TMP_ALLOC(N*sizeof(ulong));

    mpoly_gen_oneexp_offset_shift(one, &off, &sh, var, N, A->bits, ctx->minfo);
    for (i = 0; i < A->length; i++)
    {
        mpoly_monomial_msub(A->exps + N*i, A->exps + N*i, amount, one, N);
        FLINT_ASSERT(!mpoly_monomial_overflows(A->exps + N*i, N, mask));
    }

    TMP_END;
}

void _fmpz_mpoly_gen_shift_left(fmpz_mpoly_t A, slong var, ulong amount,
                                                    const fmpz_mpoly_ctx_t ctx)
{
    slong i, N;
    slong off, sh; /* not actually used */
    ulong *one;
#if WANT_ASSERT
    ulong mask;
#endif
    TMP_INIT;

    FLINT_ASSERT(A->bits <= FLINT_BITS);

#if WANT_ASSERT
    mask = 0;
    for (i = 0; i < FLINT_BITS/A->bits; i++)
        mask = (mask << A->bits) + (UWORD(1) << (A->bits - 1));
#endif

    TMP_START;

    N = mpoly_words_per_exp(A->bits, ctx->minfo);
    one = (ulong *) TMP_ALLOC(N*sizeof(ulong));

    mpoly_gen_oneexp_offset_shift(one, &off, &sh, var, N, A->bits, ctx->minfo);

    for (i = 0; i < A->length; i++)
    {
        mpoly_monomial_madd(A->exps + N*i, A->exps + N*i, amount, one, N);
        FLINT_ASSERT(!mpoly_monomial_overflows(A->exps + N*i, N, mask));
    }

    TMP_END;
}
