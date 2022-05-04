/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"

/*
    When selecting expoenent ranges for the parallel division A/B,
    we can also detect if an exact division is impossible.
    The return is non zero in this case.

    This function is used by {fmpz|nmod}_mpoly_divides_heap_threaded, and
    must be placed in the fmpz_mpoly module becuase it uses an fmpz_mpoly_t
    to accumulate and sort the exponents.
*/
int mpoly_divides_select_exps(fmpz_mpoly_t S, fmpz_mpoly_ctx_t zctx,
                             slong nworkers, ulong * Aexp, slong Alen,
                                    ulong * Bexp, slong Blen, flint_bitcnt_t bits)
{
    int failure;
    ulong mask;
    ulong * Sexp;
    slong Slen;
    fmpz * Scoeff;
    slong nA = 30 + 8*nworkers;     /* number of division of A */
    slong nB = (1 + nworkers)/2;    /* number of division of B */
    slong tot;
    ulong * T0, * T1;
    slong i, j, N;
    TMP_INIT;

    TMP_START;

    N = mpoly_words_per_exp(bits, zctx->minfo);

    mask = bits <= FLINT_BITS ? mpoly_overflow_mask_sp(bits) : 0;

    FLINT_ASSERT(Alen > 0);
    FLINT_ASSERT(Blen > 0);

    tot = 16 + nA + 2*nB;
    fmpz_mpoly_fit_bits(S, bits, zctx);
    S->bits = bits;
    fmpz_mpoly_fit_length(S, tot, zctx);
    Sexp = S->exps;
    Scoeff = S->coeffs;
    Slen = 0;

    /* get a (non-strict) upper bound on the exponents of A */
    mpoly_monomial_set(Sexp + N*Slen, Aexp + N*0, N);
    fmpz_one(Scoeff + Slen);
    Slen++;

    for (i = 1; i < nA; i++)
    {
        double a = 1.0;
        double b = 0.2;
        double d = (double)(i) / (double)(nA);
        /*
            set d to p(d) where p is the cubic satisfying
            p(0) = 0, p'(0) = a
            p(1) = 1, p'(1) = b
        */
        FLINT_ASSERT(a >= 0);
        FLINT_ASSERT(b >= 0);
        d = d*(1 + (1 - d)*((2 - a - b)*d - (1 - a)));

        /* choose exponent at relative position d in A */
        j = d * Alen;
        j = FLINT_MAX(j, WORD(0));
        j = FLINT_MIN(j, Alen - 1);
        mpoly_monomial_set(Sexp + N*Slen, Aexp + N*j, N);
        fmpz_one(Scoeff + Slen);
        Slen++;
    }
    _fmpz_mpoly_set_length(S, Slen, zctx);

    T0 = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    T1 = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    mpoly_monomial_sub_mp(T0, Aexp + N*0, Bexp + N*0, N);
    mpoly_monomial_sub_mp(T1, Aexp + N*(Alen - 1), Bexp + N*(Blen - 1), N);
    if (bits <= FLINT_BITS)
    {
        if (   mpoly_monomial_overflows(T0, N, mask)
            || mpoly_monomial_overflows(T1, N, mask))
        {
            failure = 1;
            goto cleanup;
        }
    }
    else
    {
        if (   mpoly_monomial_overflows_mp(T0, N, bits)
            || mpoly_monomial_overflows_mp(T1, N, bits))
        {
            failure = 1;
            goto cleanup;
        }
    }


    for (i = 1; i < nB; i++)
    {
        double d = (double)(i) / (double)(nB);
        /* choose exponent at relative position d in B */
        j = d * Blen;
        j = FLINT_MAX(j, WORD(0));
        j = FLINT_MIN(j, Blen - 1);

        mpoly_monomial_sub_mp(Sexp + N*Slen, Aexp + N*0, Bexp + N*0, N);
        mpoly_monomial_add_mp(Sexp + N*Slen, Sexp + N*Slen, Bexp + N*j, N);
        fmpz_one(Scoeff + Slen);
        if (bits <= FLINT_BITS)
            Slen += !(mpoly_monomial_overflows(Sexp + N*Slen, N, mask));
        else
            Slen += !(mpoly_monomial_overflows_mp(Sexp + N*Slen, N, bits));

        mpoly_monomial_sub_mp(Sexp + N*Slen, Aexp + N*(Alen - 1), Bexp + N*(Blen - 1), N);
        mpoly_monomial_add_mp(Sexp + N*Slen, Sexp + N*Slen, Bexp + N*j, N);
        fmpz_one(Scoeff + Slen);
        if (bits <= FLINT_BITS)
            Slen += !(mpoly_monomial_overflows(Sexp + N*Slen, N, mask));
        else
            Slen += !(mpoly_monomial_overflows_mp(Sexp + N*Slen, N, bits));
    }

    /* get a lower bound on the exponents of A */
    mpoly_monomial_zero(Sexp + N*Slen, N);
    fmpz_one(Scoeff + Slen);
    Slen++;

    /* done constructing S */
    FLINT_ASSERT(Slen < tot);

    _fmpz_mpoly_set_length(S, Slen, zctx);

    fmpz_mpoly_sort_terms(S, zctx);
    fmpz_mpoly_combine_like_terms(S, zctx);

    failure = 0;

cleanup:
    TMP_END;
    return failure;
}
