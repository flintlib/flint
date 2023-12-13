/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mpoly.h"

void _fmpz_mod_mpoly_set_coeff_fmpz_fmpz(
    fmpz_mod_mpoly_t A,
    const fmpz_t c,
    const fmpz * exp,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    flint_bitcnt_t exp_bits;
    slong i, N, index;
    ulong * cmpmask;
    ulong * packed_exp;
    int exists;
    fmpz_t C;
    TMP_INIT;

    TMP_START;

    fmpz_init(C);
    fmpz_mod_set_fmpz(C, c, ctx->ffinfo);

    exp_bits = mpoly_exp_bits_required_ffmpz(exp, ctx->minfo);
    exp_bits = mpoly_fix_bits(exp_bits, ctx->minfo);
    fmpz_mod_mpoly_fit_length_fit_bits(A, A->length, exp_bits, ctx);

    N = mpoly_words_per_exp(A->bits, ctx->minfo);
    cmpmask = (ulong*) TMP_ALLOC(N*sizeof(ulong));
    mpoly_get_cmpmask(cmpmask, N, A->bits, ctx->minfo);

    packed_exp = (ulong *) TMP_ALLOC(N*sizeof(ulong));

    mpoly_set_monomial_ffmpz(packed_exp, exp, A->bits, ctx->minfo);
    exists = mpoly_monomial_exists(&index, A->exps,
                                  packed_exp, A->length, N, cmpmask);

    if (!exists)
    {
        if (!fmpz_is_zero(C)) /* make new term only if coeff is nonzero*/
        {
            fmpz_mod_mpoly_fit_length(A, A->length + 1, ctx);

            for (i = A->length; i >= index + 1; i--)
            {
                fmpz_set(A->coeffs + i, A->coeffs + i - 1);
                mpoly_monomial_set(A->exps + N*i, A->exps + N*(i - 1), N);
            }

            fmpz_swap(A->coeffs + index, C);
            mpoly_monomial_set(A->exps + N*index, packed_exp, N);

            A->length++;
        }
    }
    else if (fmpz_is_zero(C)) /* zero coeff, remove term */
    {
        for (i = index; i < A->length - 1; i++)
        {
            fmpz_set(A->coeffs + i, A->coeffs + i + 1);
            mpoly_monomial_set(A->exps + N*i, A->exps + N*(i + 1), N);
        }

        A->length--;
    }
    else /* term with that monomial exists, coeff is nonzero */
    {
        fmpz_swap(A->coeffs + index, C);
    }

    fmpz_clear(C);

    TMP_END;
}


void fmpz_mod_mpoly_set_coeff_fmpz_fmpz(
    fmpz_mod_mpoly_t A,
    const fmpz_t c,
    fmpz * const * exp,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong i, nvars = ctx->minfo->nvars;
    fmpz * newexp;
    TMP_INIT;

    TMP_START;

    newexp = TMP_ALLOC(nvars * sizeof(fmpz));
    for (i = 0; i < nvars; i++)
        newexp[i] = *exp[i];

    /* GCC really wants to complain about this one */
#if defined(__GNUC__) && !defined(__clang__)
# pragma GCC diagnostic push
# pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#endif
    _fmpz_mod_mpoly_set_coeff_fmpz_fmpz(A, c, newexp, ctx);
#if defined(__GNUC__) && !defined(__clang__)
# pragma GCC diagnostic pop
#endif

    TMP_END;
}

void fmpz_mod_mpoly_set_coeff_ui_fmpz(fmpz_mod_mpoly_t poly,
                   ulong c, fmpz * const * exp, const fmpz_mod_mpoly_ctx_t ctx)
{
    fmpz_t C;
    fmpz_init_set_ui(C, c);
    fmpz_mod_mpoly_set_coeff_fmpz_fmpz(poly, C, exp, ctx);
    fmpz_clear(C);
}

void fmpz_mod_mpoly_set_coeff_si_fmpz(fmpz_mod_mpoly_t poly,
                   slong c, fmpz * const * exp, const fmpz_mod_mpoly_ctx_t ctx)
{
    fmpz_t C;
    fmpz_init_set_si(C, c);
    fmpz_mod_mpoly_set_coeff_fmpz_fmpz(poly, C, exp, ctx);
    fmpz_clear(C);
}

