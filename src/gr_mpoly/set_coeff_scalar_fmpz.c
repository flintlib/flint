/*
    Copyright (C) 2018 Daniel Schultz
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "gr_mpoly.h"

int gr_mpoly_set_coeff_scalar_fmpz(
    gr_mpoly_t A,
    gr_srcptr c,
    const fmpz * exp,
    const mpoly_ctx_t mctx, gr_ctx_t cctx)
{
    flint_bitcnt_t exp_bits;
    slong i, N, index;
    ulong * cmpmask;
    ulong * packed_exp;
    int exists;
    int status = GR_SUCCESS;
    slong sz = cctx->sizeof_elem;
    TMP_INIT;

    for (i = 0; i < mctx->nvars; i++)
    {
        if (fmpz_sgn(exp + i) < 0)
            return GR_DOMAIN;
    }

    TMP_START;

    exp_bits = mpoly_exp_bits_required_ffmpz(exp, mctx);
    exp_bits = mpoly_fix_bits(exp_bits, mctx);
    gr_mpoly_fit_length_fit_bits(A, A->length, exp_bits, mctx, cctx);

    N = mpoly_words_per_exp(A->bits, mctx);
    cmpmask = (ulong*) TMP_ALLOC(N*sizeof(ulong));
    mpoly_get_cmpmask(cmpmask, N, A->bits, mctx);

    packed_exp = (ulong *) TMP_ALLOC(N*sizeof(ulong));

    mpoly_set_monomial_ffmpz(packed_exp, exp, A->bits, mctx);
    exists = mpoly_monomial_exists(&index, A->exps,
                                  packed_exp, A->length, N, cmpmask);

    if (!exists)
    {
        if (gr_is_zero(c, cctx) != T_TRUE) /* make new term only if coeff is nonzero*/
        {
            gr_mpoly_fit_length(A, A->length + 1, mctx, cctx);

            for (i = A->length; i >= index + 1; i--)
            {
                gr_swap(GR_ENTRY(A->coeffs, i, sz), GR_ENTRY(A->coeffs, i - 1, sz), cctx);
                mpoly_monomial_set(A->exps + N*i, A->exps + N*(i - 1), N);
            }

            status |= gr_set(GR_ENTRY(A->coeffs, index, sz), c, cctx);
            mpoly_monomial_set(A->exps + N*index, packed_exp, N);

            A->length++;
        }
    }
    else if (gr_is_zero(c, cctx) == T_TRUE) /* zero coeff, remove term */
    {
        /* todo: zero out removed coefficient for deallocation? */

        for (i = index; i < A->length - 1; i++)
        {
            gr_swap(GR_ENTRY(A->coeffs, i, sz), GR_ENTRY(A->coeffs, i + 1, sz), cctx);
            mpoly_monomial_set(A->exps + N*i, A->exps + N*(i + 1), N);
        }

        A->length--;
    }
    else /* term with that monomial exists, coeff is nonzero */
    {
        status |= gr_set(GR_ENTRY(A->coeffs, index, sz), c, cctx);
    }

    TMP_END;

    return status;
}

int gr_mpoly_set_coeff_ui_fmpz(
    gr_mpoly_t A,
    ulong c,
    const fmpz * exp,
    const mpoly_ctx_t mctx, gr_ctx_t cctx)
{
    int status;
    gr_ptr t;

    GR_TMP_INIT(t, cctx);
    status = gr_set_ui(t, c, cctx);
    status |= gr_mpoly_set_coeff_scalar_fmpz(A, t, exp, mctx, cctx);
    GR_TMP_CLEAR(t, cctx);

    return status;
}

int gr_mpoly_set_coeff_si_fmpz(
    gr_mpoly_t A,
    slong c,
    const fmpz * exp,
    const mpoly_ctx_t mctx, gr_ctx_t cctx)
{
    int status;
    gr_ptr t;

    GR_TMP_INIT(t, cctx);
    status = gr_set_si(t, c, cctx);
    status |= gr_mpoly_set_coeff_scalar_fmpz(A, t, exp, mctx, cctx);
    GR_TMP_CLEAR(t, cctx);

    return status;
}

int gr_mpoly_set_coeff_fmpz_fmpz(
    gr_mpoly_t A,
    const fmpz_t c,
    const fmpz * exp,
    const mpoly_ctx_t mctx, gr_ctx_t cctx)
{
    int status;
    gr_ptr t;

    GR_TMP_INIT(t, cctx);
    status = gr_set_fmpz(t, c, cctx);
    status |= gr_mpoly_set_coeff_scalar_fmpz(A, t, exp, mctx, cctx);
    GR_TMP_CLEAR(t, cctx);

    return status;
}

int gr_mpoly_set_coeff_fmpq_fmpz(
    gr_mpoly_t A,
    const fmpq_t c,
    const fmpz * exp,
    const mpoly_ctx_t mctx, gr_ctx_t cctx)
{
    int status;
    gr_ptr t;

    GR_TMP_INIT(t, cctx);
    status = gr_set_fmpq(t, c, cctx);
    if (status == GR_SUCCESS)
        status |= gr_mpoly_set_coeff_scalar_fmpz(A, t, exp, mctx, cctx);
    GR_TMP_CLEAR(t, cctx);

    return status;
}
