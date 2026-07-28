/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_mat.h"
#include "gr.h"
#if FLINT_HAVE_FFT_SMALL
# include "fft_small.h"
#endif

/*
    Integer matrix product over a transformed-integer context: each of the
    ar*k + k*bc entries is transformed once and reused across the ar*bc*k
    pointwise products, with entry signs handled by the context's sign
    machinery (pointwise addmul/submul switching, signed conversion out).
    Returns 0 when the transformed representation is unavailable or
    judged unprofitable (small entries); C may then be partially written,
    but A and B are untouched, so a caller falls back to fmpz_mat_mul.

    Currently single-threaded at the entry level (the context itself is
    thread safe, so the two-phase pool pattern of the polynomial matrix
    drivers transfers directly; the internal transform and reconstruction
    stages already thread over large entries).
*/

#define FMPZ_MAT_MUL_FFT_SMALL_MIN_BITS 16384

/*
    Shared driver: lo_limbs = 0 gives the exact product; lo_limbs > 0
    returns, for each entry, the limbs of the magnitude above limb
    lo_limbs (with the entry's sign), which is what a working-precision
    matrix product such as arb_mat_mul needs -- the discarded low tail
    perturbs each entry by at most one unit in its lowest returned limb.
*/
static int
_fmpz_mat_mul_fft_small(fmpz_mat_t C, const fmpz_mat_t A, const fmpz_mat_t B,
                        slong lo_limbs)
{
#if FLINT_HAVE_FFT_SMALL
    slong ar = fmpz_mat_nrows(A);
    slong k = fmpz_mat_ncols(A);
    slong bc = fmpz_mat_ncols(B);
    slong abits, bbits, bits_bound;
    slong i, j, l;
    gr_ctx_t tctx;
    gr_ptr E, acc;
    nn_ptr t;
    slong tn_max;
    int ok = 1;
#define E_(i) GR_ENTRY(E, i, tctx->sizeof_elem)

    if (ar == 0 || bc == 0 || k == 0)
        return 0;
    if (C == A || C == B)
        return 0;

    abits = FLINT_ABS(fmpz_mat_max_bits(A));
    bbits = FLINT_ABS(fmpz_mat_max_bits(B));
    if (abits == 0 || bbits == 0)
        return 0;

    if (abits + bbits < FMPZ_MAT_MUL_FFT_SMALL_MIN_BITS)
        return 0;

    bits_bound = abits + bbits + FLINT_BIT_COUNT((ulong) k) + 2;

    if (gr_ctx_init_transformed_mpn(tctx, bits_bound, k, 1) != GR_SUCCESS)
        return 0;

    E = flint_malloc((ar * k + k * bc) * tctx->sizeof_elem);
    for (i = 0; i < ar * k + k * bc; i++)
        gr_init(E_(i), tctx);
    acc = gr_heap_init(tctx);

    tn_max = (abits + bbits + 128) / FLINT_BITS + 4 + k;
    t = flint_malloc(tn_max * sizeof(ulong));

    /* conversions in */
    for (i = 0; ok && i < ar; i++)
        for (l = 0; ok && l < k; l++)
        {
            const fmpz * e = fmpz_mat_entry(A, i, l);
            slong en = fmpz_size(e);
            fmpz_t ea; fmpz_init(ea); fmpz_abs(ea, e);
            fmpz_get_ui_array(t, FLINT_MAX(en, 1), ea);
            ok = gr_transformed_mpn_set(E_(i * k + l), t, en,
                    fmpz_sgn(e) < 0, tctx) == GR_SUCCESS;
            fmpz_clear(ea);
        }
    for (l = 0; ok && l < k; l++)
        for (j = 0; ok && j < bc; j++)
        {
            const fmpz * e = fmpz_mat_entry(B, l, j);
            slong en = fmpz_size(e);
            fmpz_t ea; fmpz_init(ea); fmpz_abs(ea, e);
            fmpz_get_ui_array(t, FLINT_MAX(en, 1), ea);
            ok = gr_transformed_mpn_set(E_(ar * k + l * bc + j), t, en,
                    fmpz_sgn(e) < 0, tctx) == GR_SUCCESS;
            fmpz_clear(ea);
        }

    /* accumulate and convert out */
    for (i = 0; ok && i < ar; i++)
        for (j = 0; ok && j < bc; j++)
        {
            slong need, zn; int sg;

            ok = gr_mul(acc, E_(i * k), E_(ar * k + j), tctx) == GR_SUCCESS;
            for (l = 1; ok && l < k; l++)
                ok = gr_addmul(acc, E_(i * k + l),
                        E_(ar * k + l * bc + j), tctx) == GR_SUCCESS;

            need = (lo_limbs == 0)
                    ? gr_transformed_mpn_get_limbs(tctx, acc)
                    : gr_transformed_mpn_get_limbs_trunc(tctx, acc, lo_limbs);
            if (need > tn_max)
            {
                t = flint_realloc(t, need * sizeof(ulong));
                tn_max = need;
            }
            if (lo_limbs == 0)
                ok = ok && gr_transformed_mpn_get(t, need, &zn, &sg, acc,
                        tctx) == GR_SUCCESS;
            else
                ok = ok && gr_transformed_mpn_get_trunc(t, need, &zn, &sg,
                        lo_limbs, acc, tctx) == GR_SUCCESS;
            if (ok)
            {
                fmpz * e = fmpz_mat_entry(C, i, j);
                if (zn == 0)
                    fmpz_zero(e);
                else
                {
                    fmpz_set_ui_array(e, t, zn);
                    if (sg)
                        fmpz_neg(e, e);
                }
            }
        }

    for (i = 0; i < ar * k + k * bc; i++)
        gr_clear(E_(i), tctx);
    flint_free(E);
    gr_heap_clear(acc, tctx);
    flint_free(t);
    gr_ctx_clear(tctx);
#undef E_
    return ok;
#else
    return 0;
#endif
}

int
fmpz_mat_mul_fft_small(fmpz_mat_t C, const fmpz_mat_t A, const fmpz_mat_t B)
{
    return _fmpz_mat_mul_fft_small(C, A, B, 0);
}

int
fmpz_mat_mul_fft_small_trunc(fmpz_mat_t C, const fmpz_mat_t A,
                             const fmpz_mat_t B, slong lo_limbs)
{
    if (lo_limbs < 0)
        return 0;
    return _fmpz_mat_mul_fft_small(C, A, B, lo_limbs);
}
