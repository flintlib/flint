/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gmpcompat.h"
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

/*
    Shared driver: lo_limbs = 0 gives the exact product; lo_limbs > 0
    returns, for each entry, the limbs of the magnitude above limb
    lo_limbs (with the entry's sign), which is what a working-precision
    matrix product such as arb_mat_mul needs -- the discarded low tail
    perturbs each entry by at most one unit in its lowest returned limb.

    Memory strategy. The driver budgets its transform storage against
    flint_fft_small_max_transformed_ring_size and picks the cheapest
    blocking that fits:

    - everything resident (the classical amortization: every entry
      transformed once, every output converted once); squaring (B == A)
      keeps a single transform pool here, halving the transforms and the
      memory;
    - one side resident, the other streamed in row (column) blocks --
      the windowed strategy for rectangular products: still no entry is
      transformed twice;
    - both sides blocked: the streamed inner side is re-transformed once
      per outer block, trading forward transforms (the cheap phase) for
      fitting in memory;
    - as a last resort the inner dimension itself is blocked and the
      output entries accumulated across the pieces on the fmpz side;
      this multiplies the conversions out and is unavailable for the
      truncated variant, whose one-ulp contract does not survive summing
      truncated pieces.

    Whether the transformed representation is worth using at all (entry
    sizes, dimensions) is the caller's tuning decision; the driver
    accepts any input it can compute.
*/

static int
_set_entry(gr_ptr elem, const fmpz * e, gr_ctx_t tctx)
{
    fmpz c = *e;

    if (!COEFF_IS_MPZ(c))
    {
        ulong v = FLINT_ABS(c);
        return gr_transformed_mpn_set(elem, &v, c != 0, c < 0, tctx)
                == GR_SUCCESS;
    }
    else
    {
        mpz_srcptr m = COEFF_TO_PTR(c);
        return gr_transformed_mpn_set(elem, m->_mp_d,
                FLINT_ABS(m->_mp_size), m->_mp_size < 0, tctx)
                == GR_SUCCESS;
    }
}

/* convert the accumulator out into e (add_into: accumulate); consumes
   acc, which the caller re-initializes */
static int
_export_entry(fmpz * e, gr_ptr acc, slong lo_limbs, int add_into,
              nn_ptr t, slong tn_max, gr_ctx_t tctx)
{
    slong need, zn;
    int sg, ok;

    need = (lo_limbs == 0)
            ? gr_transformed_mpn_get_limbs(tctx, acc)
            : gr_transformed_mpn_get_limbs_trunc(tctx, acc, lo_limbs);
    /* the buffer was sized from gr_transformed_mpn_get_limbs_bound, of
       which need is a per-element instance */
    FLINT_ASSERT(need <= tn_max);
    (void) tn_max;
    if (lo_limbs == 0)
        ok = gr_transformed_mpn_get_destructive(t, need, &zn, &sg, acc,
                tctx) == GR_SUCCESS;
    else
        ok = gr_transformed_mpn_get_trunc_destructive(t, need, &zn, &sg,
                lo_limbs, acc, tctx) == GR_SUCCESS;
    if (!ok)
        return 0;

    if (!add_into)
    {
        if (zn == 0)
            fmpz_zero(e);
        else
        {
            fmpz_set_ui_array(e, t, zn);
            if (sg)
                fmpz_neg(e, e);
        }
    }
    else if (zn != 0)
    {
        fmpz_t p;
        fmpz_init(p);
        fmpz_set_ui_array(p, t, zn);
        if (sg)
            fmpz_sub(e, e, p);
        else
            fmpz_add(e, e, p);
        fmpz_clear(p);
    }
    return 1;
}

static int
_fmpz_mat_mul_fft_small(fmpz_mat_t C, const fmpz_mat_t A, const fmpz_mat_t B,
                        slong lo_limbs)
{
#if FLINT_HAVE_FFT_SMALL
    slong ar = fmpz_mat_nrows(A);
    slong k = fmpz_mat_ncols(A);
    slong bc = fmpz_mat_ncols(B);
    slong abits, bbits, bits_bound;
    int signed_needed, share;
    slong i, j, l, I, J, L;
    gr_ctx_t tctx;
    gr_ptr EA, EB, acc;
    nn_ptr t;
    slong tn_max, budget, mb, nb, kb;
    int ok = 1;
#define EA_(i) GR_ENTRY(EA, i, tctx->sizeof_elem)
#define EB_(i) GR_ENTRY(EB, i, tctx->sizeof_elem)

    if (ar == 0 || bc == 0 || k == 0)
        return 0;
    if (C == A || C == B)
        return 0;

    /* fmpz_mat_max_bits encodes "any negative entry" in its sign, so
       the nonnegativity scan comes free with the bounds (and the calls
       are made once: FLINT_ABS would evaluate its argument twice) */
    {
        slong amax = fmpz_mat_max_bits(A);
        slong bmax = fmpz_mat_max_bits(B);

        abits = FLINT_ABS(amax);
        bbits = FLINT_ABS(bmax);
        signed_needed = (amax < 0 || bmax < 0);
    }
    if (abits == 0 || bbits == 0)
        return 0;

    share = (A == B);

    /* the bound describes one product's magnitude: the context
       provisions the k-fold accumulation from terms_bound and the sign
       from is_signed itself (an earlier revision added
       FLINT_BIT_COUNT(k) + 2 here, duplicating headroom the
       representation owns) */
    bits_bound = abits + bbits;

    /* the exact live count depends on the blocking decision below,
       which in turn needs the context's element size: create the
       context with plain allocation and switch to the pooled
       strategy once the count is known */
    if (gr_ctx_init_transformed_mpn(tctx, bits_bound, k, signed_needed,
            0, GR_TRANSFORMED_MPN_ALLOC_MALLOC)
            != GR_SUCCESS)
        return 0;

    /* elements of transform storage the budget admits */
    {
        ulong per = n_round_up(gr_transformed_mpn_sizeof_data(tctx),
                               FLINT_FFT_SMALL_ALIGNMENT);
        ulong nb_ = flint_fft_small_max_transformed_ring_size / per;
        budget = (nb_ > (ulong) WORD_MAX) ? WORD_MAX : (slong) nb_;
    }
    if (budget < 3)
    {
        gr_ctx_clear(tctx);
        return 0;
    }

    /* blocking decision: mb rows of A and nb columns of B resident at a
       time, over the full inner dimension (kb = k) when possible */
    kb = k;
    if (share && ar * k + 1 <= budget)
    {
        mb = ar; nb = bc;               /* single resident pool */
    }
    else if (!share && ar * k + k * bc + 1 <= budget)
    {
        mb = ar; nb = bc;               /* everything resident */
    }
    else if (k * bc + k + 1 <= budget)
    {
        nb = bc;                        /* B resident, stream A rows */
        mb = FLINT_MIN(ar, (budget - 1 - k * bc) / k);
        share = 0;
    }
    else if (ar * k + k + 1 <= budget)
    {
        mb = ar;                        /* A resident, stream B columns */
        nb = FLINT_MIN(bc, (budget - 1 - ar * k) / k);
        share = 0;
    }
    else if (2 * k + 1 <= budget)
    {
        /* both sides blocked; the inner (B) side is re-transformed once
           per outer block */
        mb = FLINT_MIN(ar, (budget - 1) / (2 * k));
        nb = FLINT_MIN(bc, (budget - 1) / (2 * k));
        share = 0;
    }
    else
    {
        /* the inner dimension itself must be blocked; entries accumulate
           across the pieces, which multiplies the conversions out and
           breaks the truncated contract */
        if (lo_limbs > 0)
        {
            gr_ctx_clear(tctx);
            return 0;
        }
        mb = 1; nb = 1;
        kb = (budget - 1) / 2;
        share = 0;
    }

    /* every element the blocking keeps live at once: the resident A
       block, the resident B block unless shared, and the accumulator;
       a refusal (a nested reservation) simply leaves plain
       allocation in place */
    (void) gr_transformed_mpn_use_fit_buffer(tctx,
        mb * kb + (share ? 0 : kb * nb) + 1);

    EA = flint_malloc((mb * kb) * tctx->sizeof_elem);
    for (i = 0; i < mb * kb; i++)
        gr_init(EA_(i), tctx);
    if (share)
        EB = EA;
    else
    {
        EB = flint_malloc((kb * nb) * tctx->sizeof_elem);
        for (i = 0; i < kb * nb; i++)
            gr_init(EB_(i), tctx);
    }
    acc = flint_malloc(tctx->sizeof_elem);
    gr_init(acc, tctx);

    /* conversion staging sized by the ring's own bound rather than a
       reconstruction of its limb requirements here */
    tn_max = gr_transformed_mpn_get_limbs_bound(tctx);
    t = flint_malloc(tn_max * sizeof(ulong));

    for (L = 0; ok && L < k; L += kb)
    {
        slong kl = FLINT_MIN(kb, k - L);

        for (I = 0; ok && I < ar; I += mb)
        {
            slong ml = FLINT_MIN(mb, ar - I);

            /* the A block, transformed once per (L, I) */
            for (i = 0; ok && i < ml; i++)
                for (l = 0; ok && l < kl; l++)
                    ok = _set_entry(EA_(i * kb + l),
                            fmpz_mat_entry(A, I + i, L + l), tctx);

            for (J = 0; ok && J < bc; J += nb)
            {
                slong nl = FLINT_MIN(nb, bc - J);

                /* the B block; when both sides fit this runs once,
                   when only A is blocked it runs once per J, and in the
                   doubly blocked regime once per (I, J) */
                /* EB already holds block (L, J) whenever B is fully
                   resident in the column direction and I has advanced;
                   otherwise the J (or J and I) traversal overwrote it */
                if (!share && (I == 0 || nb < bc))
                    for (l = 0; ok && l < kl; l++)
                        for (j = 0; ok && j < nl; j++)
                            ok = _set_entry(EB_(l * nb + j),
                                    fmpz_mat_entry(B, L + l, J + j), tctx);

                for (i = 0; ok && i < ml; i++)
                    for (j = 0; ok && j < nl; j++)
                    {
                        gr_srcptr b0 = share
                            ? EA_((L + 0) * kb + (J + j))
                            : EB_(0 * nb + j);

                        ok = gr_mul(acc, EA_(i * kb + 0), b0, tctx)
                                == GR_SUCCESS;
                        for (l = 1; ok && l < kl; l++)
                        {
                            gr_srcptr bl = share
                                ? EA_((L + l) * kb + (J + j))
                                : EB_(l * nb + j);
                            ok = gr_addmul(acc, EA_(i * kb + l), bl, tctx)
                                    == GR_SUCCESS;
                        }

                        if (ok)
                            ok = _export_entry(fmpz_mat_entry(C, I + i, J + j),
                                    acc, lo_limbs, L > 0, t, tn_max, tctx);

                        gr_clear(acc, tctx);
                        gr_init(acc, tctx);
                    }
            }
        }
    }

    for (i = 0; i < mb * kb; i++)
        gr_clear(EA_(i), tctx);
    flint_free(EA);
    if (!share)
    {
        for (i = 0; i < kb * nb; i++)
            gr_clear(EB_(i), tctx);
        flint_free(EB);
    }
    gr_clear(acc, tctx);
    flint_free(acc);
    flint_free(t);
    gr_ctx_clear(tctx);
#undef EA_
#undef EB_
    return ok;
#else
    (void) C; (void) A; (void) B; (void) lo_limbs;
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
