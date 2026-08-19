/*
    Copyright (C) 2026 Fredrik Johansson
    Developed using Claude Fable 5

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gmpcompat.h"
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "fft_small.h"
#include "gr.h"

static nn_srcptr
_coeff_limbs(const fmpz * f, ulong * scratch, slong * n, int * sgn)
{
    fmpz c = *f;

    if (!COEFF_IS_MPZ(c))
    {
        *scratch = FLINT_ABS(c);
        *n = 1;
        *sgn = c < 0;
        return scratch;
    }
    else
    {
        mpz_srcptr m = COEFF_TO_PTR(c);

        *n = FLINT_ABS(m->_mp_size);
        *sgn = m->_mp_size < 0;
        return m->_mp_d;
    }
}

static int
_set_coeff(gr_ptr elem, const fmpz * f, gr_ctx_t ctx)
{
    ulong scratch;
    slong n;
    int sgn;
    nn_srcptr p = _coeff_limbs(f, &scratch, &n, &sgn);

    return gr_transformed_mpn_set(elem, p, n, sgn, ctx);
}

/* the accumulated element, converted out directly into an fmpz */
/* converts the accumulator out destructively -- the caller re-initializes
   it (borrowed storage, so that is free) before the next output */
static int
_get_coeff(fmpz * f, gr_ptr elem, gr_ctx_t ctx)
{
    slong w = gr_transformed_mpn_get_limbs(ctx, elem);
    slong zl;
    int sg, status;
    mpz_ptr m;

    if (w <= 0)
        w = 1;

    m = _fmpz_promote(f);
    status = gr_transformed_mpn_get_destructive(FLINT_MPZ_REALLOC(m, w), w,
                                                &zl, &sg, elem, ctx);
    if (status != GR_SUCCESS)
    {
        _fmpz_demote(f);
        fmpz_zero(f);
        return status;
    }

    m->_mp_size = sg ? -zl : zl;
    _fmpz_demote_val(f);
    return GR_SUCCESS;
}

int
_fmpz_poly_mulmid_classical_fft_small(fmpz * res,
    const fmpz * poly1, slong len1, const fmpz * poly2, slong len2,
    slong nlo, slong nhi)
{
#if !FLINT_HAVE_FFT_SMALL
    return 0;
#else
    gr_ctx_t ctx;
    gr_ptr B, X, acc, acc2, E;
    slong * xidx;                 /* which poly1 index each X slot holds */
    slong k = len2, i, j, abits = 0, bbits = 0, sz;
    int share = (poly1 == poly2);
    int signed_needed;
    int status = GR_SUCCESS;

    if (len2 < 1 || len1 < 1 || nlo >= nhi)
        return 0;

    /* the product window is symmetric in the operands; keep the
       shorter operand as the transformed-and-kept side */
    if (len1 < len2)
    {
        const fmpz * tp = poly1;
        slong tl = len1;
        poly1 = poly2;
        len1 = len2;
        poly2 = tp;
        len2 = tl;
        k = len2;
    }

    /* _fmpz_vec_max_bits encodes "any negative coefficient" in its
       sign, so the bounds and the signedness scan are one pass; an
       all-nonnegative product takes the cheaper unsigned representation
       (as in the fmpz_mat driver) */
    {
        slong amax = _fmpz_vec_max_bits(poly1, len1);
        slong bmax = share ? amax : _fmpz_vec_max_bits(poly2, len2);
        abits = FLINT_ABS(amax);
        bbits = FLINT_ABS(bmax);
        signed_needed = (amax < 0 || bmax < 0);
    }

    if (abits == 0 || bbits == 0)
    {
        /* a zero operand: the window of the product is zero */
        _fmpz_vec_zero(res, nhi - nlo);
        return 1;
    }

    /* Decline if too sparse; to be moved to the algorithm selection in
       _fmpz_poly_mulmid when this  function has been integrated.
    {
        slong wa = _fmpz_vec_weight_bits(poly1, len1);
        slong wb = share ? wa : _fmpz_vec_weight_bits(poly2, len2);
        double density = (double) (wa + wb) / ((double) len1 * abits + (double) len2 * bbits);
        double need = (k <= 2) ? 0.90
                    : (k <= 4) ? 0.62
                    : (k <= 8) ? 0.50
                    : 0.30;
        if (density < need)
            return 0;
    } */

    /* one product's magnitude; the k-term accumulation and the sign
       are the context's own provisioning */
    if (gr_ctx_init_transformed_mpn(ctx,
            abits + bbits, k, signed_needed, 2 * k + 2,
            GR_TRANSFORMED_MPN_ALLOC_FIT_BUFFER) != GR_SUCCESS)
        return 0;

    sz = ctx->sizeof_elem;

    /* k kept transforms of poly2, k rolling slots for poly1 indices >= the
       kept range (all indices when not sharing); the accumulators live on
       borrowed slabs so they can be converted out destructively and
       re-initialized for free */
    GR_TMP_INIT_VEC(E, 2 * k + 2, ctx);
    B = E;
    X = GR_ENTRY(E, k, sz);
    acc = GR_ENTRY(E, 2 * k, sz);
    acc2 = GR_ENTRY(E, 2 * k + 1, sz);

    xidx = flint_malloc(k * sizeof(slong));
    for (i = 0; i < k; i++)
        xidx[i] = -1;

    /* the kept transforms; when sharing they serve both roles */
    for (j = 0; j < k && status == GR_SUCCESS; j++)
        status |= _set_coeff(GR_ENTRY(B, j, sz),
                             (share ? poly1 : poly2) + j, ctx);

#define A_ELEM(t, out) \
    do { \
        if (share && (t) < k) \
            (out) = GR_ENTRY(B, (t), sz); \
        else \
        { \
            slong _s = share ? ((t) - k) % k : (t) % k; \
            gr_ptr _e = GR_ENTRY(X, _s, sz); \
            if (xidx[_s] != (t)) \
            { \
                status |= _set_coeff(_e, poly1 + (t), ctx); \
                xidx[_s] = (t); \
            } \
            (out) = _e; \
        } \
    } while (0)

    for (i = nlo; i < nhi && status == GR_SUCCESS; i++)
    {
        slong jlo = FLINT_MAX(0, i - len1 + 1);
        slong jhi = FLINT_MIN(i, k - 1);
        gr_ptr ea, eb;

        if (!share)
        {
            for (j = jlo; j <= jhi; j++)
            {
                A_ELEM(i - j, ea);
                eb = GR_ENTRY(B, j, sz);
                if (status != GR_SUCCESS)
                    break;
                status |= gr_addmul(acc, ea, eb, ctx);
            }
        }
        else
        {
            /* pair j with i - j: within [jlo, jhi] both orders of each
               product appear, so multiply once and add the pair sum in
               twice; the diagonal term (i even) is a plain square, and
               terms whose partner falls outside the range stay single */
            for (j = jlo; j <= jhi; j++)
            {
                slong jp = i - j;         /* the partner index */
                int paired = (jp >= jlo && jp <= jhi);

                if (paired && j > jp)
                    continue;             /* counted at the partner */

                A_ELEM(i - j, ea);
                A_ELEM(j, eb);
                if (status != GR_SUCCESS)
                    break;

                if (paired && j < jp)
                {
                    status |= gr_addmul(acc2, ea, eb, ctx);
                }
                else
                {
                    status |= gr_addmul(acc, ea, eb, ctx);
                }
            }

            /* fold the paired sum in twice; when nothing paired,
            acc2 is the ring's zero and the additions are
            short-circuited bookkeeping */
            status |= gr_add(acc, acc, acc2, ctx);
            status |= gr_add(acc, acc, acc2, ctx);
        }

        if (status != GR_SUCCESS)
            break;

        status |= _get_coeff(res + (i - nlo), acc, ctx);
            /* the destructive get consumed acc: bring it back for the
               next output; under the fit_buffer strategy this is a
               slab-pool round trip */
        gr_clear(acc, ctx);
        gr_init(acc, ctx);
        gr_clear(acc2, ctx);
        gr_init(acc2, ctx);
    }

#undef A_ELEM

    GR_TMP_CLEAR_VEC(E, 2 * k + 2, ctx);
    flint_free(xidx);
    gr_ctx_clear(ctx);

    return status == GR_SUCCESS;
#endif
}

int
fmpz_poly_mulmid_classical_fft_small(fmpz_poly_t res,
    const fmpz_poly_t poly1, const fmpz_poly_t poly2, slong nlo,
    slong nhi)
{
    slong len1 = poly1->length, len2 = poly2->length;
    int ok;

    if (len1 == 0 || len2 == 0 || nlo >= nhi)
    {
        fmpz_poly_zero(res);
        return 1;
    }

    if (res == poly1 || res == poly2)
    {
        fmpz_poly_t t;
        fmpz_poly_init2(t, nhi - nlo);
        ok = fmpz_poly_mulmid_classical_fft_small(t, poly1, poly2,
                                                  nlo, nhi);
        if (ok)
            fmpz_poly_swap(res, t);
        fmpz_poly_clear(t);
        return ok;
    }

    fmpz_poly_fit_length(res, nhi - nlo);
    ok = _fmpz_poly_mulmid_classical_fft_small(res->coeffs,
            poly1->coeffs, len1, poly2->coeffs, len2, nlo, nhi);
    if (ok)
    {
        _fmpz_poly_set_length(res, nhi - nlo);
        _fmpz_poly_normalise(res);
    }
    return ok;
}
