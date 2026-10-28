/*
    From GMP's mpn/generic/gcd.c:
    Copyright 1991, 1993, 1994, 2000-2005, 2008, 2010, 2012, 2019 Free
    Software Foundation, Inc.

    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"
#include "impl.h"

/*
    Uses the HGCD operation described in

        N. Möller, On Schönhage's algorithm and subquadratic integer gcd
        computation, Math. Comp. 77 (2008), 589-607.

    to reduce inputs until they are of size below FLINT_MPN_GCD_DC_CUTOFF,
    and then uses Lehmer's algorithm.

    Some reasonable choices are n / 2 (same as in hgcd), and p = (n + 2)/3,
    which gives a balanced multiplication in hgcd_matrix_adjust. However,
    p = 2 n/3 gives slightly better performance. The matrix-vector
    multiplication is then 4:1-unbalanced, with matrix elements of size n/6,
    and vector elements of size p = 2n/3.

    From analysis of the theoretical running time, it appears that when
    multiplication takes time O(n^alpha), p should be chosen so that the
    ratio of the time for the hgcd call, and the time for the multiplication
    in hgcd_matrix_adjust, is roughly 1/(alpha - 1).
*/

#ifndef FLINT_MPN_GCD_CHOOSE_P
# define FLINT_MPN_GCD_CHOOSE_P(n) (2 * (n) / 3)
#endif

#define EXTRACT_LIMB(count, xh, xl) \
    (((xh) << (count)) | (((xl) >> 1) >> ((FLINT_BITS - 1 - (count)))))

typedef struct
{
    mp_ptr gp;
    mp_size_t gn;
}
gcd_ctx_struct;

static void
gcd_hook(void * p, mp_srcptr gp, mp_size_t gn, mp_srcptr qp, mp_size_t qn, int d)
{
    gcd_ctx_struct * ctx = (gcd_ctx_struct *) p;
    (void) qp; (void) qn; (void) d;
    flint_mpn_copyi(ctx->gp, gp, gn);
    ctx->gn = gn;
}

/*
    Sets (gp, retval) to the gcd of (up, un) and (vp, vn). Requires
    un >= vn > 0, vp[vn - 1] != 0, and at least one of the operands odd.
    Both inputs are destroyed. gp needs room for vn limbs; it may be equal
    to up or vp but must not otherwise overlap them.
*/
mp_size_t
_flint_mpn_gcd_inplace(mp_ptr gp, mp_ptr up, mp_size_t usize, mp_ptr vp, mp_size_t n)
{
    mp_size_t talloc;
    mp_size_t scratch;
    mp_size_t matrix_scratch;

    gcd_ctx_struct ctx;
    mp_ptr tp;
    TMP_INIT;

    FLINT_ASSERT(usize >= n);
    FLINT_ASSERT(n > 0);
    FLINT_ASSERT(vp[n - 1] > 0);

    TMP_START;

    talloc = n;     /* subdiv step */

    /* For initial division */
    scratch = usize - n + 1;
    if (scratch > talloc)
        talloc = scratch;

    if (n >= FLINT_MPN_GCD_DC_CUTOFF)
    {
        mp_size_t hgcd_scratch;
        mp_size_t update_scratch;
        mp_size_t p = FLINT_MPN_GCD_CHOOSE_P(n);

        matrix_scratch = FLINT_MPN_HGCD_MATRIX_INIT_ITCH(n - p);
        hgcd_scratch = _flint_mpn_hgcd_itch(n - p);
        update_scratch = p + n - 1;

        scratch = matrix_scratch + FLINT_MAX(hgcd_scratch, update_scratch);
        if (scratch > talloc)
            talloc = scratch;
    }

    tp = TMP_ALLOC(talloc * sizeof(mp_limb_t));

    if (usize > n)
    {
        flint_mpn_tdiv_qr(tp, up, up, usize, vp, n);

        if (flint_mpn_zero_p(up, n))
        {
            flint_mpn_copyi(gp, vp, n);
            ctx.gn = n;
            goto done;
        }
    }

    ctx.gp = gp;

    while (n >= FLINT_MPN_GCD_DC_CUTOFF)
    {
        flint_mpn_hgcd_matrix_struct M;
        mp_size_t p = FLINT_MPN_GCD_CHOOSE_P(n);
        mp_size_t nn;

        matrix_scratch = FLINT_MPN_HGCD_MATRIX_INIT_ITCH(n - p);
        _flint_mpn_hgcd_matrix_init(&M, n - p, tp);
        nn = _flint_mpn_hgcd(up + p, vp + p, n - p, &M, tp + matrix_scratch);
        if (nn > 0)
        {
            FLINT_ASSERT(M.n <= (n - p - 1) / 2);
            FLINT_ASSERT(M.n + p <= (p + n - 1) / 2);
            /* Temporary storage 2 (p + M->n) <= p + n - 1. */
            n = _flint_mpn_hgcd_matrix_adjust(&M, p + nn, up, vp, p, tp + matrix_scratch);
        }
        else
        {
            /* Temporary storage n */
            n = _flint_mpn_gcd_subdiv_step(up, vp, n, 0, gcd_hook, &ctx, tp);
            if (n == 0)
                goto done;
        }
    }

    while (n > 2)
    {
        flint_mpn_hgcd_matrix1_struct M;
        mp_limb_t uh, ul, vh, vl;
        mp_limb_t mask;

        mask = up[n - 1] | vp[n - 1];
        FLINT_ASSERT(mask > 0);

        if (mask & (UWORD(1) << (FLINT_BITS - 1)))
        {
            uh = up[n - 1]; ul = up[n - 2];
            vh = vp[n - 1]; vl = vp[n - 2];
        }
        else
        {
            int shift = flint_clz(mask);
            uh = EXTRACT_LIMB(shift, up[n - 1], up[n - 2]);
            ul = EXTRACT_LIMB(shift, up[n - 2], up[n - 3]);
            vh = EXTRACT_LIMB(shift, vp[n - 1], vp[n - 2]);
            vl = EXTRACT_LIMB(shift, vp[n - 2], vp[n - 3]);
        }

        /* Try an hgcd2 step */
        if (_flint_mpn_hgcd2(uh, ul, vh, vl, &M))
        {
            n = _flint_mpn_hgcd_mul_matrix1_inverse_vector(&M, up, up, vp, n);
        }
        else
        {
            /* hgcd2 has failed. Then either one of a or b is very small, or
               the difference is very small. Perform one subtraction followed
               by one division. */

            /* Temporary storage n */
            n = _flint_mpn_gcd_subdiv_step(up, vp, n, 0, gcd_hook, &ctx, tp);
            if (n == 0)
                goto done;
        }
    }

    FLINT_ASSERT(up[n - 1] | vp[n - 1]);

    /* Due to the calling convention for _flint_mpn_gcd_inplace, at most one can be
       even. */
    if ((up[0] & 1) == 0)
        FLINT_SWAP(mp_ptr, up, vp);
    FLINT_ASSERT((up[0] & 1) != 0);

    {
        mp_limb_t u0, u1, v0, v1;

        u0 = up[0];
        v0 = vp[0];

        if (n == 1)
        {
            int cnt = flint_ctz(v0);
            *gp = mpn_gcd_11(u0, v0 >> cnt);
            ctx.gn = 1;
            goto done;
        }

        v1 = vp[1];
        if (FLINT_UNLIKELY(v0 == 0))
        {
            v0 = v1;
            v1 = 0;
        }
        if ((v0 & 1) == 0)
        {
            int cnt = flint_ctz(v0);
            v0 = (v1 << (FLINT_BITS - cnt)) | (v0 >> cnt);
            v1 >>= cnt;
        }

        u1 = up[1];
        _flint_mpn_gcd_22(gp, u1, u0, v1, v0);
        ctx.gn = 1 + (gp[1] > 0);
    }
done:
    TMP_END;
    return ctx.gn;
}
