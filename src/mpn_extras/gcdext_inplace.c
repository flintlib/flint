/*
    From GMP's mpn/generic/gcdext.c, gcdext_lehmer.c:
    Copyright 1996, 1998, 2000-2005, 2008, 2009, 2012 Free Software
    Foundation, Inc.

    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"
#include "impl.h"

#define MUL(rp, ap, an, bp, bn) \
    do { \
        if ((an) >= (bn)) \
            flint_mpn_mul(rp, ap, an, bp, bn); \
        else \
            flint_mpn_mul(rp, bp, bn, ap, an); \
    } while (0)

#define EXTRACT_LIMB(count, xh, xl) \
    (((xh) << (count)) | (((xl) >> 1) >> ((FLINT_BITS - 1 - (count)))))

/* Here, d is the index of the cofactor to update. */
void
_flint_mpn_gcdext_hook(void * p, mp_srcptr gp, mp_size_t gn,
    mp_srcptr qp, mp_size_t qn, int d)
{
    _flint_mpn_gcdext_ctx_struct * ctx = (_flint_mpn_gcdext_ctx_struct *) p;
    mp_size_t un = ctx->un;

    if (gp)
    {
        mp_srcptr up;

        FLINT_ASSERT(gn > 0);
        FLINT_ASSERT(gp[gn - 1] > 0);

        flint_mpn_copyi(ctx->gp, gp, gn);
        ctx->gn = gn;

        if (d < 0)
        {
            int c;

            /* Must return the smallest cofactor, +u1 or -u0 */
            c = mpn_cmp(ctx->u0, ctx->u1, un);
            FLINT_ASSERT(c != 0 || (un == 1 && ctx->u0[0] == 1 && ctx->u1[0] == 1));

            d = c < 0;
        }

        up = d ? ctx->u0 : ctx->u1;

        MPN_NORM(up, un);
        flint_mpn_copyi(ctx->up, up, un);

        *ctx->usize = d ? -un : un;
    }
    else
    {
        mp_limb_t cy;
        mp_ptr u0 = ctx->u0;
        mp_ptr u1 = ctx->u1;

        FLINT_ASSERT(d >= 0);

        if (d)
            FLINT_SWAP(mp_ptr, u0, u1);

        qn -= (qp[qn - 1] == 0);

        /* Update u0 += q  * u1 */
        if (qn == 1)
        {
            mp_limb_t q = qp[0];

            if (q == 1)
                /* A common case. */
                cy = mpn_add_n(u0, u0, u1, un);
            else
                cy = mpn_addmul_1(u0, u1, un, q);
        }
        else
        {
            mp_size_t u1n;
            mp_ptr tp;

            u1n = un;
            MPN_NORM(u1, u1n);

            if (u1n == 0)
                return;

            /* Should always have u1n == un here, and u1 >= u0. The reason
               is that we alternate adding u0 to u1 and u1 to u0
               (corresponding to subtractions a - b and b - a), and we can
               get a large quotient only just after a switch, which means
               that we'll add (a multiple of) the larger u to the smaller. */

            tp = ctx->tp;

            MUL(tp, qp, qn, u1, u1n);

            u1n += qn;
            u1n -= tp[u1n - 1] == 0;

            if (u1n >= un)
            {
                cy = mpn_add(u0, tp, u1n, u0, un);
                un = u1n;
            }
            else
                /* Note: Unlikely case, maybe never happens? */
                cy = mpn_add(u0, u0, un, tp, u1n);
        }
        u0[un] = cy;
        ctx->un = un + (cy > 0);
    }
}

/* Temporary storage: 3*(n+1) for u. If hgcd2 succeeds, we need n for the
   matrix-vector multiplication adjusting a, b. If hgcd fails, we need at
   most n for the quotient and n+1 for the u update (reusing the extra u).
   In all, 4n + 3. */
mp_size_t
_flint_mpn_gcdext_lehmer_n(mp_ptr gp, mp_ptr up, mp_size_t * usize,
    mp_ptr ap, mp_ptr bp, mp_size_t n, mp_ptr tp)
{
    mp_size_t ualloc = n + 1;

    /* Keeps track of the second row of the reduction matrix

       M = (v0, v1 ; u0, u1)

       which correspond to the first column of the inverse

       M^{-1} = (u1, -v1; -u0, v0)

       This implies that

       a =  u1 A (mod B)
       b = -u0 A (mod B)

       where A, B denotes the input values. */

    _flint_mpn_gcdext_ctx_struct ctx;
    mp_size_t un;
    mp_ptr u0;
    mp_ptr u1;
    mp_ptr u2;

    flint_mpn_zero(tp, 3 * ualloc);
    u0 = tp; tp += ualloc;
    u1 = tp; tp += ualloc;
    u2 = tp; tp += ualloc;

    u1[0] = 1; un = 1;

    ctx.gp = gp;
    ctx.up = up;
    ctx.usize = usize;

    while (n >= 2)
    {
        flint_mpn_hgcd_matrix1_struct M;
        mp_limb_t ah, al, bh, bl;
        mp_limb_t mask;

        mask = ap[n - 1] | bp[n - 1];
        FLINT_ASSERT(mask > 0);

        if (mask & (UWORD(1) << (FLINT_BITS - 1)))
        {
            ah = ap[n - 1]; al = ap[n - 2];
            bh = bp[n - 1]; bl = bp[n - 2];
        }
        else if (n == 2)
        {
            /* We use the full inputs without truncation, so we can safely
               shift left. */
            int shift = flint_clz(mask);
            ah = EXTRACT_LIMB(shift, ap[1], ap[0]);
            al = ap[0] << shift;
            bh = EXTRACT_LIMB(shift, bp[1], bp[0]);
            bl = bp[0] << shift;
        }
        else
        {
            int shift = flint_clz(mask);
            ah = EXTRACT_LIMB(shift, ap[n - 1], ap[n - 2]);
            al = EXTRACT_LIMB(shift, ap[n - 2], ap[n - 3]);
            bh = EXTRACT_LIMB(shift, bp[n - 1], bp[n - 2]);
            bl = EXTRACT_LIMB(shift, bp[n - 2], bp[n - 3]);
        }

        /* Try an hgcd2 step */
        if (_flint_mpn_hgcd2(ah, al, bh, bl, &M))
        {
            n = _flint_mpn_hgcd_mul_matrix1_inverse_vector(&M, ap, ap, bp, n);
            un = _flint_mpn_hgcd_mul_matrix1_vector(&M, u0, u0, u1, un);
        }
        else
        {
            /* hgcd2 has failed. Then either one of a or b is very small, or
               the difference is very small. Perform one subtraction followed
               by one division. */
            ctx.u0 = u0;
            ctx.u1 = u1;
            ctx.tp = u2;
            ctx.un = un;

            /* Temporary storage n for the quotient and ualloc for the new
               cofactor. */
            n = _flint_mpn_gcd_subdiv_step(ap, bp, n, 0, _flint_mpn_gcdext_hook, &ctx, tp);
            if (n == 0)
                return ctx.gn;

            un = ctx.un;
        }
    }

    FLINT_ASSERT(ap[0] > 0);
    FLINT_ASSERT(bp[0] > 0);

    if (ap[0] == bp[0])
    {
        int c;

        /* Which cofactor to return now? Candidates are +u1 and -u0,
           depending on which of a and b was most recently reduced, which we
           don't keep track of. So compare and get the smallest one. */

        gp[0] = ap[0];

        c = mpn_cmp(u0, u1, un);
        FLINT_ASSERT(c != 0 || (un == 1 && u0[0] == 1 && u1[0] == 1));
        if (c < 0)
        {
            MPN_NORM(u0, un);
            flint_mpn_copyi(up, u0, un);
            *usize = -un;
        }
        else
        {
            MPN_NORM(u1, un);
            flint_mpn_copyi(up, u1, un);
            *usize = un;
        }
        return 1;
    }
    else
    {
        mp_limb_t uh, vh;
        mp_limb_signed_t u;
        mp_limb_signed_t v;
        int negate;

        gp[0] = mpn_gcdext_1(&u, &v, ap[0], bp[0]);

        /* Set up = u u1 - v u0. Keep track of size, un grows by one or two
           limbs. */

        if (u == 0)
        {
            FLINT_ASSERT(v == 1);
            MPN_NORM(u0, un);
            flint_mpn_copyi(up, u0, un);
            *usize = -un;
            return 1;
        }
        else if (v == 0)
        {
            FLINT_ASSERT(u == 1);
            MPN_NORM(u1, un);
            flint_mpn_copyi(up, u1, un);
            *usize = un;
            return 1;
        }
        else if (u > 0)
        {
            negate = 0;
            FLINT_ASSERT(v < 0);
            v = -v;
        }
        else
        {
            negate = 1;
            FLINT_ASSERT(v > 0);
            u = -u;
        }

        uh = mpn_mul_1(up, u1, un, u);
        vh = mpn_addmul_1(up, u0, un, v);

        if ((uh | vh) > 0)
        {
            uh += vh;
            up[un++] = uh;
            if (uh < vh)
                up[un++] = 1;
        }

        MPN_NORM(up, un);

        *usize = negate ? -un : un;
        return 1;
    }
}

/* Computes (r;b) = (a; b) M. Result is of size n + M->n +/- 1, and the size
   is returned (if inputs are non-normalized, result may be non-normalized
   too). Temporary space needed is M->n + n. */
static mp_size_t
hgcd_mul_matrix_vector(flint_mpn_hgcd_matrix_struct * M,
    mp_ptr rp, mp_srcptr ap, mp_ptr bp, mp_size_t n, mp_ptr tp)
{
    mp_limb_t ah, bh;

#if FLINT_HAVE_FFT_SMALL
    {
        mp_size_t rn = _flint_mpn_hgcd_mul_matrix_vector_fft(M, rp, ap, bp, n);
        if (rn != 0)
            return rn;
    }
#endif

    /* Compute (r,b) <-- (u00 a + u10 b, u01 a + u11 b) as

       t  = u00 * a
       r  = u10 * b
       r += t;

       t  = u11 * b
       b  = u01 * a
       b += t;
    */

    MUL(tp, M->p[0][0], M->n, ap, n);
    MUL(rp, M->p[1][0], M->n, bp, n);

    ah = mpn_add_n(rp, rp, tp, n + M->n);

    MUL(tp, M->p[1][1], M->n, bp, n);
    MUL(bp, M->p[0][1], M->n, ap, n);

    bh = mpn_add_n(bp, bp, tp, n + M->n);

    n += M->n;
    if ((ah | bh) > 0)
    {
        rp[n] = ah;
        bp[n] = bh;
        n++;
    }
    else
    {
        /* Normalize */
        while ((rp[n - 1] | bp[n - 1]) == 0)
            n--;
    }

    return n;
}

/* Computes |v| = |(g - u a)| / b, where u may be positive or negative, and
   v is of the opposite sign. max(a, b) is of size n, u and v at most size n,
   and v must have space for n+1 limbs. */
static mp_size_t
compute_v(mp_ptr vp, mp_srcptr ap, mp_srcptr bp, mp_size_t n,
    mp_srcptr gp, mp_size_t gn, mp_srcptr up, mp_size_t usize, mp_ptr tp)
{
    mp_size_t size;
    mp_size_t an;
    mp_size_t bn;
    mp_size_t vn;

    FLINT_ASSERT(n > 0);
    FLINT_ASSERT(gn > 0);
    FLINT_ASSERT(usize != 0);

    size = FLINT_ABS(usize);
    FLINT_ASSERT(size <= n);
    FLINT_ASSERT(up[size - 1] > 0);

    an = n;
    MPN_NORM(ap, an);
    FLINT_ASSERT(gn <= an);

    MUL(tp, ap, an, up, size);

    size += an;

    if (usize > 0)
    {
        /* |v| = -v = (u a - g) / b */

        mpn_sub(tp, tp, size, gp, gn);
        MPN_NORM(tp, size);
        if (size == 0)
            return 0;
    }
    else
    {
        /* |v| = v = (g - u a) / b = (g + |u| a) / b. Since g <= a,
           (g + |u| a) always fits in (|usize| + an) limbs. */

        mpn_add(tp, tp, size, gp, gn);
        size -= (tp[size - 1] == 0);
    }

    /* Now divide t / b. There must be no remainder */
    bn = n;
    MPN_NORM(bp, bn);
    FLINT_ASSERT(size >= bn);

    vn = size + 1 - bn;
    FLINT_ASSERT(vn <= n + 1);

    flint_mpn_divexact(vp, tp, size, bp, bn);
    vn -= (vp[vn - 1] == 0);

    return vn;
}

/*
    Cofactor accumulation. GMP updates the cofactors after each half-gcd
    step as (u0, u1) <- (u0, u1) M_k, so every step costs products with the
    full current cofactor length, which grows to the input length, while
    the matrices M_k shrink geometrically (by about 2/3 per step). We
    instead collect the matrices and form Q = M_j ... M_k right to left,
    where each product costs about the size of the matrices involved, and
    apply Q to the cofactors once: (u0, u1) <- (u0, u1) Q. This is flushed
    before any step that needs the current cofactors (a subdivision step
    with the cofactor hook) and at the end of the main loop.
*/

typedef struct
{
    flint_mpn_hgcd_matrix_struct * M;
    mp_ptr * mem;
    slong len;
    slong alloc;
}
matrix_list_struct;

static void
matrix_list_init(matrix_list_struct * L)
{
    L->M = NULL;
    L->mem = NULL;
    L->len = 0;
    L->alloc = 0;
}

static void
matrix_list_clear(matrix_list_struct * L)
{
    slong i;
    for (i = 0; i < L->len; i++)
        flint_free(L->mem[i]);
    flint_free(L->M);
    flint_free(L->mem);
}

/* a matrix with entries of alloc limbs, zeroed, with M->n = 1 */
static mp_ptr
matrix_alloc(flint_mpn_hgcd_matrix_struct * M, mp_size_t alloc)
{
    mp_ptr p = flint_calloc(4 * alloc, sizeof(mp_limb_t));
    M->alloc = alloc;
    M->n = 1;
    M->p[0][0] = p;
    M->p[0][1] = p + alloc;
    M->p[1][0] = p + 2 * alloc;
    M->p[1][1] = p + 3 * alloc;
    return p;
}

static void
matrix_list_push(matrix_list_struct * L, const flint_mpn_hgcd_matrix_struct * M)
{
    mp_size_t n = M->n;
    int i, j;

    while (n > 1 && (M->p[0][0][n - 1] | M->p[0][1][n - 1] | M->p[1][0][n - 1] | M->p[1][1][n - 1]) == 0)
        n--;

    if (L->len == L->alloc)
    {
        L->alloc = FLINT_MAX(4, 2 * L->alloc);
        L->M = flint_realloc(L->M, L->alloc * sizeof(flint_mpn_hgcd_matrix_struct));
        L->mem = flint_realloc(L->mem, L->alloc * sizeof(mp_ptr));
    }

    L->mem[L->len] = matrix_alloc(L->M + L->len, n + 1);
    for (i = 0; i < 2; i++)
        for (j = 0; j < 2; j++)
            flint_mpn_copyi(L->M[L->len].p[i][j], M->p[i][j], n);
    L->M[L->len].n = n;
    L->len++;
}

static mp_size_t
hgcd_mul_matrix_vector(flint_mpn_hgcd_matrix_struct * M,
    mp_ptr rp, mp_srcptr ap, mp_ptr bp, mp_size_t n, mp_ptr tp);

/* (u0, u1) <- (u0, u1) M_1 ... M_k for the pending matrices, emptying
   the list; returns the new length of the cofactors, which have room for
   ualloc limbs. */
static mp_size_t
cofactors_flush(mp_ptr u0, mp_ptr u1, mp_size_t un, mp_size_t ualloc,
    matrix_list_struct * L)
{
    flint_mpn_hgcd_matrix_struct Q;
    mp_ptr Qmem, t0, t1, scratch;
    mp_size_t rn;
    slong i;

    if (L->len == 0)
        return un;

    /* Q = M_1 ... M_k, right to left */
    Q = L->M[L->len - 1];
    Qmem = L->mem[L->len - 1];

    for (i = L->len - 2; i >= 0; i--)
    {
        flint_mpn_hgcd_matrix_struct R;
        mp_ptr Rmem;
        const flint_mpn_hgcd_matrix_struct * A = L->M + i;
        int j, k;

        Rmem = matrix_alloc(&R, A->n + Q.n + 2);
        for (j = 0; j < 2; j++)
            for (k = 0; k < 2; k++)
                flint_mpn_copyi(R.p[j][k], A->p[j][k], A->n);
        R.n = A->n;

        scratch = flint_malloc(_flint_mpn_hgcd_matrix22_mul_strassen_itch(R.n, Q.n) * sizeof(mp_limb_t));
        _flint_mpn_hgcd_matrix_mul(&R, &Q, scratch);
        flint_free(scratch);

        flint_free(L->mem[i]);
        flint_free(Qmem);
        Q = R;
        Qmem = Rmem;
    }

    L->len = 0;

    /* (u0, u1) Q, through temporaries since the intermediate lengths may
       exceed that of the result */
    rn = un + Q.n + 1;
    t0 = flint_malloc((3 * rn + un + Q.n) * sizeof(mp_limb_t));
    t1 = t0 + rn;
    scratch = t1 + rn;

    flint_mpn_copyi(scratch, u0, un);
    flint_mpn_copyi(t1, u1, un);
    rn = hgcd_mul_matrix_vector(&Q, t0, scratch, t1, un, scratch + un);

    while ((t0[rn - 1] | t1[rn - 1]) == 0)
        rn--;
    FLINT_ASSERT(rn < ualloc);
    (void) ualloc;

    flint_mpn_copyi(u0, t0, rn);
    flint_mpn_copyi(u1, t1, rn);
    if (rn < un)
    {
        flint_mpn_zero(u0 + rn, un - rn);
        flint_mpn_zero(u1 + rn, un - rn);
    }

    flint_free(t0);
    flint_free(Qmem);

    return rn;
}

/*
    Temporary storage:

    Initial division: Quotient of at most an - n + 1 <= an limbs.

    Storage for u0 and u1: 2(n+1).

    Storage for hgcd matrix M, with input ceil(n/2): 5 * ceil(n/4)

    Storage for hgcd, input (n + 1)/2: 9 n/4 plus some.

    When hgcd succeeds: 1 + floor(3n/2) for adjusting a and b, and 2(n+1) for
    the cofactors.

    When hgcd fails: 2n + 1 for mpn_gcdext_subdiv_step, which is less.

    For the lehmer call after the loop, Let T denote
    FLINT_MPN_GCDEXT_DC_CUTOFF. For the gcdext_lehmer call, we need T each
    for u, a and b, and 4T+3 scratch space. Next, for compute_v, we need T
    for u, T+1 for v and 2T scratch space. In all, 7T + 3 is sufficient for
    both operations.
*/

/* Optimal choice of p seems difficult. In each iteration the division of
   work between hgcd and the updates of u0 and u1 depends on the current
   size of the u. It may be desirable to use a different choice of p in each
   iteration. Also the input size seems to matter; choosing p = n / 3 in the
   first iteration seems to improve performance slightly for input size just
   above the threshold, but degrade performance for larger inputs. */
#ifndef FLINT_MPN_GCDEXT_CHOOSE_P_1
# define FLINT_MPN_GCDEXT_CHOOSE_P_1(n) ((n) / 2)
#endif
#ifndef FLINT_MPN_GCDEXT_CHOOSE_P_2
# define FLINT_MPN_GCDEXT_CHOOSE_P_2(n) ((n) / 3)
#endif

/*
    Sets (gp, retval) to the gcd g of (ap, an) and (bp, n), and (up, *usizep)
    to a signed cofactor u with g = u a + v b for some v. Requires
    an >= n > 0 and bp[n - 1] != 0. Both inputs are destroyed, and ap needs
    room for an + 1 limbs. gp needs room for n limbs and up for n + 1 limbs.
*/
mp_size_t
_flint_mpn_gcdext_inplace(mp_ptr gp, mp_ptr up, mp_size_t * usizep,
    mp_ptr ap, mp_size_t an, mp_ptr bp, mp_size_t n)
{
    mp_size_t talloc;
    mp_size_t scratch;
    mp_size_t matrix_scratch;
    mp_size_t ualloc = n + 1;

    _flint_mpn_gcdext_ctx_struct ctx;
    mp_size_t un;
    mp_ptr u0;
    mp_ptr u1;
    matrix_list_struct pending;

    mp_ptr tp;

    TMP_INIT;

    FLINT_ASSERT(an >= n);
    FLINT_ASSERT(n > 0);
    FLINT_ASSERT(bp[n - 1] > 0);

    TMP_START;

    talloc = FLINT_MPN_GCDEXT_LEHMER_N_ITCH(n);

    /* For initial division */
    scratch = an - n + 1;
    if (scratch > talloc)
        talloc = scratch;

    matrix_scratch = 0;

    if (n >= FLINT_MPN_GCDEXT_DC_CUTOFF)
    {
        /* For hgcd loop. */
        mp_size_t hgcd_scratch;
        mp_size_t update_scratch;
        mp_size_t p1 = FLINT_MPN_GCDEXT_CHOOSE_P_1(n);
        mp_size_t p2 = FLINT_MPN_GCDEXT_CHOOSE_P_2(n);
        mp_size_t min_p = FLINT_MIN(p1, p2);
        mp_size_t max_p = FLINT_MAX(p1, p2);
        matrix_scratch = FLINT_MPN_HGCD_MATRIX_INIT_ITCH(n - min_p);
        hgcd_scratch = _flint_mpn_hgcd_itch(n - min_p);
        update_scratch = max_p + n - 1;

        scratch = matrix_scratch + FLINT_MAX(hgcd_scratch, update_scratch);
        if (scratch > talloc)
            talloc = scratch;

        /* Final gcdext_lehmer_n call. Need space for u and for copies of a
           and b. */
        scratch = FLINT_MPN_GCDEXT_LEHMER_N_ITCH(FLINT_MPN_GCDEXT_DC_CUTOFF)
            + 3 * FLINT_MPN_GCDEXT_DC_CUTOFF;

        if (scratch > talloc)
            talloc = scratch;

        /* Cofactors u0 and u1 */
        talloc += 2 * (n + 1);
    }

    tp = TMP_ALLOC(talloc * sizeof(mp_limb_t));

    if (an > n)
    {
        flint_mpn_tdiv_qr(tp, ap, ap, an, bp, n);

        if (flint_mpn_zero_p(ap, n))
        {
            flint_mpn_copyi(gp, bp, n);
            *usizep = 0;
            TMP_END;
            return n;
        }
    }

    if (n < FLINT_MPN_GCDEXT_DC_CUTOFF)
    {
        mp_size_t gn = _flint_mpn_gcdext_lehmer_n(gp, up, usizep, ap, bp, n, tp);

        TMP_END;
        return gn;
    }

    flint_mpn_zero(tp, 2 * ualloc);
    u0 = tp; tp += ualloc;
    u1 = tp; tp += ualloc;

    ctx.gp = gp;
    ctx.up = up;
    ctx.usize = usizep;

    {
        /* For the first hgcd call, there are no u updates, and it makes some
           sense to use a different choice for p. */

        flint_mpn_hgcd_matrix_struct M;
        mp_size_t p = FLINT_MPN_GCDEXT_CHOOSE_P_1(n);
        mp_size_t nn;

        _flint_mpn_hgcd_matrix_init(&M, n - p, tp);
        nn = _flint_mpn_hgcd(ap + p, bp + p, n - p, &M, tp + matrix_scratch);
        if (nn > 0)
        {
            FLINT_ASSERT(M.n <= (n - p - 1) / 2);
            FLINT_ASSERT(M.n + p <= (p + n - 1) / 2);

            /* Temporary storage 2 (p + M->n) <= p + n - 1 */
            n = _flint_mpn_hgcd_matrix_adjust(&M, p + nn, ap, bp, p, tp + matrix_scratch);

            flint_mpn_copyi(u0, M.p[1][0], M.n);
            flint_mpn_copyi(u1, M.p[1][1], M.n);
            un = M.n;
            while ((u0[un - 1] | u1[un - 1]) == 0)
                un--;
        }
        else
        {
            /* hgcd has failed. Then either one of a or b is very small, or
               the difference is very small. Perform one subtraction followed
               by one division. */
            u1[0] = 1;

            ctx.u0 = u0;
            ctx.u1 = u1;
            ctx.tp = tp + n; /* ualloc */
            ctx.un = 1;

            /* Temporary storage n */
            n = _flint_mpn_gcd_subdiv_step(ap, bp, n, 0, _flint_mpn_gcdext_hook, &ctx, tp);
            if (n == 0)
            {
                TMP_END;
                return ctx.gn;
            }

            un = ctx.un;
            FLINT_ASSERT(un < ualloc);
        }
    }

    matrix_list_init(&pending);

    while (n >= FLINT_MPN_GCDEXT_DC_CUTOFF)
    {
        flint_mpn_hgcd_matrix_struct M;
        mp_size_t p = FLINT_MPN_GCDEXT_CHOOSE_P_2(n);
        mp_size_t nn;

        _flint_mpn_hgcd_matrix_init(&M, n - p, tp);
        nn = _flint_mpn_hgcd(ap + p, bp + p, n - p, &M, tp + matrix_scratch);
        if (nn > 0)
        {
            mp_ptr t0;

            t0 = tp + matrix_scratch;
            FLINT_ASSERT(M.n <= (n - p - 1) / 2);
            FLINT_ASSERT(M.n + p <= (p + n - 1) / 2);

            /* Temporary storage 2 (p + M->n) <= p + n - 1 */
            n = _flint_mpn_hgcd_matrix_adjust(&M, p + nn, ap, bp, p, t0);

            /* The cofactor update (u0, u1) <- (u0, u1) M is deferred. */
            matrix_list_push(&pending, &M);
        }
        else
        {
            /* hgcd has failed. Then either one of a or b is very small, or
               the difference is very small. Perform one subtraction followed
               by one division. */
            un = cofactors_flush(u0, u1, un, ualloc, &pending);

            ctx.u0 = u0;
            ctx.u1 = u1;
            ctx.tp = tp + n; /* ualloc */
            ctx.un = un;

            /* Temporary storage n */
            n = _flint_mpn_gcd_subdiv_step(ap, bp, n, 0, _flint_mpn_gcdext_hook, &ctx, tp);
            if (n == 0)
            {
                matrix_list_clear(&pending);
                TMP_END;
                return ctx.gn;
            }

            un = ctx.un;
            FLINT_ASSERT(un < ualloc);
        }
    }

    un = cofactors_flush(u0, u1, un, ualloc, &pending);
    matrix_list_clear(&pending);

    FLINT_ASSERT(un < ualloc);
    FLINT_ASSERT((u0[un - 1] | u1[un - 1]) > 0);

    /* We have A = ... a + ... b
               B =  u0 a +  u1 b

               a = u1  A + ... B
               b = -u0 A + ... B

       with bounds

         |u0|, |u1| <= B / min(a, b)

       We always have u1 > 0, and u0 == 0 is possible only if u1 == 1, in
       which case the only reduction done so far is a = A - k B for some k.

       Compute g = u a + v b = (u u1 - v u0) A + (...) B
       Here, u, v are bounded by

         |u| <= b,
         |v| <= a
    */

    FLINT_ASSERT((ap[n - 1] | bp[n - 1]) > 0);

    if (FLINT_UNLIKELY(mpn_cmp(ap, bp, n) == 0))
    {
        /* Must return the smallest cofactor, +u1 or -u0 */
        int c;

        flint_mpn_copyi(gp, ap, n);

        c = mpn_cmp(u0, u1, un);
        /* c == 0 can happen only when A = (2k+1) G, B = 2 G. And in this
           case we choose the cofactor + 1, corresponding to G = A - k B,
           rather than -1, corresponding to G = - A + (k+1) B. */
        FLINT_ASSERT(c != 0 || (un == 1 && u0[0] == 1 && u1[0] == 1));
        if (c < 0)
        {
            MPN_NORM(u0, un);
            flint_mpn_copyi(up, u0, un);
            *usizep = -un;
        }
        else
        {
            MPN_NORM(u1, un);
            flint_mpn_copyi(up, u1, un);
            *usizep = un;
        }

        TMP_END;
        return n;
    }
    else if (FLINT_UNLIKELY(u0[0] == 0) && un == 1)
    {
        mp_size_t gn;
        FLINT_ASSERT(u1[0] == 1);

        /* g = u a + v b = (u u1 - v u0) A + (...) B = u A + (...) B */
        gn = _flint_mpn_gcdext_lehmer_n(gp, up, usizep, ap, bp, n, tp);

        TMP_END;
        return gn;
    }
    else
    {
        mp_size_t u0n;
        mp_size_t u1n;
        mp_size_t lehmer_un;
        mp_size_t lehmer_vn;
        mp_size_t gn;

        mp_ptr lehmer_up;
        mp_ptr lehmer_vp;
        int negate;

        lehmer_up = tp; tp += n;

        /* Call gcdext_lehmer_n with copies of a and b. */
        flint_mpn_copyi(tp, ap, n);
        flint_mpn_copyi(tp + n, bp, n);
        gn = _flint_mpn_gcdext_lehmer_n(gp, lehmer_up, &lehmer_un, tp, tp + n, n, tp + 2 * n);

        u0n = un;
        MPN_NORM(u0, u0n);
        FLINT_ASSERT(u0n > 0);

        if (lehmer_un == 0)
        {
            /* u == 0  ==>  v = g / b == 1  ==> g = - u0 A + (...) B */
            flint_mpn_copyi(up, u0, u0n);
            *usizep = -u0n;

            TMP_END;
            return gn;
        }

        lehmer_vp = tp;
        /* Compute v = (g - u a) / b */
        lehmer_vn = compute_v(lehmer_vp, ap, bp, n, gp, gn, lehmer_up, lehmer_un, tp + n + 1);

        if (lehmer_un > 0)
            negate = 0;
        else
        {
            lehmer_un = -lehmer_un;
            negate = 1;
        }

        u1n = un;
        MPN_NORM(u1, u1n);
        FLINT_ASSERT(u1n > 0);

        FLINT_ASSERT(lehmer_un + u1n <= ualloc);
        FLINT_ASSERT(lehmer_vn + u0n <= ualloc);

        /* We may still have v == 0 */

        /* Compute u u0 */
        MUL(up, u1, u1n, lehmer_up, lehmer_un);

        un = u1n + lehmer_un;
        un -= (up[un - 1] == 0);

        if (lehmer_vn > 0)
        {
            mp_limb_t cy;

            /* Overwrites old u1 value */
            MUL(u1, u0, u0n, lehmer_vp, lehmer_vn);

            u1n = u0n + lehmer_vn;
            u1n -= (u1[u1n - 1] == 0);

            if (u1n <= un)
            {
                cy = mpn_add(up, up, un, u1, u1n);
            }
            else
            {
                cy = mpn_add(up, u1, u1n, up, un);
                un = u1n;
            }
            up[un] = cy;
            un += (cy != 0);

            FLINT_ASSERT(un < ualloc);
        }
        *usizep = negate ? -un : un;

        TMP_END;
        return gn;
    }
}
