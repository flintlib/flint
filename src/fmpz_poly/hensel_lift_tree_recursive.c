/*
    Copyright (C) 2011 Andy Novocin
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "gr.h"
#include "gr_vec.h"
#include "gr_poly.h"
#include "mpn_mod.h"

/*
    One node of the Hensel tree.

    On entry f = g h and a g + b h = 1 modulo p0, with g and h monic; on
    exit the same holds modulo P = p0 p1.  Writing C = ((f - g h) mod P)/p0
    and C2 = ((1 - a G - b H) mod P)/p0, the updates are

        G  = g + p0 ((C  b) mod g mod p1)
        H  = f / G                              (exact in Z/P)
        B  = b + p0 ((C2 b) mod g mod p1)
        A  = a + p0 ((C2 a) mod h mod p1).

    This differs from calling fmpz_poly_hensel_lift() once per node in
    three ways.

    (1) The two halves are fused.  Since p1 | p0 we have G = g and H = h
        modulo p1, so the reduced moduli g mod p1 and h mod p1 and their
        Newton preinverses are shared by all of the divisions above.  The
        split version recomputes each of them four times.

    (2) The arithmetic runs over a gr ring chosen by the size of the
        modulus: nmod for a single limb, mpn_mod for a few limbs and
        fmpz_mod above that.  Quadratic lifting walks a chain of moduli
        p, p^2, p^4, ..., so most levels are small even when the final
        precision is large.

    (3) H is obtained by an exact division instead of a second lift.  This
        requires the constant coefficient of G to be a unit, which holds
        whenever p does not divide f(0); the caller tests this once for
        the whole tree and falls back to the lift otherwise.
*/

/* Beyond this many limbs mpn_mod loses to fmpz_mod at the polynomial
   lengths that arise here.  This is deliberately below
   MPN_MOD_MAX_LIMBS. */
#define HENSEL_MPN_MOD_MAX_LIMBS 12

/* Nodes with lenG*lenH below this run the original code; see the
   comment at the dispatch in _hensel_lift_node(). */
#define HENSEL_MIN_AREA 36

/* Below this many bits in the modulus, a shared Newton preinverse only
   pays off for long polynomials. */
#define HENSEL_PREINV_BIG_BITS 768
#define HENSEL_PREINV_MIN_LEN 16
#define HENSEL_PREINV_MIN_LEN_SMALL 64
#define HENSEL_PREINV_MIN_QUOT 8

typedef struct
{
    gr_ctx_struct R1[1];    /* Z/p1: reductions modulo g and h             */
    gr_ctx_struct RP[1];    /* Z/(p0 p1): products and the exact division  */
    fmpz_t p0;
    fmpz_t p1;
    fmpz_t P;
    slong sz1;
    slong szP;
    flint_bitcnt_t p1bits;
    int use_divexact;   /* -1 until resolved */
    int R1_ready;
    int RP_ready;
    const fmpz_poly_struct * f;
}
hensel_ctx_struct;

typedef hensel_ctx_struct hensel_ctx_t[1];

static void
_hensel_ring_init(gr_ctx_t R, const fmpz_t m)
{
    flint_bitcnt_t bits = fmpz_bits(m);

    if (bits <= FLINT_BITS)
        gr_ctx_init_nmod(R, fmpz_get_ui(m));
    else if (bits > HENSEL_MPN_MOD_MAX_LIMBS * FLINT_BITS ||
             gr_ctx_init_mpn_mod(R, m) != GR_SUCCESS)
        gr_ctx_init_fmpz_mod(R, m);
}

static void
_hensel_ctx_init(hensel_ctx_t ctx, const fmpz_t p0, const fmpz_t p1,
                 const fmpz_poly_t f)
{
    fmpz_init_set(ctx->p0, p0);
    fmpz_init_set(ctx->p1, p1);
    fmpz_init(ctx->P);
    fmpz_mul(ctx->P, p0, p1);

    /*
        Both rings are initialised lazily.  Building a context for a
        modulus of tens of thousands of bits precomputes an inverse,
        which is wasted work on a tree whose nodes are all small enough
        to take the original code path below.
    */
    ctx->R1_ready = 0;
    ctx->RP_ready = 0;
    ctx->sz1 = 0;
    ctx->szP = 0;
    ctx->f = f;
    ctx->p1bits = fmpz_bits(p1);

    /* Resolved on first use; see _hensel_ctx_use_divexact(). */
    ctx->use_divexact = -1;
}

static void
_hensel_ctx_ensure_RP(hensel_ctx_struct * ctx)
{
    if (!ctx->RP_ready)
    {
        _hensel_ring_init(ctx->RP, ctx->P);
        ctx->szP = ctx->RP->sizeof_elem;
        ctx->RP_ready = 1;
    }
}

static void
_hensel_ctx_ensure_R1(hensel_ctx_struct * ctx)
{
    if (!ctx->R1_ready)
    {
        _hensel_ring_init(ctx->R1, ctx->p1);
        ctx->sz1 = ctx->R1->sizeof_elem;
        ctx->R1_ready = 1;
    }
}

/*
    Whether the sibling factor may be recovered by an exact division,
    which needs the constant coefficient of every G in the tree to be a
    unit.  Each factor divides f modulo p, so testing the root suffices.
*/
static int
_hensel_ctx_use_divexact(hensel_ctx_struct * ctx)
{
    if (ctx->use_divexact < 0)
    {
        const fmpz_poly_struct * f = ctx->f;

        if (f == NULL || f->length == 0)
            ctx->use_divexact = 0;
        else
        {
            fmpz_t t;
            fmpz_init(t);
            fmpz_gcd(t, f->coeffs, ctx->p1);
            ctx->use_divexact = fmpz_is_one(t);
            fmpz_clear(t);
        }
    }

    return ctx->use_divexact;
}

static void
_hensel_ctx_clear(hensel_ctx_t ctx)
{
    if (ctx->R1_ready)
        gr_ctx_clear(ctx->R1);
    if (ctx->RP_ready)
        gr_ctx_clear(ctx->RP);
    fmpz_clear(ctx->p0);
    fmpz_clear(ctx->p1);
    fmpz_clear(ctx->P);
}

static gr_ptr
_hensel_vec_init(slong len, gr_ctx_t R)
{
    gr_ptr v = flint_malloc(len * R->sizeof_elem);
    _gr_vec_init(v, len, R);
    return v;
}

static void
_hensel_vec_clear(gr_ptr v, slong len, gr_ctx_t R)
{
    _gr_vec_clear(v, len, R);
    flint_free(v);
}

static void
_hensel_set_fmpz_vec(gr_ptr res, const fmpz * v, slong len, gr_ctx_t R)
{
    slong i, sz = R->sizeof_elem;

    for (i = 0; i < len; i++)
        GR_MUST_SUCCEED(gr_set_fmpz(GR_ENTRY(res, i, sz), v + i, R));
}

static void
_hensel_get_fmpz_vec(fmpz * res, gr_srcptr v, slong len, gr_ctx_t R)
{
    slong i, sz = R->sizeof_elem;

    for (i = 0; i < len; i++)
        GR_MUST_SUCCEED(gr_get_fmpz(res + i, GR_ENTRY(v, i, sz), R));
}

/*
    Whether to precompute a Newton inverse of a reduced modulus of length
    lenM, for divisions with quotients of length up to lenQ.  The inverse
    is shared by every division against that modulus in the node, so it
    pays off well before the generic divrem cutoff, and earlier still as
    the coefficients grow.
*/
static int
_hensel_want_preinv(slong lenM, slong lenQ, flint_bitcnt_t bits)
{
    if (lenQ < HENSEL_PREINV_MIN_QUOT || lenM < HENSEL_PREINV_MIN_LEN)
        return 0;

    return (bits >= HENSEL_PREINV_BIG_BITS) ? 1
                                            : (lenM >= HENSEL_PREINV_MIN_LEN_SMALL);
}

/* R = A mod M, of length lenM - 1.  Q is scratch of length >= lenA - lenM + 1. */
static void
_hensel_rem(gr_ptr R, gr_srcptr A, slong lenA, gr_srcptr M, slong lenM,
            gr_srcptr Minv, slong lenMinv, gr_ptr Q, gr_ctx_t R1)
{
    slong sz = R1->sizeof_elem;

    if (lenA < lenM)
    {
        GR_MUST_SUCCEED(_gr_vec_set(R, A, lenA, R1));
        GR_MUST_SUCCEED(_gr_vec_zero(GR_ENTRY(R, lenA, sz), lenM - 1 - lenA, R1));
    }
    else if (Minv != NULL && lenMinv >= lenA - lenM + 1)
        GR_MUST_SUCCEED(_gr_poly_divrem_newton_n_preinv(Q, R, A, lenA, M, lenM,
                                                        Minv, lenMinv, R1));
    else
        GR_MUST_SUCCEED(_gr_poly_divrem(Q, R, A, lenA, M, lenM, R1));
}

/*
    res[0, lenM-1) = base[0, lenBase) + p0 ((C cof) mod M), over Z, where
    lenBase <= lenM - 1 and res may alias base.  D, E, Q are scratch over
    R1 and t is scratch over Z.
*/
static void
_hensel_correct(fmpz * res, const fmpz * base, slong lenBase,
                gr_srcptr C, slong lenC, gr_srcptr cof, slong lenCof,
                gr_srcptr M, slong lenM, gr_srcptr Minv, slong lenMinv,
                gr_ptr D, gr_ptr E, gr_ptr Q, fmpz * t,
                const hensel_ctx_struct * ctx)
{
    gr_ctx_struct * R1 = (gr_ctx_struct *) ctx->R1;

    _hensel_rem(D, C, lenC, M, lenM, Minv, lenMinv, Q, R1);

    if (lenCof > 1)
    {
        if (lenM - 1 >= lenCof)
            GR_MUST_SUCCEED(_gr_poly_mul(E, D, lenM - 1, cof, lenCof, R1));
        else
            GR_MUST_SUCCEED(_gr_poly_mul(E, cof, lenCof, D, lenM - 1, R1));

        _hensel_rem(D, E, lenM + lenCof - 2, M, lenM, Minv, lenMinv, Q, R1);
    }
    else
        GR_MUST_SUCCEED(_gr_vec_mul_scalar(D, D, lenM - 1, cof, R1));

    _hensel_get_fmpz_vec(t, D, lenM - 1, R1);
    _fmpz_vec_scalar_mul_fmpz(t, t, lenM - 1, ctx->p0);
    _fmpz_poly_add(res, t, lenM - 1, base, lenBase);
}

static void
_hensel_lift_node(fmpz_poly_t G, fmpz_poly_t H, fmpz_poly_t A, fmpz_poly_t B,
                  const fmpz_poly_t f, int inv, hensel_ctx_struct * ctx)
{
    gr_ctx_struct * R1 = ctx->R1;
    gr_ctx_struct * RP = ctx->RP;

    slong sz1;
    slong szP;
    const slong lenG = G->length;
    const slong lenH = H->length;
    const slong lenA = A->length;
    const slong lenB = B->length;
    const slong lenF = lenG + lenH - 1;
    const slong lenC2 = FLINT_MAX(lenA + lenG - 1, lenB + lenH - 1);
    const slong Lg = FLINT_MAX(lenH, lenG - 1);
    const slong Lh = FLINT_MAX(lenG, lenH - 1);
    const slong lenmax = FLINT_MAX(lenF, lenC2) + FLINT_MAX(lenG, lenH) + 2;

    gr_ptr Mg, Mh, ca, cb, Mginv, Mhinv, C, D, E, Q;
    gr_ptr xP, yP, zP, wP;
    fmpz * t;
    fmpz * zv;
    slong i, n1, nP;
    int use_ringP, use_dx, need_h, newtG, newtH;

    /* The products f - g h and 1 - a G - b H are needed only modulo P.
       Computing them in Z/P pays for itself once the polynomials are
       long enough to amortise the coefficient reductions, which is the
       same regime in which the exact division below is worthwhile; for
       short polynomials with huge coefficients we keep the plain integer
       arithmetic instead. */
    /*
        For a node this small the work is dominated by per-coefficient
        overhead rather than by the polynomial arithmetic, and none of
        the optimisations here apply: the preinverse is not worth
        building, and the exact division costs one modular inverse
        (about six modular multiplications at any size) to save only
        O(lenG lenH) coefficient operations.  Run the original node code
        unchanged, which also avoids converting coefficients into and
        out of a gr ring.
    */
    if (lenG * lenH < HENSEL_MIN_AREA)
    {
        if (inv == 1)
            fmpz_poly_hensel_lift(G, H, A, B, f, G, H, A, B, ctx->p0, ctx->p1);
        else if (inv == -1)
            fmpz_poly_hensel_lift_only_inverse(A, B, G, H, A, B,
                                               ctx->p0, ctx->p1);
        else
            fmpz_poly_hensel_lift_without_inverse(G, H, f, G, H, A, B,
                                                  ctx->p0, ctx->p1);
        return;
    }

    _hensel_ctx_ensure_R1(ctx);
    sz1 = ctx->sz1;

    use_ringP = (inv >= 0) && _hensel_ctx_use_divexact(ctx);
    use_dx = use_ringP;

    if (use_ringP)
        _hensel_ctx_ensure_RP(ctx);
    szP = ctx->szP;

    n1 = lenG + lenH + lenA + lenB + Lg + Lh + 4 * lenmax + 8;
    Mg = _hensel_vec_init(n1, R1);
    Mh    = GR_ENTRY(Mg, lenG, sz1);
    ca    = GR_ENTRY(Mh, lenH, sz1);
    cb    = GR_ENTRY(ca, lenA, sz1);
    Mginv = GR_ENTRY(cb, lenB, sz1);
    Mhinv = GR_ENTRY(Mginv, Lg + 4, sz1);
    C     = GR_ENTRY(Mhinv, Lh + 4, sz1);
    D     = GR_ENTRY(C, lenmax, sz1);
    E     = GR_ENTRY(D, lenmax, sz1);
    Q     = GR_ENTRY(E, lenmax, sz1);

    if (use_ringP)
    {
        nP = 4 * lenmax;
        xP = _hensel_vec_init(nP, RP);
        yP = GR_ENTRY(xP, lenmax, szP);
        zP = GR_ENTRY(xP, 2 * lenmax, szP);
        wP = GR_ENTRY(xP, 3 * lenmax, szP);
        zv = _fmpz_vec_init(lenmax);
    }
    else
    {
        nP = 0;
        xP = yP = zP = wP = NULL;
        zv = _fmpz_vec_init(3 * lenmax);
    }

    t = zv + ((use_ringP) ? 0 : 2 * lenmax);

    FLINT_ASSERT(lenG >= 2 && lenH >= 2);
    FLINT_ASSERT(lenA >= 1 && lenB >= 1);
    FLINT_ASSERT(inv < 0 || f->length == lenF);

    fmpz_poly_fit_length(G, lenG);
    fmpz_poly_fit_length(H, lenH);
    if (inv != 0)
    {
        fmpz_poly_fit_length(A, FLINT_MAX(lenH - 1, lenA));
        fmpz_poly_fit_length(B, FLINT_MAX(lenG - 1, lenB));
    }

    /* Reduced moduli and cofactors.  These are also correct for the
       second half, since p1 | p0 gives G = g and H = h modulo p1. */
    _hensel_set_fmpz_vec(Mg, G->coeffs, lenG, R1);
    _hensel_set_fmpz_vec(Mh, H->coeffs, lenH, R1);
    _hensel_set_fmpz_vec(ca, A->coeffs, lenA, R1);
    _hensel_set_fmpz_vec(cb, B->coeffs, lenB, R1);

    need_h = (inv != 0) || !use_dx;

    newtG = _hensel_want_preinv(lenG, Lg, ctx->p1bits);
    newtH = need_h && _hensel_want_preinv(lenH, Lh, ctx->p1bits);

    if (newtG)
    {
        GR_MUST_SUCCEED(_gr_poly_reverse(D, Mg, lenG, lenG, R1));
        GR_MUST_SUCCEED(_gr_poly_inv_series(Mginv, D, lenG, Lg, R1));
    }

    if (newtH)
    {
        GR_MUST_SUCCEED(_gr_poly_reverse(D, Mh, lenH, lenH, R1));
        GR_MUST_SUCCEED(_gr_poly_inv_series(Mhinv, D, lenH, Lh, R1));
    }

    if (inv >= 0)
    {
        if (use_ringP)
        {
            _hensel_set_fmpz_vec(xP, G->coeffs, lenG, RP);
            _hensel_set_fmpz_vec(yP, H->coeffs, lenH, RP);

            if (lenG >= lenH)
                GR_MUST_SUCCEED(_gr_poly_mul(zP, xP, lenG, yP, lenH, RP));
            else
                GR_MUST_SUCCEED(_gr_poly_mul(zP, yP, lenH, xP, lenG, RP));

            _hensel_set_fmpz_vec(xP, f->coeffs, lenF, RP);
            GR_MUST_SUCCEED(_gr_vec_sub(zP, xP, zP, lenF, RP));

            _hensel_get_fmpz_vec(t, zP, lenF, RP);
            for (i = 0; i < lenF; i++)
                fmpz_divexact(t + i, t + i, ctx->p0);
        }
        else
        {
            if (lenG >= lenH)
                _fmpz_poly_mul(zv, G->coeffs, lenG, H->coeffs, lenH);
            else
                _fmpz_poly_mul(zv, H->coeffs, lenH, G->coeffs, lenG);

            _fmpz_vec_sub(zv, f->coeffs, zv, lenF);
            _fmpz_vec_scalar_divexact_fmpz(t, zv, lenF, ctx->p0);
        }

        /* C = ((f - g h) mod P)/p0, reduced modulo p1. */
        _hensel_set_fmpz_vec(C, t, lenF, R1);

        _hensel_correct(G->coeffs, G->coeffs, lenG - 1, C, lenF, cb, lenB,
                        Mg, lenG, newtG ? Mginv : NULL, Lg, D, E, Q, t, ctx);
        fmpz_one(G->coeffs + lenG - 1);
        _fmpz_poly_set_length(G, lenG);

        if (use_dx)
        {
            /* xP still holds f modulo P. */
            _hensel_set_fmpz_vec(yP, G->coeffs, lenG, RP);

            if (_gr_poly_divexact(zP, xP, lenF, yP, lenG, RP) == GR_SUCCESS)
                _hensel_get_fmpz_vec(H->coeffs, zP, lenH, RP);
            else
                use_dx = 0;
        }

        if (!use_dx)
        {
            _hensel_correct(H->coeffs, H->coeffs, lenH - 1, C, lenF, ca, lenA,
                            Mh, lenH, newtH ? Mhinv : NULL, Lh, D, E, Q, t, ctx);
            fmpz_one(H->coeffs + lenH - 1);
        }

        _fmpz_poly_set_length(H, lenH);
    }

    if (inv != 0)
    {
        if (use_ringP)
        {
            GR_MUST_SUCCEED(_gr_vec_zero(xP, nP, RP));

            _hensel_set_fmpz_vec(xP, G->coeffs, lenG, RP);
            _hensel_set_fmpz_vec(yP, A->coeffs, lenA, RP);

            if (lenG >= lenA)
                GR_MUST_SUCCEED(_gr_poly_mul(zP, xP, lenG, yP, lenA, RP));
            else
                GR_MUST_SUCCEED(_gr_poly_mul(zP, yP, lenA, xP, lenG, RP));

            GR_MUST_SUCCEED(_gr_vec_zero(xP, lenmax, RP));
            GR_MUST_SUCCEED(_gr_vec_zero(yP, lenmax, RP));
            _hensel_set_fmpz_vec(xP, H->coeffs, lenH, RP);
            _hensel_set_fmpz_vec(yP, B->coeffs, lenB, RP);

            if (lenH >= lenB)
                GR_MUST_SUCCEED(_gr_poly_mul(wP, xP, lenH, yP, lenB, RP));
            else
                GR_MUST_SUCCEED(_gr_poly_mul(wP, yP, lenB, xP, lenH, RP));

            GR_MUST_SUCCEED(_gr_vec_add(zP, zP, wP, lenC2, RP));
            GR_MUST_SUCCEED(_gr_vec_neg(zP, zP, lenC2, RP));
            GR_MUST_SUCCEED(gr_add_ui(zP, zP, 1, RP));

            _hensel_get_fmpz_vec(t, zP, lenC2, RP);
            for (i = 0; i < lenC2; i++)
                fmpz_divexact(t + i, t + i, ctx->p0);
        }
        else
        {
            fmpz * wv = zv + lenmax;

            _fmpz_vec_zero(zv, 2 * lenmax);

            if (lenG >= lenA)
                _fmpz_poly_mul(zv, G->coeffs, lenG, A->coeffs, lenA);
            else
                _fmpz_poly_mul(zv, A->coeffs, lenA, G->coeffs, lenG);

            if (lenH >= lenB)
                _fmpz_poly_mul(wv, H->coeffs, lenH, B->coeffs, lenB);
            else
                _fmpz_poly_mul(wv, B->coeffs, lenB, H->coeffs, lenH);

            _fmpz_vec_add(zv, zv, wv, lenC2);
            fmpz_sub_ui(zv, zv, 1);
            _fmpz_vec_neg(zv, zv, lenC2);
            _fmpz_vec_scalar_divexact_fmpz(t, zv, lenC2, ctx->p0);
        }

        _hensel_set_fmpz_vec(C, t, lenC2, R1);

        _hensel_correct(B->coeffs, B->coeffs, lenB, C, lenC2, cb, lenB,
                        Mg, lenG, newtG ? Mginv : NULL, Lg, D, E, Q, t, ctx);
        _fmpz_poly_set_length(B, lenG - 1);
        _fmpz_poly_normalise(B);

        _hensel_correct(A->coeffs, A->coeffs, lenA, C, lenC2, ca, lenA,
                        Mh, lenH, newtH ? Mhinv : NULL, Lh, D, E, Q, t, ctx);
        _fmpz_poly_set_length(A, lenH - 1);
        _fmpz_poly_normalise(A);
    }

    _fmpz_vec_clear(zv, use_ringP ? lenmax : 3 * lenmax);
    if (use_ringP)
        _hensel_vec_clear(xP, nP, RP);
    _hensel_vec_clear(Mg, n1, R1);
}

static void
_hensel_lift_tree_recursive(slong * link, fmpz_poly_t * v, fmpz_poly_t * w,
                            fmpz_poly_t f, slong j, int inv,
                            hensel_ctx_struct * ctx)
{
    if (j >= 0)
    {
        _hensel_lift_node(v[j], v[j + 1], w[j], w[j + 1], f, inv, ctx);
        _hensel_lift_tree_recursive(link, v, w, v[j], link[j], inv, ctx);
        _hensel_lift_tree_recursive(link, v, w, v[j + 1], link[j + 1], inv, ctx);
    }
}

void fmpz_poly_hensel_lift_tree_recursive(slong *link,
    fmpz_poly_t *v, fmpz_poly_t *w, fmpz_poly_t f, slong j, slong inv,
    const fmpz_t p0, const fmpz_t p1)
{
    hensel_ctx_t ctx;

    if (j < 0)
        return;

    _hensel_ctx_init(ctx, p0, p1, f);
    _hensel_lift_tree_recursive(link, v, w, f, j, (int) inv, ctx);
    _hensel_ctx_clear(ctx);
}
