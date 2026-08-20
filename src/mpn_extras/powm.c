/*
    Copyright (C) 2026 Fredrik Johansson
    Developed using Claude Fable 5

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"
#include "ulong_extras.h"
#include "nmod.h"

#define POWM_PROFILE 0
#if FLINT_HAVE_FFT_SMALL
#include "fft_small.h"
#endif

/* below this, mulders mulhigh reduction; above, the wraparound fold
   (provisional, tuned below) */
#define FLINT_MPN_POWM_REDC_FOLD_THRESHOLD FLINT_MPN_POWM_BNM1_FOLD_THRESHOLD

/* above this, the quotient q = X_lo * minv mod B^mn also runs in the
   transform domain against a cached transform of minv */
#define FLINT_MPN_POWM_REDC_QSTEP_FFT_THRESHOLD 480

/* moduli below the fold tiers with large exponents go to mpz_powm
   for now: there the reduction is plain mulhigh + mullow against
   GMP's redc_1/redc_2 assembly and low-threshold mulmod_bnm1 */
#define FLINT_MPN_POWM_MPZ_FALLBACK_LIMBS 110

/* in the basecase, replace mulmod_preinvn's mullow by the wraparound
   fold in this limb range: below, mulders mullow wins; above,
   flint_mpn_mullow_n dispatches to fft_small and beats the Toom-based
   bnm1 recursion */
#define FLINT_MPN_POWM_BASECASE_FOLD_MIN 120
#define FLINT_MPN_POWM_BASECASE_FOLD_MAX 2000

/* below this many exponent bits, a caller-supplied preinvn inverse
   makes the division-based basecase cheaper than paying the redc
   setup (binvert, Montgomery conversions, fold initialization) */
#define FLINT_MPN_POWM_PREINVN_REDC_EBITS 24

/* fold tiers, from in-situ per-reduction timings: mulders mulhigh,
   then flint_mpn_mulmod_bnm1 (small working set, no fixed transform
   costs), then the cached chain or cyclic-plan fold */
#define FLINT_MPN_POWM_BNM1_FOLD_THRESHOLD 120
#define FLINT_MPN_POWM_LARGE_FOLD_THRESHOLD 480
#define FLINT_MPN_POWM_CHAIN_THRESHOLD 120
#define FLINT_MPN_POWM_CHAIN_NEG_H 128
#define FLINT_MPN_POWM_CHAIN_MAX_LEVELS 40
#include "longlong.h"

#define TSTBIT(x, b) ((x[(b) / FLINT_BITS] >> ((b) % FLINT_BITS)) & UWORD(1))

/* window size for an ebits-bit exponent, as in _gr_pow_mpn_sliding,
   capped so that the table of 2^(k-1) n-limb entries stays modest */
static int
powm_select_k(flint_bitcnt_t ebits, mp_size_t n)
{
    int k;

    if      (ebits <= 8)     k = 1;
    else if (ebits <= 24)    k = 2;
    else if (ebits <= 69)    k = 3;
    else if (ebits <= 196)   k = 4;
    else if (ebits <= 538)   k = 5;
    else if (ebits <= 1433)  k = 6;
    else if (ebits <= 3714)  k = 7;
    else if (ebits <= 9399)  k = 8;
    else if (ebits <= 23290) k = 9;
    else if (ebits <= 56651) k = 10;
    else                     k = 11;

    while (k > 1 && (n << (k - 1)) > (WORD(1) << 23))
        k--;

    return k;
}

/* value of the bits [j, j + len) of e */
static ulong
extract_bits(nn_srcptr e, slong j, slong len)
{
    ulong w = 0;
    slong h;

    for (h = len - 1; h >= 0; h--)
        w = (w << 1) | TSTBIT(e, j + h);

    return w;
}

/* the window at the top bit i: the longest chain of at most k bits
   ending in a 1; returns its value and moves *i below the chain */
static ulong
next_window(nn_srcptr e, slong * i, int k, slong * nsqr)
{
    slong j = FLINT_MAX(*i - k + 1, 0);
    ulong w;

    while (TSTBIT(e, j) == 0)
        j++;

    *nsqr = *i - j + 1;
    w = extract_bits(e, j, *i - j + 1);
    *i = j - 1;
    return w;
}

/*
    r = b^e mod m via flint_mpn_mulmod_preinvn on residues shifted left
    by norm. Requires m > 1 normalized (mn >= 1, top limb nonzero),
    b < m (mn limbs), e nonzero with top limb nonzero.
    No aliasing of r with the inputs.
*/
#if FLINT_HAVE_FFT_SMALL
/*
    Drop-in for flint_mpn_mulmod_preinvn with the mullow replaced by a
    wraparound product: with the same quotient estimate
    q = X_hi + mulhigh(X_hi, dinv), the remainder r0 = X - q*d lies in
    [0, C*d) for a small C (the correction loop below is the same as
    mulmod_preinvn's), and since C*d < B^rn - 1 for rn > n it is the
    canonical value of (X - q*d) mod (B^rn - 1), with q*d computed by
    flint_mpn_mulmod_bnm1 at about 0.6 multiplications. Needs no
    precomputation beyond dinv. scratch: 5n + itch(rn) + 4rn limbs.
*/
static void
mulmod_preinvn_fold(nn_ptr r, nn_srcptr a, nn_srcptr b, mp_size_t n,
                    nn_srcptr d, nn_srcptr dinv, ulong norm,
                    mp_size_t rn, nn_ptr scratch)
{
    nn_ptr X = scratch, q = X + 2 * n, S, FX, tp;
    ulong cy;
    mp_size_t c;

    S = q + n + 1;
    FX = S + rn;
    tp = FX + rn;         /* itch(rn) limbs */

    if (a == b)
        flint_mpn_sqr(X, a, n);
    else
        flint_mpn_mul_n(X, a, b, n);

    if (norm)
        mpn_rshift(X, X, 2 * n, norm);

    /* q = X_hi + mulhigh(X_hi, dinv) < B^n */
    flint_mpn_mulhigh_n(q, X + n, dinv, n);
    mpn_add_n(q, q, X + n, n);

    /* S = q*d mod (B^rn - 1) */
    flint_mpn_mulmod_bnm1(S, rn, q, n, d, n, tp);

    /* FX = X mod (B^rn - 1) */
    {
        nn_srcptr v = X;
        mp_size_t vn = 2 * n;
        c = FLINT_MIN(vn, rn);
        flint_mpn_copyi(FX, v, c);
        flint_mpn_zero(FX + c, rn - c);
        v += c;
        vn -= c;

        while (vn > 0)
        {
            c = FLINT_MIN(vn, rn);
            cy = mpn_add(FX, FX, rn, v, c);
            while (cy)
                cy = mpn_add_1(FX, FX, rn, cy);
            v += c;
            vn -= c;
        }
    }

    /* r0 = (FX - S) mod (B^rn - 1), canonical since r0 < C*d < B^rn - 1 */
    cy = mpn_sub_n(FX, FX, S, rn);
    while (cy)
        cy = mpn_sub_1(FX, FX, rn, cy);

    FLINT_ASSERT(flint_mpn_zero_p(FX + n + 1, rn - n - 1));

    while (FX[n] != 0)
        FX[n] -= mpn_sub_n(FX, FX, d, n);

    if (mpn_cmp(FX, d, n) >= 0)
        mpn_sub_n(FX, FX, d, n);

    flint_mpn_copyi(r, FX, n);
}
#endif

static void
powm_basecase_core(nn_ptr r, nn_srcptr b, nn_srcptr e, mp_size_t en,
                   nn_srcptr msh, mp_size_t mn, nn_srcptr dinv,
                   flint_bitcnt_t norm)
{
    flint_bitcnt_t ebits = (en - 1) * FLINT_BITS + FLINT_BIT_COUNT(e[en - 1]);
    int k = powm_select_k(ebits, mn);
    slong tlen = WORD(1) << (k - 1);
    nn_ptr tab, acc, t;
    slong i, nsqr, h;
    ulong w;
    TMP_INIT;

#if FLINT_HAVE_FFT_SMALL
    mp_size_t bf_rn = 0;
    nn_ptr bf_scratch = NULL;
#endif

    TMP_START;
    acc = TMP_ALLOC((2 + tlen) * mn * sizeof(ulong));
    tab = acc + mn;
    t = tab + tlen * mn;

#if FLINT_HAVE_FFT_SMALL
    if (mn >= FLINT_MPN_POWM_BASECASE_FOLD_MIN
        && mn < FLINT_MPN_POWM_BASECASE_FOLD_MAX)
    {
        bf_rn = flint_mpn_mulmod_bnm1_next_size(mn + 1);
        bf_scratch = TMP_ALLOC((3 * mn + 1 + 2 * bf_rn
                        + flint_mpn_mulmod_bnm1_itch(bf_rn))
                        * sizeof(ulong));
    }
#define BC_MULMOD(rr, aa, bb) \
    do { \
        if (bf_rn) \
            mulmod_preinvn_fold(rr, aa, bb, mn, msh, dinv, norm, \
                                bf_rn, bf_scratch); \
        else \
            flint_mpn_mulmod_preinvn(rr, aa, bb, mn, msh, dinv, norm); \
    } while (0)
#else
#define BC_MULMOD(rr, aa, bb) \
    flint_mpn_mulmod_preinvn(rr, aa, bb, mn, msh, dinv, norm)
#endif

    /* tab[i] = b^(2i+1) << norm */
    if (norm)
        mpn_lshift(tab, b, mn, norm);
    else
        flint_mpn_copyi(tab, b, mn);
    if (tlen > 1)
    {
        BC_MULMOD(t, tab, tab);
        for (i = 1; i < tlen; i++)
            BC_MULMOD(tab + i * mn, tab + (i - 1) * mn, t);
    }

    i = ebits - 1;
    w = next_window(e, &i, k, &nsqr);
    flint_mpn_copyi(acc, tab + (w >> 1) * mn, mn);

    while (i >= 0)
    {
        if (TSTBIT(e, i) == 0)
        {
            BC_MULMOD(acc, acc, acc);
            i--;
        }
        else
        {
            w = next_window(e, &i, k, &nsqr);
            for (h = 0; h < nsqr; h++)
                BC_MULMOD(acc, acc, acc);
            BC_MULMOD(acc, acc, tab + (w >> 1) * mn);
        }
    }
#undef BC_MULMOD

    if (norm)
        mpn_rshift(r, acc, mn, norm);
    else
        flint_mpn_copyi(r, acc, mn);

    TMP_END;
}

void
_flint_mpn_powm_basecase(nn_ptr r, nn_srcptr b, nn_srcptr e, mp_size_t en,
                         nn_srcptr m, mp_size_t mn)
{
    flint_bitcnt_t norm = flint_clz(m[mn - 1]);
    nn_ptr msh, dinv;
    TMP_INIT;

    TMP_START;
    msh = TMP_ALLOC(2 * mn * sizeof(ulong));
    dinv = msh + mn;
    if (norm)
        mpn_lshift(msh, m, mn, norm);
    else
        flint_mpn_copyi(msh, m, mn);
    flint_mpn_preinvn(dinv, msh, mn);
    powm_basecase_core(r, b, e, en, msh, mn, dinv, norm);
    TMP_END;
}

/* acc = acc * s mod m for a scalar s of sl in {1, 2} limbs,
   via a full product and division; scratch: n + 3 limbs */
static void
mul_scalar_mod(nn_ptr acc, nn_srcptr s, int sl, nn_srcptr m, mp_size_t n,
               nn_ptr scratch)
{
    ulong q[3];

    if (sl == 2 && s[1] == 0)
        sl = 1;

    if (sl == 1 && s[0] == 2)
    {
        ulong cy = mpn_lshift(acc, acc, n, 1);
        if (cy != 0 || mpn_cmp(acc, m, n) >= 0)
            mpn_sub_n(acc, acc, m, n);
    }
    else
    {
        flint_mpn_mul(scratch, acc, n, s, sl);
        mpn_tdiv_qr(q, acc, 0, scratch, n + sl, m, n);
    }
}

#if FLINT_HAVE_FFT_SMALL
/*
    Montgomery reduction via a wraparound product, as in GMP's redc_n:
    with q = X_lo * minv mod B^mn the low half of q*m is exactly
    L = B^mn - X_lo, so the high half H = floor(q*m / B^mn) < m is
    recovered from S = q*m mod (B^nnc - 1) as

        H = (S - L) * B^(nnc - mn) mod (B^nnc - 1),

    a cyclic limb rotation, unique since H < m <= B^nnc - 1. The
    wraparound product runs against a cached cyclic transform of m, so
    each reduction costs one mullow plus about two thirds of a
    half-length cyclic multiplication.
*/
/*
    All-cached wraparound chain: S = x * m mod (B^rn - 1) for rn a
    power of two >= mn, splitting as in flint_mpn_mulmod_bnm1 but with
    every m-dependent quantity precomputed per modulus: the residues of
    m down the -1 spine, and for each +1 half of size h either the
    cached negacyclic transform of m mod B^h + 1 (h >= 128, b = 64,
    np = 3: 0.3-0.6 of a multiplication) or the basecase operand.
*/
typedef struct
{
    mp_size_t rn;
    slong nlev;                        /* levels with h = rn >> (i+1) >= 64 */
    sd_fft_mpn_mulmod_2expp1_ctx_struct * ctx;   /* nlev entries; h >= NEG_H used */
    sd_fft_mpn_mulmod_2expp1_scratch_struct * scr;
    double ** Fm;                      /* cached transforms, or NULL */
    nn_ptr * m2;                       /* m mod B^h + 1, h + 1 limbs (basecase levels) */
    nn_ptr mbase;                      /* m mod B^base - 1 */
    mp_size_t base;
    nn_ptr work;                       /* recursion scratch */
    nn_ptr storage;
}
powm_chain_struct;

/* residue of (xp, xn) mod B^h + 1 into (x2, h + 1 limbs, top in {0,1}) */
static void
chain_res_p1(nn_ptr x2, nn_srcptr xp, mp_size_t xn, mp_size_t h)
{
    mp_size_t lo = FLINT_MIN(xn, h), hi = xn - lo;
    ulong br, cy;

    flint_mpn_copyi(x2, xp, lo);
    flint_mpn_zero(x2 + lo, h + 1 - lo);
    if (hi > 0)
    {
        br = mpn_sub(x2, x2, h, xp + h, hi);
        if (br)
        {
            cy = mpn_add_1(x2, x2, h, 1);
            x2[h] = cy;
        }
    }
}

/* residue of (xp, xn) mod B^h - 1 into (x1, h limbs), representative */
static void
chain_res_m1(nn_ptr x1, nn_srcptr xp, mp_size_t xn, mp_size_t h)
{
    mp_size_t c = FLINT_MIN(xn, h);
    ulong cy;

    flint_mpn_copyi(x1, xp, c);
    flint_mpn_zero(x1 + c, h - c);
    xp += c; xn -= c;
    while (xn > 0)
    {
        c = FLINT_MIN(xn, h);
        cy = mpn_add(x1, x1, h, xp, c);
        while (cy)
            cy = mpn_add_1(x1, x1, h, cy);
        xp += c; xn -= c;
    }
}

/* S (s limbs) = x * m mod (B^s - 1) at chain level lev, s = rn >> lev;
   x: xn <= s limbs; work: private to this level and below */
static void
chain_mul(powm_chain_struct * H, slong lev, nn_ptr S, nn_srcptr x,
          mp_size_t xn, nn_ptr work)
{
    mp_size_t s = H->rn >> lev, h = s / 2;
    nn_ptr x1, x2, v, sb;
    ulong br, cy;
    int vtop, stop;
    mp_size_t i;

    if (s <= H->base)
    {
        /* base: plain multiplication and fold */
        mp_size_t xs = xn;
        while (xs > 0 && x[xs - 1] == 0)
            xs--;
        if (xs == 0)
        {
            flint_mpn_zero(S, s);
            return;
        }
        if (xs >= s)
            flint_mpn_mul(work, x, xs, H->mbase, s);
        else
            flint_mpn_mul(work, H->mbase, s, x, xs);
        chain_res_m1(S, work, xs + s, s);
        return;
    }

    x1 = work;              /* h */
    x2 = x1 + h;            /* h + 1 */
    v = x2 + h + 1;         /* h + 1 */
    sb = v + h + 1;         /* h + 1 */
    work = sb + h + 1;

    chain_res_m1(x1, x, xn, h);
    chain_res_p1(x2, x, xn, h);

    /* u into the low half of S */
    chain_mul(H, lev + 1, S, x1, h, work);

    /* v = x2 * (m mod B^h + 1) mod B^h + 1, in [0, B^h] */
    if (H->Fm[lev] != NULL)
        sd_fft_mpn_mulmod_2expp1_mul_cached(H->ctx + lev, v, x2,
                                            H->Fm[lev], H->scr + lev);
    else
    {
        int c = (x2[h] != 0) | ((H->m2[lev][h] != 0) << 1);
        v[h] = (ulong) flint_mpn_mulmod_2expp1_basecase(v, x2,
                    H->m2[lev], c, FLINT_BITS * h, work);
    }
    vtop = (int) v[h];

    /* s = (u - v)/2 mod B^h + 1, then S = u + s*B^h - s mod B^s - 1;
       identical to flint_mpn_mulmod_bnm1 */
    flint_mpn_copyi(sb, S, h);
    sb[h] = 0;
    br = mpn_sub(sb, sb, h + 1, v, h);
    if (vtop)
        br += mpn_sub_1(sb + h, sb + h, 1, 1);
    if (br)
    {
        cy = mpn_add_1(sb, sb, h + 1, 1);
        sb[h] += 1;
        (void) cy;
    }
    if (sb[0] & 1)
    {
        cy = mpn_add_1(sb, sb, h + 1, 1);
        sb[h] += 1;
        (void) cy;
    }
    mpn_rshift(sb, sb, h + 1, 1);
    stop = (int) sb[h];

    flint_mpn_copyi(S + h, sb, h);
    if (stop)
    {
        cy = mpn_add_1(S, S, s, 1);
        while (cy)
            cy = mpn_add_1(S, S, s, cy);
    }
    br = mpn_sub(S, S, s, sb, h + 1);
    while (br)
        br = mpn_sub_1(S, S, s, 1);

    (void) i;
}

typedef struct
{
    fft_small_plan_t P;
    fft_small_op_t Fm, FC;
    ulong nnc, zc;
    nn_ptr w, S, V;      /* zc, nnc, nnc limbs */
    int q_fft;           /* quotient via cached transform of minv */
    fft_small_plan_t PL;
    fft_small_op_t Fminv, FQ;
    int method;          /* 0: cyclic plan, 1: chain, 2: mulmod_bnm1 */
    powm_chain_struct H[1];
    mp_size_t rg;        /* method 2 fold length */
    nn_ptr tp;           /* method 2 scratch */
}
redc_fold_struct;

/* returns 1 on success; rn = smallest power of two >= max(mn, 128) */
static int
powm_chain_init(powm_chain_struct * H, nn_srcptr m, mp_size_t mn)
{
    mp_size_t rn = 128, s, h;
    slong lev, nlev;
    nn_ptr mres, p;

    while (rn < mn)
        rn *= 2;
    H->rn = rn;
    H->base = 64;

    nlev = 0;
    for (s = rn; s > H->base; s /= 2)
        nlev++;
    H->nlev = nlev;
    FLINT_ASSERT(nlev <= FLINT_MPN_POWM_CHAIN_MAX_LEVELS);

    H->ctx = flint_malloc(nlev * sizeof(sd_fft_mpn_mulmod_2expp1_ctx_struct));
    H->scr = flint_malloc(nlev * sizeof(sd_fft_mpn_mulmod_2expp1_scratch_struct));
    H->Fm = flint_malloc(nlev * sizeof(double *));
    H->m2 = flint_malloc(nlev * sizeof(nn_ptr));
    /* m residues (h+1 each level) + mbase + recursion work (~5rn) */
    H->storage = flint_malloc((3 * rn + H->base + 6 * rn + 64
                               + 8 * H->base) * sizeof(ulong));
    p = H->storage;

    mres = p;                     /* -1-spine residue of m, reused */
    p += rn;
    flint_mpn_copyi(mres, m, FLINT_MIN(mn, rn));
    if (mn < rn)
        flint_mpn_zero(mres + mn, rn - mn);

    for (lev = 0, s = rn; s > H->base; s /= 2, lev++)
    {
        h = s / 2;
        H->m2[lev] = p;
        p += h + 1;
        chain_res_p1(H->m2[lev], mres, s, h);

        if (h >= FLINT_MPN_POWM_CHAIN_NEG_H)
        {
            sd_fft_mpn_mulmod_2expp1_ctx_init(H->ctx + lev,
                get_default_mpn_ctx(), FLINT_BITS * h, h);
            sd_fft_mpn_mulmod_2expp1_scratch_init(H->scr + lev,
                H->ctx + lev);
            H->Fm[lev] = flint_aligned_alloc(4096, n_round_up(
                sd_fft_mpn_mulmod_2expp1_transformed_size(H->ctx + lev)
                * sizeof(double), 4096));
            sd_fft_mpn_mulmod_2expp1_transform(H->ctx + lev,
                H->Fm[lev], H->m2[lev], H->scr + lev);
        }
        else
            H->Fm[lev] = NULL;

        /* -1 spine: fold m's residue down to h limbs in place */
        chain_res_m1(mres, mres, s, h);
    }
    H->mbase = p;
    p += H->base;
    flint_mpn_copyi(H->mbase, mres, H->base);
    H->work = p;
    return 1;
}

static void
powm_chain_clear(powm_chain_struct * H)
{
    for (slong lev = 0; lev < H->nlev; lev++)
    {
        if (H->Fm[lev] != NULL)
        {
            flint_aligned_free(H->Fm[lev]);
            sd_fft_mpn_mulmod_2expp1_scratch_clear(H->scr + lev);
            sd_fft_mpn_mulmod_2expp1_ctx_clear(H->ctx + lev);
        }
    }
    flint_free(H->ctx);
    flint_free(H->scr);
    flint_free(H->Fm);
    flint_free(H->m2);
    flint_free(H->storage);
}

/* t = X * B^(-mn) mod m, canonical in [0, m); X: 2*mn limbs < m*B^mn;
   q: mn limbs of scratch; t must not alias X, m, minv */
static void
redc_fold(nn_ptr t, nn_srcptr X, nn_srcptr m, nn_srcptr minv, mp_size_t mn,
          redc_fold_struct * F, nn_ptr q)
{
    mp_size_t rlen = (F->method == 2) ? F->rg :
                     (F->method == 1) ? F->H->rn : (mp_size_t) F->nnc;
    mp_size_t nnc = rlen, k = rlen - mn, i;
    nn_ptr S = F->S, V = F->V;
    ulong br, cy;

    for (i = 0; i < mn; i++)
        if (X[i] != 0)
            break;
    if (i == mn)
    {
        flint_mpn_copyi(t, X + mn, mn);
        goto reduce;
    }

    if (F->q_fft)
    {
        fft_small_fft_mpn(F->FQ, X, mn, F->PL);
        fft_small_op_mul(F->FQ, F->FQ, F->Fminv, F->PL);
        fft_small_ifft(F->FQ, F->PL);
        fft_small_export_mpn(q, mn, F->FQ, F->PL);
    }
    else
        flint_mpn_mullow_n(q, X, minv, mn);

    /* S = q*m mod (B^rlen - 1), some representative < B^rlen */
    if (F->method == 2)
        flint_mpn_mulmod_bnm1(S, rlen, q, mn, m, mn, F->tp);
    else if (F->method == 1)
        chain_mul(F->H, 0, S, q, mn, F->H->work);
    else
    {
        fft_small_fft_mpn(F->FC, q, mn, F->P);
        fft_small_op_mul(F->FC, F->FC, F->Fm, F->P);
        fft_small_ifft(F->FC, F->P);
        fft_small_export_mpn(F->w, F->zc, F->FC, F->P);
        flint_mpn_zero(S, nnc);
        {
            nn_srcptr v = F->w;
            mp_size_t vn = F->zc;
            while (vn > 0)
            {
                mp_size_t c = FLINT_MIN(vn, nnc);
                cy = mpn_add(S, S, nnc, v, c);
                while (cy)
                    cy = mpn_add_1(S, S, nnc, cy);
                v += c;
                vn -= c;
            }
        }
    }

    /* U = S - L mod (B^nnc - 1) with L = B^mn - X_lo = -X_lo mod B^mn */
    mpn_neg(q, X, mn);               /* q now holds L; X_lo != 0 */
    br = mpn_sub(S, S, nnc, q, mn);
    if (br)
        mpn_sub_1(S, S, nnc, 1);     /* wrap: -B^nnc = -1 mod B^nnc - 1 */

    /* H = U * B^k mod (B^nnc - 1): rotate left by k limbs */
    if (k > 0)
    {
        flint_mpn_copyi(V, S + nnc - k, k);
        flint_mpn_copyi(V + k, S, nnc - k);
    }
    else
        flint_mpn_copyi(V, S, nnc);

    /* V = H exactly, or V = B^nnc - 1 (the H = 0 representative,
       since H <= m - 1 < B^nnc - 1). Any nonzero top limb, or an
       all-ones value when nnc == mn, forces the latter. */
    if (nnc > mn ? !flint_mpn_zero_p(V + mn, nnc - mn)
                 : V[mn - 1] == UWORD_MAX)
    {
        for (i = 0; i < nnc; i++)
            if (V[i] != UWORD_MAX)
                break;

        if (i == nnc)
        {
            flint_mpn_zero(V, nnc);  /* H = 0 */
        }
        else
        {
            FLINT_ASSERT(nnc == mn); /* a genuine H with huge top limb */
        }
    }

    cy = mpn_add_n(t, X + mn, V, mn);
    cy += mpn_add_1(t, t, mn, 1);
    if (cy)
        mpn_sub_n(t, t, m, mn);

reduce:
    if (mpn_cmp(t, m, mn) >= 0)
        mpn_sub_n(t, t, m, mn);
}

/* returns 1 and fills F on success */
static int
redc_fold_init(redc_fold_struct * F, nn_srcptr m, nn_srcptr minv,
               mp_size_t mn)
{
    ulong nnc = mn, rlen, rn_chain = 128;

    if (mn < FLINT_MPN_POWM_LARGE_FOLD_THRESHOLD)
    {
        F->method = 2;
        F->rg = flint_mpn_mulmod_bnm1_next_size(mn);
        F->tp = flint_malloc(flint_mpn_mulmod_bnm1_itch(F->rg)
                             * sizeof(ulong));
        F->nnc = 0;
        rlen = F->rg;
        F->zc = 0;
        F->w = flint_malloc(2 * rlen * sizeof(ulong));
        F->S = F->w;
        F->V = F->S + rlen;
        goto qstep;
    }

    if (!fft_small_plan_init_mpn_cyclic(F->P, get_default_mpn_ctx(),
                                        &nnc, 1))
        return 0;

    while ((mp_size_t) rn_chain < mn)
        rn_chain *= 2;

    /* measured per-limb costs: chain ~ 27, cyclic-plan fold ~ 40
       (arbitrary common unit), so the chain wins unless its
       power-of-two padding exceeds 40/27 of the cyclic length */
    F->method = (27 * rn_chain <= 33 * nnc) ? 1 : 0;
    if (F->method == 1 && !powm_chain_init(F->H, m, mn))
        F->method = 0;

    F->nnc = nnc;
    rlen = (F->method == 1) ? (ulong) F->H->rn : nnc;
    F->zc = nnc + (2 * F->P->bits + F->P->depth) / FLINT_BITS + 2;
    F->w = flint_malloc((F->zc + 2 * rlen) * sizeof(ulong));
    F->S = F->w + F->zc;
    F->V = F->S + rlen;
    if (F->method == 0)
    {
        fft_small_op_init(F->Fm, F->P);
        fft_small_op_init(F->FC, F->P);
        fft_small_fft_mpn(F->Fm, m, mn, F->P);
    }

qstep:

    F->q_fft = (mn >= FLINT_MPN_POWM_REDC_QSTEP_FFT_THRESHOLD) &&
               fft_small_plan_init_mpn(F->PL, get_default_mpn_ctx(),
                                       mn, mn, 1, 0);
    if (F->q_fft)
    {
        fft_small_op_init(F->Fminv, F->PL);
        fft_small_op_init(F->FQ, F->PL);
        fft_small_fft_mpn(F->Fminv, minv, mn, F->PL);
    }
    return 1;
}

static void
redc_fold_clear(redc_fold_struct * F)
{
    if (F->method == 2)
    {
        flint_free(F->tp);
        flint_free(F->w);
        goto qclear;
    }
    if (F->method == 1)
        powm_chain_clear(F->H);
    else
    {
        fft_small_op_clear(F->Fm);
        fft_small_op_clear(F->FC);
    }
    fft_small_plan_clear(F->P);
    flint_free(F->w);
qclear:
    if (F->q_fft)
    {
        fft_small_op_clear(F->Fminv);
        fft_small_op_clear(F->FQ);
        fft_small_plan_clear(F->PL);
    }
}
#endif

/*
    r = b^e mod m for odd m > 1, via Montgomery multiplication.
    Same conventions as _flint_mpn_powm_basecase.
*/
void
_flint_mpn_powm_redc(nn_ptr r, nn_srcptr b, nn_srcptr e, mp_size_t en,
                     nn_srcptr m, mp_size_t mn)
{
    flint_bitcnt_t ebits = (en - 1) * FLINT_BITS + FLINT_BIT_COUNT(e[en - 1]);
    int k = powm_select_k(ebits, mn);
    mp_size_t bn = mn;
    int small_base;
    nn_ptr minv, acc, X, scratch, tab = NULL;
    ulong sc[64];   /* scalar window table, 2 limbs each; k <= 6 below */
    slong i, h, nsqr, tlen;
    ulong w;
    TMP_INIT;

    FLINT_ASSERT(m[0] & 1);

    while (bn > 0 && b[bn - 1] == 0)
        bn--;

    /* scalars b^(2i+1) must fit two limbs */
    small_base = (bn == 1 && b[0] >= 2);
    if (small_base)
    {
        flint_bitcnt_t bbits = FLINT_BIT_COUNT(b[0]);
        while (k > 1 && ((UWORD(1) << k) - 1) * bbits > 2 * FLINT_BITS)
            k--;
        k = FLINT_MIN(k, 6);
    }

    TMP_START;
    minv = TMP_ALLOC(mn * sizeof(ulong));
    acc = TMP_ALLOC((mn + 3) * sizeof(ulong));
    X = TMP_ALLOC((2 * mn + 3) * sizeof(ulong));
    scratch = TMP_ALLOC(2 * mn * sizeof(ulong));

    _flint_mpn_binvert(minv, m, mn);
    mpn_neg(minv, minv, mn);

#if FLINT_HAVE_FFT_SMALL
    {
        int use_fold = (mn >= FLINT_MPN_POWM_REDC_FOLD_THRESHOLD);
        redc_fold_struct F[1];
        if (use_fold)
            use_fold = redc_fold_init(F, m, minv, mn);
#define REDC_STEP(dst, XX) \
        do { \
            if (use_fold) \
                redc_fold(dst, XX, m, minv, mn, F, scratch); \
            else \
                _flint_mpn_redc_n(dst, XX, m, minv, mn, scratch); \
        } while (0)
#else
    {
#define REDC_STEP(dst, XX) _flint_mpn_redc_n(dst, XX, m, minv, mn, scratch)
#endif

    if (small_base)
    {
        /* sc[2i], sc[2i+1] = b^(2i+1) as a two-limb value */
        ulong s0 = b[0], s1 = 0, b2hi, b2lo;

        tlen = WORD(1) << (k - 1);
        umul_ppmm(b2hi, b2lo, b[0], b[0]);
        for (i = 0; ; i++)
        {
            sc[2 * i] = s0;
            sc[2 * i + 1] = s1;
            if (i + 1 == tlen)
                break;
            /* (s1, s0) *= (b2hi, b2lo), guaranteed to fit by the k cap */
            {
                ulong t1, t0;
                umul_ppmm(t1, t0, s0, b2lo);
                t1 += s1 * b2lo + s0 * b2hi;
                s0 = t0;
                s1 = t1;
            }
        }

        /* acc = 1 * B^mn mod m */
        flint_mpn_zero(X, mn);
        X[mn] = 1;
        mpn_tdiv_qr(scratch, acc, 0, X, mn + 1, m, mn);
    }
    else
    {
        /* Montgomery table of b^(2i+1) */
        nn_ptr msq;

        tlen = WORD(1) << (k - 1);
        tab = TMP_ALLOC((tlen + 1) * mn * sizeof(ulong));
        msq = tab + tlen * mn;

        /* tab[0] = b * B^mn mod m */
        flint_mpn_zero(X, mn);
        flint_mpn_copyi(X + mn, b, mn);
        mpn_tdiv_qr(scratch, tab, 0, X, 2 * mn, m, mn);

        if (tlen > 1)
        {
            flint_mpn_sqr(X, tab, mn);
            REDC_STEP(msq, X);
            for (i = 1; i < tlen; i++)
            {
                flint_mpn_mul_n(X, tab + (i - 1) * mn, msq, mn);
                REDC_STEP(tab + i * mn, X);
            }
        }
    }

    i = ebits - 1;
    if (small_base)
    {
        /* start from mont(1) and apply the first window as a scalar */
        w = next_window(e, &i, k, &nsqr);
        mul_scalar_mod(acc, sc + 2 * (w >> 1), 2, m, mn, X);
    }
    else
    {
        w = next_window(e, &i, k, &nsqr);
        flint_mpn_copyi(acc, tab + (w >> 1) * mn, mn);
    }

    while (i >= 0)
    {
        if (TSTBIT(e, i) == 0)
        {
            flint_mpn_sqr(X, acc, mn);
            REDC_STEP(acc, X);
            i--;
        }
        else
        {
            w = next_window(e, &i, k, &nsqr);
            for (h = 0; h < nsqr; h++)
            {
                flint_mpn_sqr(X, acc, mn);
                REDC_STEP(acc, X);
            }
            if (small_base)
                mul_scalar_mod(acc, sc + 2 * (w >> 1), 2, m, mn, X);
            else
            {
                flint_mpn_mul_n(X, acc, tab + (w >> 1) * mn, mn);
                REDC_STEP(acc, X);
            }
        }
    }

    /* leave the Montgomery domain */
    flint_mpn_copyi(X, acc, mn);
    flint_mpn_zero(X + mn, mn);
    REDC_STEP(r, X);

#undef REDC_STEP
#if FLINT_HAVE_FFT_SMALL
        if (use_fold)
            redc_fold_clear(F);
    }
#else
    }
#endif

    TMP_END;
}

/* x = x^e mod 2^t2 on t2l = cdiv(t2, 64) limbs, b2 = base mod 2^t2 */
static void
powm_2exp(nn_ptr x, nn_srcptr b2, nn_srcptr e, mp_size_t en,
          flint_bitcnt_t t2, mp_size_t t2l, nn_ptr scratch)
{
    flint_bitcnt_t ebits = (en - 1) * FLINT_BITS + FLINT_BIT_COUNT(e[en - 1]);
    ulong mask = (t2 % FLINT_BITS) ?
                 ((UWORD(1) << (t2 % FLINT_BITS)) - 1) : UWORD_MAX;
    slong i;

    flint_mpn_zero(x, t2l);
    x[0] = 1;

    for (i = ebits - 1; i >= 0; i--)
    {
        /* TODO: flint_mpn_sqrlow */
        flint_mpn_mullow_n(scratch, x, x, t2l);
        flint_mpn_copyi(x, scratch, t2l);
        if (TSTBIT(e, i))
        {
            flint_mpn_mullow_n(scratch, x, b2, t2l);
            flint_mpn_copyi(x, scratch, t2l);
        }
        x[t2l - 1] &= mask;
    }
}

/*
    r = b^e mod m. r, b, m: mn limbs with m normalized and b < m;
    e: en limbs. No aliasing of r with the inputs.
*/
void
flint_mpn_powm(nn_ptr r, nn_srcptr b, nn_srcptr e, mp_size_t en,
               nn_srcptr m, mp_size_t mn)
{
    flint_bitcnt_t val2;
    mp_size_t i0, on, t2l;

    FLINT_ASSERT(mn >= 1);
    FLINT_ASSERT(m[mn - 1] != 0);

    while (en > 0 && e[en - 1] == 0)
        en--;

    if (mn == 1 && m[0] == 1)
    {
        r[0] = 0;
        return;
    }

    if (en == 0)
    {
        flint_mpn_zero(r, mn);
        r[0] = 1;
        return;
    }

    if (flint_mpn_zero_p(b + 1, mn - 1))
    {
        ulong b0 = b[0];
        /* TODO: when there is an flint_mpn_pow, more small powers could
           be special-cased here. More ambitiously, certain we could split
           powers with a small base and a relatively small exponent into
           a plain integer power followed by modular powering. */

        if (b0 == 0 || b0 == 1)
        {
            flint_mpn_zero(r, mn);
            r[0] = b0;
            return;
        }
        else if (b0 == 2 && en == 1 &&
            e[0] < (mn - 1) * FLINT_BITS + FLINT_BIT_COUNT(m[mn - 1]) - 1)
        {
            flint_mpn_zero(r, mn);
            r[e[0] / FLINT_BITS] = UWORD(1) << (e[0] % FLINT_BITS);
            return;
        }
    }

    if (en == 1 && e[0] == 1)
    {
        flint_mpn_copyi(r, b, mn);
        return;
    }

    if (mn == 1 && en == 1)
    {
        nmod_t mod;

        /* the nmod inverse is not the n = 1 flint_mpn_preinvn inverse,
           so it is computed on the fly in both entry points */
        nmod_init(&mod, m[0]);
        r[0] = nmod_pow_ui(b[0], e[0], mod);
        return;
    }

    /* Todo: for a few-limb modulus and a large exponent, GMP's
       redc_1/redc_2 assembly currently wins; fall back to mpz_powm
       until we have dedicated small-n Montgomery code. mn == 1 with a
       multi-limb exponent takes this path too. */
    if (mn <= FLINT_MPN_POWM_MPZ_FALLBACK_LIMBS)
    {
        flint_bitcnt_t ebits2 = (en - 1) * FLINT_BITS
                                + FLINT_BIT_COUNT(e[en - 1]);
        if (ebits2 > 16)
        {
            mpz_t rz;
            mpz_t bz, ez, mz;
            mp_size_t c;

            mpz_roinit_n(bz, b, mn);
            mpz_roinit_n(ez, e, en);
            mpz_roinit_n(mz, m, mn);
            mpz_init2(rz, (mn + 1) * FLINT_BITS);
            mpz_powm(rz, bz, ez, mz);
            c = FLINT_MIN((mp_size_t) mpz_size(rz), mn);
            flint_mpn_copyi(r, mpz_limbs_read(rz), c);
            flint_mpn_zero(r + c, mn - c);
            mpz_clear(rz);
            return;
        }
    }

    /* small exponents: the division-based basecase (computing its own
       preinvn inverse) beats the Montgomery setup; see
       flint_mpn_powm_preinvn for the caller-supplied-inverse variant */
    {
        flint_bitcnt_t ebits3 = (en - 1) * FLINT_BITS
                                + FLINT_BIT_COUNT(e[en - 1]);
        if (ebits3 <= FLINT_MPN_POWM_PREINVN_REDC_EBITS)
        {
            _flint_mpn_powm_basecase(r, b, e, en, m, mn);
            return;
        }
    }

    for (i0 = 0; m[i0] == 0; i0++)
        ;
    val2 = i0 * FLINT_BITS + flint_ctz(m[i0]);

    if (val2 == 0)
    {
        _flint_mpn_powm_redc(r, b, e, en, m, mn);
    }
    else
    {
        /* m = modd * 2^val2: CRT of the odd and 2-power parts */
        nn_ptr modd, bo, r2, b2, inv2, u, z, scratch;
        flint_bitcnt_t sh = val2 % FLINT_BITS;
        ulong mask = sh ? ((UWORD(1) << sh) - 1) : UWORD_MAX;
        TMP_INIT;

        on = mn - i0;
        t2l = (val2 + FLINT_BITS - 1) / FLINT_BITS;

        TMP_START;
        modd = TMP_ALLOC((3 * on + mn + 1 + 6 * t2l) * sizeof(ulong));
        bo = modd + on;
        r2 = bo + on;
        b2 = r2 + t2l;
        inv2 = b2 + t2l;
        u = inv2 + t2l;
        z = u + t2l;                       /* on + mn + 1 (>= on + t2l + 1) */
        scratch = z + on + mn + 1;         /* 2*t2l */

        if (sh)
            mpn_rshift(modd, m + i0, on, sh);
        else
            flint_mpn_copyi(modd, m + i0, on);
        if (modd[on - 1] == 0)
            on--;

        if (on == 1 && modd[0] == 1)
        {
            flint_mpn_zero(bo, 1);         /* r_odd = 0 */
            bo[0] = 0;
            flint_mpn_zero(modd, 1);
            modd[0] = 1;
            flint_mpn_zero(r2, t2l);       /* keep sizes uniform below */
            flint_mpn_zero(z, on);
            z[0] = 0;
            /* r = b^e mod 2^val2 */
            flint_mpn_copyi(b2, b, t2l);
            b2[t2l - 1] &= mask;
            powm_2exp(r2, b2, e, en, val2, t2l, scratch);
            flint_mpn_zero(r, mn);
            flint_mpn_copyi(r, r2, t2l);
            TMP_END;
            return;
        }

        /* odd part */
        {
            nn_ptr q = z;                  /* scratch for the division */
            if (mn >= on)
                mpn_tdiv_qr(q, bo, 0, b, mn, modd, on);
        }
        _flint_mpn_powm_redc(z, bo, e, en, modd, on);   /* r_odd in z */

        /* 2-power part */
        flint_mpn_copyi(b2, b, t2l);
        b2[t2l - 1] &= mask;
        powm_2exp(r2, b2, e, en, val2, t2l, scratch);

        /* u = (r2 - r_odd) * modd^(-1) mod 2^val2 */
        {
            mp_size_t c = FLINT_MIN(on, t2l);
            flint_mpn_copyi(u, z, c);
            flint_mpn_zero(u + c, t2l - c);
            mpn_sub_n(u, r2, u, t2l);

            flint_mpn_copyi(inv2, modd, c);
            flint_mpn_zero(inv2 + c, t2l - c);
            flint_mpn_copyi(b2, inv2, t2l);            /* reuse b2 */
            _flint_mpn_binvert(inv2, b2, t2l);

            flint_mpn_mullow_n(scratch, u, inv2, t2l);
            flint_mpn_copyi(u, scratch, t2l);
            u[t2l - 1] &= mask;
        }

        /* r = r_odd + modd * u  (< m) */
        {
            nn_ptr prod = TMP_ALLOC((on + t2l) * sizeof(ulong));
            ulong cy;

            if (on >= t2l)
                flint_mpn_mul(prod, modd, on, u, t2l);
            else
                flint_mpn_mul(prod, u, t2l, modd, on);
            cy = mpn_add(prod, prod, on + t2l, z, on);
            (void) cy;
            FLINT_ASSERT(cy == 0);

            flint_mpn_copyi(r, prod, mn);
            FLINT_ASSERT(on + t2l == mn ||
                flint_mpn_zero_p(prod + mn, on + t2l - mn));
        }

        TMP_END;
    }
}

/*
    r = b^e mod m with dinv the flint_mpn_preinvn inverse of m shifted
    left by norm = clz(m[mn-1]), as stored for instance in an
    fmpz_preinvn_t. Conventions as flint_mpn_powm. Small exponents run
    on the supplied inverse with no per-call precomputation; large
    exponents fall through to the Montgomery machinery, whose setup
    then amortizes.
*/
void
flint_mpn_powm_preinvn(nn_ptr r, nn_srcptr b, nn_srcptr e, mp_size_t en,
                       nn_srcptr m, mp_size_t mn, nn_srcptr dinv,
                       flint_bitcnt_t norm)
{
    flint_bitcnt_t ebits;

    FLINT_ASSERT(mn >= 1);
    FLINT_ASSERT(m[mn - 1] != 0);
    FLINT_ASSERT(norm == (flint_bitcnt_t) flint_clz(m[mn - 1]));

    while (en > 0 && e[en - 1] == 0)
        en--;

    if ((mn == 1 && m[0] == 1) || en == 0 || flint_mpn_zero_p(b, mn)
        || (en == 1 && e[0] == 1) || mn == 1)
    {
        flint_mpn_powm(r, b, e, en, m, mn);
        return;
    }

    ebits = (en - 1) * FLINT_BITS + FLINT_BIT_COUNT(e[en - 1]);

    /* for tiny moduli the mpz fallback of flint_mpn_powm takes over
       from fewer exponent bits */
    if (ebits <= (mn <= 12 ? (flint_bitcnt_t) 8
                           : (flint_bitcnt_t)
                             FLINT_MPN_POWM_PREINVN_REDC_EBITS))
    {
        nn_ptr msh;
        TMP_INIT;

        TMP_START;
        msh = TMP_ALLOC(mn * sizeof(ulong));
        if (norm)
            mpn_lshift(msh, m, mn, norm);
        else
            flint_mpn_copyi(msh, m, mn);
        powm_basecase_core(r, b, e, en, msh, mn, dinv, norm);
        TMP_END;
    }
    else
        flint_mpn_powm(r, b, e, en, m, mn);
}

