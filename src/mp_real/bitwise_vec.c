/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "mpn_extras.h"
#include "machine_vectors.h"
#include "mp_real.h"
#include "thread_support.h"
#include "impl.h"

#if FLINT_BITS == 64

/* (vec4n values never cross a function boundary that is not inlined
   here; see machine_vectors.h on -Wpsabi) */
#if defined(__clang__)
# if defined(__has_warning)
#  if __has_warning("-Wpsabi")
#   pragma clang diagnostic ignored "-Wpsabi"
#  endif
# endif
#elif defined(__GNUC__)
# pragma GCC diagnostic ignored "-Wpsabi"
#endif

/*
    Bitwise reductions with signed digits, the full-precision arithmetic
    in radix 2^56 on vec4n lanes: exp, log1p, atan, and sin, cos, tan by
    the tangent half-angle reconstruction.  The exponential is described
    first; the others follow it with their own notes below.

    exp(x), 0 <= x < 1.

    REDUCTION.  With a_i = log(1 + 2^-i) and b_i = -log(1 - 2^-i),

        x = d_0 a_0 + sum_(i = 1)^(r + 1) (d_i = 1 ? a_i : -b_i) + t,

    d_0 in {0, 1} (greedy: x < 1 < a_0 + 2 a_1) and d_i in {-1, 0, 1}
    for i >= 1, then exp(x) = 2^(d_0) prod (1 + d_i 2^-i) exp(t).  The
    digits keep |x| < 2 a_i before step i: with the threshold
    0.75 2^-i, d_i = 1 if x >= 0.75 2^-i, -1 if x <= -0.75 2^-i, else
    0, which holds the invariant (a_i < 2 a_(i+1), b_i - 2 a_(i+1) <
    2^(-2i) and so on) for any approximation of x within 0.14 2^-i
    (at i = 1, where b_1 = log 2 is tightest; 0.18 2^-i at i = 2 and
    more beyond, the limit being 2 a_(i+1) - 0.75 2^-i): so the
    decisions need no exact comparison and no branch on the full
    operand.  They run on a signed 64-bit model X = x 2^P, P = i0 + 62,
    of one window of VW = 52 indices i0 <= i < i0 + 52, X updated by the
    model values floor(a_i 2^P), floor(b_i 2^P) (one unit of error per
    step, at most 54 units by the window's end with that of reading X,
    against 0.18 2^(P - i) > 2^8 units there).  Between windows the window's digits are
    applied to a head of h = ceil((r + 128) / 64) ordinary limbs of x,
    from which the next window's model is read.  Finally all of them
    are applied to the full operand in one vector pass (in passes of
    VK entries, the digits' headroom), the remainder t then lying in
    (-2^-(r+1), 2^-(r+1)); a negative t gets one more -b_(r+1), so
    that 0 <= t < b_(r+1) < 2^-r, the series contract.  About r/3 of
    the digits are nonzero (the unsigned greedy reduction uses r/2).

    RECONSTRUCTION.  exp(t) at n + 1 limbs is converted to radix 2^56
    (F = ceil(64 n / 56) + 1 digits, one guard digit) and multiplied by
    the factors 1 + d_i 2^-i in increasing order of i, each one pass
    y += d_i (y >> i): with i = 56 q + b, digit k of y >> i is
    (y_(k+q) >> b) + ((y_(k+q+1) mod 2^b) << (56 - b)), an arithmetic
    shift (emulated by a bias of 2^63) of a signed carry-save digit.
    The digits are carry-normalized (a vector pass) only when their
    tracked bound would pass 2^62, about every 64 factors.

    ERRORS.  The table entries are floors (or lower bounds within a
    few guard units) at least one limb below the output's ulp, so t
    differs from the true remainder by less than 2^-(64 n) r (above it
    through the a_i, below through the b_i); each
    reconstruction pass truncates below one unit of the guard digit,
    and the conversions and the series add their own few ulps: the
    result stays within the bound 9 r + 100 of the bitwise form (it
    measures below 10 ulps).
*/

#define VB 56
#define VM ((UWORD(1) << VB) - 1)
#define VW 52
/* entries per pass of the full-precision application: the digits stay
   below 2^62 */
#define VK ((WORD(1) << (62 - VB)) - 2)

#define P_OF(i) (VW * ((i) / VW) + 62)

/* table storage **************************************************************/

/* for i <= r: the entries a_i (and b_i) in radix 2^56 at F digits
   (floors of the values scaled by 2^(56 F), rows of stride digits), their
   top hb limbs (floors scaled by 2^(64 hb), for the head of the
   reduction) and the window models; the logarithmic tables
   a_i = log(1 + 2^-i), b_i = -log(1 - 2^-i), and the angles
   a_i = A_i = atan(2^-i) (b aliasing a).  One version of a shared cache
   (thread_support.h) per family, built at the largest F, r and hb
   requested so far by any thread (F grown geometrically, like the limb
   tables it derives from) and read by smaller requests through the top
   digits and limbs of the rows. */
typedef struct
{
    slong F, r, stride, hb;
    ulong * a;
    ulong * b;
    ulong * ah;
    ulong * bh;
    slong * amod;
    slong * bmod;
}
_vx_tab_struct;

static void
_vx_tab_free(void * version)
{
    _vx_tab_struct * T = version;

    if (T->b != T->a)
    {
        flint_free(T->b);
        flint_free(T->bh);
        flint_free(T->bmod);
    }
    flint_free(T->a);
    flint_free(T->ah);
    flint_free(T->amod);
    flint_free(T);
}

static flint_shared_cache_struct _vx_tab_cache[2] =
    { FLINT_SHARED_CACHE_INIT(_vx_tab_free),
      FLINT_SHARED_CACHE_INIT(_vx_tab_free) };

/* the table b_i = -log(1 - 2^-i) */

/* acc (nc limbs) += 2^e, e < 64 nc */
static void
_vx_addbit(nn_ptr acc, slong nc, slong e)
{
    slong q = e / FLINT_BITS;
    mpn_add_1(acc + q, acc + q, nc - q, UWORD(1) << (e % FLINT_BITS));
}

/*
    b_i = sum_(k >= 1) 2^(-ik) / k at the precision of the cached
    log(1 + 2^-i) table, for i <= rc.  For i >= s, s the larger of
    rc / 2 + 1 and min(sqrt(prec), rc + 1), by the fixed-point
    multi-summation of _mp_real_exp_logs_ensure with every sign positive
    (one reciprocal floor(2^wp / k) per odd k serving every i, the
    multiples k 2^l by shifts); below, by b_i = a_i + b_(2i)
    ((1 - u)(1 + u) = 1 - u^2), from i = s - 1 down, the table extended
    to 2 s - 2 for the purpose.  Every term is a floor, so the entries
    are lower bounds, short by fewer guard units than terms were
    summed: far below one ulp of the value limbs.
*/
/* the table (b_i, i <= R) at the layout and precision of the cached
   (a_i) = L (entries of nc = L->nc limbs, the value scaled by
   2^(64 nc), the bottom limb a guard); index 0 unused; R >= rc is
   returned */
static nn_ptr
_vx_nlogs(const _mp_real_tab_struct * L, slong rc, slong * Rp)
{
    slong nc = L->nc, R, s, i, k, l, wp, prec;
    nn_ptr B, t, scratch, num;

    FLINT_ASSERT(rc < L->num);

    prec = FLINT_BITS * (nc - 1);
    wp = FLINT_BITS * nc;
    /* direct sums from s, at least sqrt(prec) (or past rc) where they
       are cheap, and at least rc / 2 + 1: the recursion fills the rest
       with one addition per entry */
    s = FLINT_MIN(n_sqrt(prec), rc + 1);
    s = FLINT_MAX(s, rc / 2 + 1);
    s = FLINT_MAX(s, 2);
    R = FLINT_MAX(rc, 2 * s - 2);

    B = flint_calloc((R + 1) * nc, sizeof(ulong));

    t = flint_malloc((nc + 2) * sizeof(ulong));
    scratch = flint_malloc(nc * sizeof(ulong));
    num = flint_malloc((nc + 1) * sizeof(ulong));

    /* k = 2^l: 2^(-i 2^l) / 2^l, pure bits */
    for (i = s; i <= R; i++)
    {
        nn_ptr acc = B + i * nc;

        if (i < wp)
            _vx_addbit(acc, nc, wp - i);
        for (l = 1; (i << l) + l < wp; l++)
            _vx_addbit(acc, nc, wp - ((i << l) + l));
    }

    /* odd k >= 3 and their multiples k 2^l */
    for (k = 3; k * s < wp; k += 2)
    {
        flint_mpn_zero(num, nc);
        num[nc] = 1;
        mpn_divrem_1(t, 0, num, nc + 1, (ulong) k);

        for (i = s; i <= R && i * k < wp; i++)
        {
            nn_ptr acc = B + i * nc;

            _mp_real_tab_addshift(acc, t, nc, i * k, 0, scratch);
            for (l = 1; ((i * k) << l) + l < wp; l++)
                _mp_real_tab_addshift(acc, t, nc, ((i * k) << l) + l, 0, scratch);
        }
    }

    /* b_i = a_i + b_(2i) */
    for (i = s - 1; i >= 1; i--)
        mpn_add_n(B + i * nc, L->tab + i * nc, B + 2 * i * nc, nc);

    flint_free(t);
    flint_free(scratch);
    flint_free(num);

    *Rp = R;
    return B;
}

/* conversions ****************************************************************/

/* D[0..len) = bits [56 k + sh, 56 k + sh + 56) of (X, xn), zero outside,
   for sh >= -128; tmp has room for (56 len + sh + 128) / 64 + 3 limbs and
   at least xn + 2 */
static void
_vx_from_limbs(ulong * D, slong len, nn_srcptr X, slong xn, slong sh,
    nn_ptr tmp)
{
    slong k, p, tn = (VB * len + sh + 128) / FLINT_BITS + 3;

    tmp[0] = 0;
    tmp[1] = 0;
    flint_mpn_copyi(tmp + 2, X, FLINT_MIN(xn, tn - 2));
    if (tn - 2 > xn)
        flint_mpn_zero(tmp + 2 + xn, tn - 2 - xn);

    for (k = 0, p = sh + 128; k < len; k++, p += VB)
    {
        slong q = p / FLINT_BITS, b = p % FLINT_BITS;
        ulong v = (tmp[q] >> b) | ((tmp[q + 1] << 1) << (FLINT_BITS - 1 - b));
        D[k] = v & VM;
    }
}

/* R[0..rn) = bits [64 m + sh, ...) of the normalized digits (D, len),
   0 <= sh; D has two readable digits past len */
static void
_vx_to_limbs(nn_ptr R, slong rn, ulong * D, slong len, slong sh)
{
    slong m, pos;

    D[len] = 0;
    D[len + 1] = 0;
    D[len + 2] = 0;
    for (m = 0, pos = sh; m < rn; m++, pos += FLINT_BITS)
    {
        slong k = pos / VB, off = pos - VB * k;
        ulong v = (D[k] >> off) | (D[k + 1] << (VB - off));
        if (off > 2 * VB - FLINT_BITS)
            v |= D[k + 2] << (2 * VB - off);
        R[m] = v;
    }
}

/* full signed carry normalization; returns the carry out (0 or -1 for a
   value in (-1, 1)) */
static slong
_vx_normalize(ulong * D, slong len)
{
    slong c = 0, k;

    for (k = 0; k < len; k++)
    {
        slong v = (slong) D[k] + c;
        D[k] = ((ulong) v) & VM;
        c = v >> VB;
    }
    return c;
}

/* arithmetic shift right of signed lanes: ((a + 2^63) >> b) - 2^(63 - b) */
FLINT_FORCE_INLINE vec4n
_vx_sra(vec4n a, ulong b, vec4n bias, vec4n biasb)
{
    return vec4n_sub(vec4n_bit_shift_right(vec4n_add(a, bias), b), biasb);
}

/* partial signed normalization of (D, len4), len4 a multiple of 4 that
   leaves room for the carry out of the top digit: digits into
   [-2^7, 2^56 + 2^7] */
static void
_vx_pnormalize(ulong * D, slong len4)
{
    vec4n M = vec4n_set_n(VM), cprev = vec4n_zero();
    vec4n bias = vec4n_set_n(UWORD(1) << 63);
    vec4n biasb = vec4n_set_n((UWORD(1) << 63) >> VB);
    slong k;

    for (k = 0; k < len4; k += 4)
    {
        vec4n v = vec4n_load_aligned(D + k);
        vec4n c = _vx_sra(v, VB, bias, biasb);
        v = vec4n_add(vec4n_bit_and(v, M), vec4n_shift_lanes_up(c, cprev));
        vec4n_store_aligned(D + k, v);
        cprev = c;
    }
}

/* the cached tables **********************************************************/

/* the 64 bits at pos of (E, en), zero above */
static slong
_vx_bits(nn_srcptr E, slong en, slong pos)
{
    slong q = pos / FLINT_BITS, b = pos % FLINT_BITS;
    ulong v = E[q] >> b;
    if (b != 0 && q + 1 < en)
        v |= E[q + 1] << (FLINT_BITS - b);
    return (slong) v;
}

#define PAD4(m) (((m) + 3) & ~(slong) 3)

/* the rows of one table (the entries E of nc limbs, i <= rc) */
static void
_vx_tab_fill(ulong * D, ulong * H, slong * M, nn_srcptr E, slong nc,
    slong rc, slong F, slong st, slong hb, nn_ptr tmp)
{
    slong i, j;

    for (i = 0; i <= rc; i++)
    {
        nn_srcptr e = E + i * nc;
        slong pos = FLINT_BITS * nc - P_OF(i);

        /* floor(v 2^(56 F)): the bits from 64 nc - 56 F of the entry */
        _vx_from_limbs(D + i * st, F, e, nc, FLINT_BITS * nc - VB * F, tmp);
        for (j = F; j < st; j++)
            D[i * st + j] = 0;

        /* the top hb limbs */
        flint_mpn_copyi(H + i * hb, e + nc - hb, hb);

        /* floor(v 2^P), the 64 bits at pos (v 2^P < 2^62) */
        FLINT_ASSERT(pos >= 0);
        M[i] = _vx_bits(e, nc, pos);
    }
}

/* a request: F digits, i <= r, heads of hb limbs */
typedef struct
{
    slong F, r, hb;
    int atan;
}
_vx_tab_req_struct;

static int
_vx_tab_sufficient(const void * version, void * args)
{
    const _vx_tab_struct * T = version;
    const _vx_tab_req_struct * R = args;

    return R->F <= T->F && R->r <= T->r && R->hb <= T->hb;
}

/* the version for the request, grown from the current one: F
   geometrically (the limb tables grow with it), r and hb to the larger */
static void *
_vx_tab_build(const void * old, void * args)
{
    const _vx_tab_struct * O = old;
    const _vx_tab_req_struct * R = args;
    _vx_tab_struct * T;
    const _mp_real_tab_struct * L;
    slong F = R->F, rc = R->r, hb = R->hb, nc, st, nR, nv;
    int atan = R->atan;
    nn_ptr tmp, B = NULL;
    nn_srcptr E;

    if (O != NULL)
    {
        F = _mp_real_tab_grow_nc(O->F, F);
        rc = FLINT_MAX(rc, O->r);
        hb = FLINT_MAX(hb, O->hb);
    }

    /* the limb tables, at least one limb beyond both the digits and the
       heads */
    nv = FLINT_MAX((VB * F) / FLINT_BITS + 1, hb);
    L = atan ? _mp_real_atans_ensure(nv, rc) : _mp_real_exp_logs_ensure(nv, rc);
    nc = L->nc;
    E = L->tab;
    if (!atan)
        B = _vx_nlogs(L, rc, &nR);

    st = PAD4(F) + 4;
    T = flint_malloc(sizeof(_vx_tab_struct));
    T->a = flint_malloc((rc + 1) * st * sizeof(ulong));
    T->ah = flint_malloc((rc + 1) * hb * sizeof(ulong));
    T->amod = flint_malloc((rc + 1) * sizeof(slong));
    if (atan)
    {
        T->b = T->a;
        T->bh = T->ah;
        T->bmod = T->amod;
    }
    else
    {
        T->b = flint_malloc((rc + 1) * st * sizeof(ulong));
        T->bh = flint_malloc((rc + 1) * hb * sizeof(ulong));
        T->bmod = flint_malloc((rc + 1) * sizeof(slong));
    }
    T->F = F;
    T->r = rc;
    T->stride = st;
    T->hb = hb;

    tmp = flint_malloc((nc + (VB * F) / FLINT_BITS + 8) * sizeof(ulong));
    _vx_tab_fill(T->a, T->ah, T->amod, E, nc, rc, F, st, hb, tmp);
    if (!atan)
        _vx_tab_fill(T->b, T->bh, T->bmod, B, nc, rc, F, st, hb, tmp);
    flint_free(tmp);
    flint_free(B);

    return T;
}

/* the cached tables (logarithms, or angles if atan) for F digits,
   i <= rc and heads of hb limbs */
static const _vx_tab_struct *
_vx_tab_ensure(int atan, slong F, slong rc, slong hb)
{
    _vx_tab_req_struct R;

    R.F = F;
    R.r = rc;
    R.hb = hb;
    R.atan = atan;
    return flint_shared_cache_ensure(&_vx_tab_cache[atan],
        _vx_tab_sufficient, _vx_tab_build, &R);
}

/* free the radix tables (tuning and testing: no other thread may be
   using them) */
void
_mp_real_bitwise_vec_clear(void)
{
    flint_shared_cache_clear(&_vx_tab_cache[0]);
    flint_shared_cache_clear(&_vx_tab_cache[1]);
}

/* the reduction **************************************************************/

/* x (F digits, aligned, two readable past PAD4(F + 1)) -= sum of the rows
   ta[sub[j]], += sum of the rows tb[add[j]], VK rows per vector pass (the
   digits' headroom), partially normalized between passes */
static void
_vx_apply(ulong * x, slong F, const ulong * ta, const ulong * tb, slong st,
    const slong * sub, slong ns, const slong * add, slong na)
{
    slong F4 = PAD4(F), ja = 0, jb = 0, j, k;

    while (ja < ns || jb < na)
    {
        slong ea = FLINT_MIN(ns, ja + VK);
        slong eb = FLINT_MIN(na, jb + VK - (ea - ja));

        for (k = 0; k < F4; k += 4)
        {
            vec4n v = vec4n_load_aligned(x + k);
            for (j = ja; j < ea; j++)
                v = vec4n_sub(v, vec4n_load_unaligned(ta + sub[j] * st + k));
            for (j = jb; j < eb; j++)
                v = vec4n_add(v, vec4n_load_unaligned(tb + add[j] * st + k));
            vec4n_store_aligned(x + k, v);
        }
        ja = ea;
        jb = eb;

        if (ja < ns || jb < na)
        {
            _vx_pnormalize(x, PAD4(F + 1));
            /* the carry out of the top digit, back into it */
            x[F - 1] += x[F] << VB;
            x[F] = 0;
        }
    }
}

/* signed 64-bit view of bits [pos, pos + 64) of a two's complement head */
FLINT_FORCE_INLINE slong
_vx_head_bits(nn_srcptr H, slong pos)
{
    slong q = pos / FLINT_BITS, b = pos % FLINT_BITS;
    return (slong) ((H[q] >> b) | ((H[q + 1] << 1) << (FLINT_BITS - 1 - b)));
}

/* x (F digits, normalized, below 1; xl its n limbs, or an approximation
   within an ulp) becomes t; the digits 1 go to plus, -1 to minus
   (indices increasing but for the last); rt = r + 1 is the last index,
   after which t in (-2^-rt, 2^-rt), and a negative t gets one more digit
   -1 at index fix */
static void
_vx_reduce(ulong * x, nn_srcptr xl, slong n, slong F, slong rt, slong fix,
    const _vx_tab_struct * T, slong * plus, slong * np_, slong * minus,
    slong * nm_, nn_ptr head, slong h)
{
    slong st = T->stride, hb = T->hb;
    slong np = 0, nm = 0, c, i, j, k;
    nn_srcptr A = T->ah + (hb - h), B = T->bh + (hb - h);
    const ulong * ta = T->a + (T->F - F), * tb = T->b + (T->F - F);

    /* head = floor(x 2^(64 h)), h + 1 limbs (two's complement, the top
       one integral) */
    if (h <= n)
        flint_mpn_copyi(head, xl + n - h, h);
    else
    {
        flint_mpn_zero(head, h - n);
        flint_mpn_copyi(head + h - n, xl, n);
    }
    head[h] = 0;
    head[h + 1] = 0;

    /* i = 0, greedy */
    if (mpn_cmp(head, A, h) >= 0)
    {
        mpn_sub_n(head, head, A, h);
        plus[np++] = 0;
    }

    for (c = 0; VW * c <= rt; c++)
    {
        slong i0 = VW * c, i1 = FLINT_MIN(VW * (c + 1), rt + 1);
        slong P = i0 + 62, sp = np, sm = nm;
        slong X = _vx_head_bits(head, FLINT_BITS * h - P);
        const slong * am = T->amod, * bm = T->bmod;

        for (i = FLINT_MAX(i0, 1); i < i1; i++)
        {
            slong th = WORD(3) << (P - i - 2);
            slong dp = (X >= th), dm = (X <= -th);
            X -= am[i] & -dp;
            X += bm[i] & -dm;
            plus[np] = i;
            np += dp;
            minus[nm] = i;
            nm += dm;
        }

        /* the head, for the next window */
        if (i1 <= rt)
        {
            for (j = sp; j < np; j++)
                mpn_sub(head, head, h + 1, A + plus[j] * hb, h);
            for (j = sm; j < nm; j++)
                mpn_add(head, head, h + 1, B + minus[j] * hb, h);
        }
    }

    /* the full operand */
    _vx_apply(x, F, ta, tb, st, plus, np, minus, nm);

    if (_vx_normalize(x, F) < 0)
    {
        /* t in (-2^-rt, 0) */
        for (k = 0; k < F; k++)
            x[k] += tb[fix * st + k];
        _vx_normalize(x, F);
        minus[nm++] = fix;
    }

    *np_ = np;
    *nm_ = nm;
}

/* the reconstruction *********************************************************/

FLINT_FORCE_INLINE void
_vx_factor(ulong * y, slong len, slong len4, slong i, int neg, ulong * bound)
{
    slong q = i / VB, b = i % VB, k;
    vec4n mask = vec4n_set_n((UWORD(1) << b) - 1);
    vec4n bias = vec4n_set_n(UWORD(1) << 63);
    vec4n biasb = vec4n_set_n((UWORD(1) << 63) >> b);
    ulong bl = VB - b;

    /* the digit bound in units of 2^(56 - 12): y + (y >> i) grows by the
       shifted bound and one digit */
    *bound += (*bound >> b) + (UWORD(1) << 12) + 1;
    if (*bound > (UWORD(1) << (62 - VB + 12)))
    {
        /* the carry out of the integral digit (the value is positive
           and below 5, so there is none once the digits are exact, but
           a partial normalization may leave one) back into it */
        _vx_pnormalize(y, len4);
        y[len - 1] += y[len] << VB;
        y[len] = 0;
        *bound = (UWORD(1) << 12) + 2;
        *bound += (*bound >> b) + (UWORD(1) << 12) + 1;
    }

    if (neg)
    {
        for (k = 0; k + q < len; k += 4)
        {
            vec4n a = vec4n_load_unaligned(y + k + q);
            vec4n a1 = vec4n_load_unaligned(y + k + q + 1);
            vec4n t = vec4n_add(_vx_sra(a, b, bias, biasb),
                vec4n_bit_shift_left(vec4n_bit_and(a1, mask), bl));
            vec4n_store_aligned(y + k, vec4n_sub(vec4n_load_aligned(y + k), t));
        }
    }
    else
    {
        for (k = 0; k + q < len; k += 4)
        {
            vec4n a = vec4n_load_unaligned(y + k + q);
            vec4n a1 = vec4n_load_unaligned(y + k + q + 1);
            vec4n t = vec4n_add(_vx_sra(a, b, bias, biasb),
                vec4n_bit_shift_left(vec4n_bit_and(a1, mask), bl));
            vec4n_store_aligned(y + k, vec4n_add(vec4n_load_aligned(y + k), t));
        }
    }
}

/* y (len digits, zero beyond) *= prod (1 + 2^-plus[j]) (1 - 2^-minus[j]),
   in increasing order of the index */
static void
_vx_recon(ulong * y, slong len, const slong * plus, slong np,
    const slong * minus, slong nm)
{
    slong len4 = PAD4(len + 1), a = 0, b = 0;
    ulong bound = (UWORD(1) << 12) + 1;

    while (a < np || b < nm)
    {
        if (b >= nm || (a < np && plus[a] <= minus[b]))
            _vx_factor(y, len, len4, plus[a++], 0, &bound);
        else
            _vx_factor(y, len, len4, minus[b++], 1, &bound);
    }
}

/* (xr, xi) (len digits each, zero beyond, aligned, room as for
   _vx_recon) *= 1 + s i 2^-i, s = 1 or -1 (neg): xr -= s (xi >> i),
   xi += s (xr >> i), from the old values; bound as in _vx_factor, for
   both */
FLINT_FORCE_INLINE void
_vx_rot(ulong * xr, ulong * xi, slong len, slong len4, slong i, int neg,
    ulong * bound)
{
    slong q = i / VB, b = i % VB, k;
    vec4n mask = vec4n_set_n((UWORD(1) << b) - 1);
    vec4n bias = vec4n_set_n(UWORD(1) << 63);
    vec4n biasb = vec4n_set_n((UWORD(1) << 63) >> b);
    ulong bl = VB - b;

    *bound += (*bound >> b) + (UWORD(1) << 12) + 1;
    if (*bound > (UWORD(1) << (62 - VB + 12)))
    {
        /* signed values: the carry out of the top digit folded back */
        _vx_pnormalize(xr, len4);
        xr[len - 1] += xr[len] << VB;
        xr[len] = 0;
        _vx_pnormalize(xi, len4);
        xi[len - 1] += xi[len] << VB;
        xi[len] = 0;
        *bound = (UWORD(1) << 12) + 2;
        *bound += (*bound >> b) + (UWORD(1) << 12) + 1;
    }

    for (k = 0; k + q < len; k += 4)
    {
        vec4n a = vec4n_load_unaligned(xr + k + q);
        vec4n a1 = vec4n_load_unaligned(xr + k + q + 1);
        vec4n c = vec4n_load_unaligned(xi + k + q);
        vec4n c1 = vec4n_load_unaligned(xi + k + q + 1);
        vec4n tr = vec4n_add(_vx_sra(a, b, bias, biasb),
            vec4n_bit_shift_left(vec4n_bit_and(a1, mask), bl));
        vec4n ti = vec4n_add(_vx_sra(c, b, bias, biasb),
            vec4n_bit_shift_left(vec4n_bit_and(c1, mask), bl));
        vec4n vr = vec4n_load_aligned(xr + k), vi = vec4n_load_aligned(xi + k);

        if (neg)
        {
            vr = vec4n_add(vr, ti);
            vi = vec4n_sub(vi, tr);
        }
        else
        {
            vr = vec4n_sub(vr, ti);
            vi = vec4n_add(vi, tr);
        }
        vec4n_store_aligned(xr + k, vr);
        vec4n_store_aligned(xi + k, vi);
    }
}

/* (xr, xi) *= prod (1 + s i 2^-plus[j]) (1 - s i 2^-minus[j]), s = 1, or
   -1 if flip */
static void
_vx_rotate(ulong * xr, ulong * xi, slong len, const slong * plus, slong np,
    const slong * minus, slong nm, int flip)
{
    slong len4 = PAD4(len + 1), a = 0, b = 0;
    ulong bound = (UWORD(1) << 12) + 1;

    while (a < np || b < nm)
    {
        if (b >= nm || (a < np && plus[a] <= minus[b]))
            _vx_rot(xr, xi, len, len4, plus[a++], flip, &bound);
        else
            _vx_rot(xr, xi, len, len4, minus[b++], !flip, &bound);
    }
}

/* exp ************************************************************************/

/* align to 32 bytes */
#define ALIGN32(p) ((ulong *) (((uintptr_t) (p) + 31) & ~(uintptr_t) 31))

void
_mp_real_exp_bitwise_vec(nn_ptr res, nn_srcptr x, slong n, int r)
{
    slong F, len, rt, h, np, nm, q;
    slong * plus, * minus;
    ulong * X, * Y;
    nn_ptr t, y, head, tmp;
    const _vx_tab_struct * T;
    TMP_INIT;

    FLINT_ASSERT(n >= 1);
    FLINT_ASSERT(r >= 32 && r <= FLINT_BITS * n - 16);

    rt = r + 1;
    F = (FLINT_BITS * n + VB - 1) / VB + 1;
    len = F + 1;
    h = (rt + 127 + FLINT_BITS - 1) / FLINT_BITS;
    q = rt / VB + 1;

    T = _vx_tab_ensure(0, F, rt, h);

    TMP_START;
    plus = TMP_ALLOC(2 * (rt + 2) * sizeof(slong));
    minus = plus + (rt + 2);
    X = TMP_ALLOC((PAD4(F + 1) + 8 + 4) * sizeof(ulong));
    Y = TMP_ALLOC((PAD4(len + 1) + q + 16 + 4) * sizeof(ulong));
    X = ALIGN32(X);
    Y = ALIGN32(Y);
    t = TMP_ALLOC((n + n + 1 + h + 2) * sizeof(ulong));
    y = t + n;
    head = y + n + 1;
    tmp = TMP_ALLOC(((VB * len + 128) / FLINT_BITS + n + 8) * sizeof(ulong));

    /* x -> F digits (floor(x 2^(56 F)), exact) */
    _vx_from_limbs(X, F, x, n, FLINT_BITS * n - VB * F, tmp);
    flint_mpn_zero(X + F, PAD4(F + 1) + 8 - F);

    /* t in [0, b_rt), b_rt = -log(1 - 2^-rt) < 2^-r */
    _vx_reduce(X, x, n, F, rt, rt, T, plus, &np, minus, &nm, head, h);

    /* t < 2^-r */
    _vx_to_limbs(t, n, X, F, VB * F - FLINT_BITS * n);
    _mp_real_exp_reduced(y, NULL, t, n, r, 0);

    /* y (n + 1 limbs) -> F fraction digits and an integral one */
    _vx_from_limbs(Y, len, y, n + 1, FLINT_BITS * n - VB * F, tmp);
    flint_mpn_zero(Y + len, PAD4(len + 1) + q + 16 - len);

    _vx_recon(Y, len, plus, np, minus, nm);

    _vx_normalize(Y, len);
    _vx_to_limbs(res, n + 1, Y, len, VB * F - FLINT_BITS * n);

    TMP_END;
}

/* log1p **********************************************************************/

/*
    log(1 + x), 0 <= x < 1.  The digits drive Y = 1 + x to 1 by factors
    1 + d_i 2^-i, i = 1..r: with w = -log Y, the step w -= log(1 + d_i 2^-i)
    is that of the exponential's reduction, so the thresholds carry over,
    transferred to u = Y - 1 (d_i = -1 if u >= c + c^2/2, 1 if
    u <= -(c - c^2/2), c = 0.75 2^-i: the thresholds on w to O(c^3), 0.009
    at i = 1 against the tolerance 0.07 there).  No step at i = 0:
    log 2 < 2 a_1.  The decisions run on U = u 2^P, P = i0 + 61 (|u| < 1 in
    the first window, below 2^(1-i0) (1 + 2^(1-i0)) later), updated by
    U += d_i (2^(P-i) + (U >> i)) with one unit of error per step; between
    windows the factors are applied to a head of Y in ordinary limbs.
    Afterwards |w| < 2 a_(r+1) < 2^-r, and Y (full precision) is formed
    by one vector pass per factor.  The residual is
    log(1 + u) = +-2 atanh(|u| / (2 + u)), the quotient below 2^-r, and

        log(1 + x) = log(1 + u) - sum_(d_i = 1) a_i + sum_(d_i = -1) b_i,

    the table sum in one batched vector pass.  The radix arithmetic
    contributes far below an ulp (each pass truncates below a unit of the
    guard digit, 2^-56 ulp); the conversions, the division, the doubled
    atanh series and the final rounding stay within 40 ulps, inside the
    3 r + 64 of _mp_real_log1p_bitwise_rs.  A result rounded below zero
    (x near 0) is clamped to zero.
*/

/* (head, hn) *= 1 + 2^-i resp. 1 - 2^-i (truncated), s scratch of hn limbs */
FLINT_FORCE_INLINE void
_vx_head_factor(nn_ptr head, slong hn, slong i, int neg, nn_ptr s)
{
    slong q = i / FLINT_BITS, m = hn - q;
    int b = (int) (i % FLINT_BITS);

    if (m <= 0)
        return;
    if (b != 0)
        mpn_rshift(s, head + q, m, b);
    else
        flint_mpn_copyi(s, head + q, m);
    if (neg)
        mpn_sub(head, head, hn, s, m);
    else
        mpn_add(head, head, hn, s, m);
}

/* the digits driving Y = (head, h + 1 limbs, the top one integral, Y in
   [1, 2)) to 1 by factors 1 + d_i 2^-i, i = 1..r: d_i = 1 to plus, -1 to
   minus (indices increasing); the head is consumed */
static void
_vx_log1p_digits(nn_ptr head, slong h, slong r, slong * plus, slong * np_,
    slong * minus, slong * nm_, nn_ptr scratch)
{
    slong np = 0, nm = 0, c, i, j;

    for (c = 0; VW * c <= r; c++)
    {
        slong i0 = VW * c, i1 = FLINT_MIN(VW * (c + 1), r + 1);
        slong P = i0 + 61, sp = np, sm = nm;
        /* U = (Y - 1) 2^P, the bits of Y 2^P modulo 2^64 less 2^P */
        slong U = _vx_head_bits(head, FLINT_BITS * h - P);
        /* the quadratic term below matters in the first window only
           (where i < VW: the shifts are in range) */
        int first = (c == 0);

        if (P < FLINT_BITS)
            U -= WORD(1) << P;

        for (i = FLINT_MAX(i0, 1); i < i1; i++)
        {
            /* u >= c + c^2 / 2 resp. u <= -(c - c^2 / 2), c = 0.75 2^-i:
               w = -log(1 + u) <= -c resp. >= c to O(c^3) */
            slong th = WORD(3) << (P - i - 2);
            slong sq = first ? (slong) ((UWORD(9) << (P - i - 5)) >> i) : 0;
            slong v = (WORD(1) << (P - i)) + (U >> FLINT_MIN(i, FLINT_BITS - 1));
            slong dp = (U <= sq - th), dm = (U >= th + sq);
            U += (v & -dp) - (v & -dm);
            plus[np] = i;
            np += dp;
            minus[nm] = i;
            nm += dm;
        }

        /* the head, for the next window */
        if (i1 <= r)
        {
            for (j = sp; j < np; j++)
                _vx_head_factor(head, h + 1, plus[j], 0, scratch);
            for (j = sm; j < nm; j++)
                _vx_head_factor(head, h + 1, minus[j], 1, scratch);
        }
    }

    *np_ = np;
    *nm_ = nm;
}

void
_mp_real_log1p_bitwise_vec(nn_ptr res, nn_srcptr x, slong n, int r)
{
    slong F, len, h, np, nm, q, wn, cz;
    slong * plus, * minus;
    ulong * Y, * Z;
    nn_ptr u, S, D, t, s, head, tmp;
    const _vx_tab_struct * T;
    const ulong * ta, * tb;
    ulong hi;
    int neg;
    TMP_INIT;

    FLINT_ASSERT(n >= 1);
    FLINT_ASSERT(r >= 32 && r <= FLINT_BITS * n - 16);

    F = (FLINT_BITS * n + VB - 1) / VB + 1;
    len = F + 1;
    h = (r + 127 + FLINT_BITS - 1) / FLINT_BITS;
    q = r / VB + 1;
    wn = n + 1;

    T = _vx_tab_ensure(0, F, r, 1);
    ta = T->a + (T->F - F);
    tb = T->b + (T->F - F);

    TMP_START;
    plus = TMP_ALLOC(2 * (r + 2) * sizeof(slong));
    minus = plus + (r + 2);
    Y = TMP_ALLOC((PAD4(len + 1) + q + 16 + 4) * sizeof(ulong));
    Z = TMP_ALLOC((PAD4(F + 1) + 8 + 4) * sizeof(ulong));
    Y = ALIGN32(Y);
    Z = ALIGN32(Z);
    u = TMP_ALLOC((wn + wn + n + n + n + 2 * h + 4) * sizeof(ulong));
    S = u + wn;
    D = S + wn;
    t = D + n;
    s = t + n;
    head = s + n;
    tmp = TMP_ALLOC(((VB * len + 128) / FLINT_BITS + n + 8) * sizeof(ulong));

    /* the digits, from a head of Y = 1 + x */
    if (h <= n)
        flint_mpn_copyi(head, x + n - h, h);
    else
    {
        flint_mpn_zero(head, h - n);
        flint_mpn_copyi(head + h - n, x, n);
    }
    head[h] = 1;
    _vx_log1p_digits(head, h, r, plus, &np, minus, &nm, tmp);

    /* Y = (1 + x) prod (1 + d_i 2^-i) = 1 + u, |u| < 2^-r (1 + 2^-r):
       F fraction digits and an integral one */
    _vx_from_limbs(Y, F, x, n, FLINT_BITS * n - VB * F, tmp);
    Y[F] = 1;
    flint_mpn_zero(Y + len, PAD4(len + 1) + q + 16 - len);
    _vx_recon(Y, len, plus, np, minus, nm);
    _vx_normalize(Y, len);
    _vx_to_limbs(u, wn, Y, len, VB * F - FLINT_BITS * n);

    /* D = |u|, S = 1 + Y = 2 + u */
    neg = (u[n] == 0);
    if (neg)
        mpn_neg(D, u, n);
    else
        flint_mpn_copyi(D, u, n);
    flint_mpn_copyi(S, u, wn);
    S[n] += 1;

    /* |u| < 2^(1-r); t = |u| / (2 + u) < 2^-r */
    FLINT_ASSERT(flint_mpn_zero_p(D + n - (r - 1) / FLINT_BITS, (r - 1) / FLINT_BITS));
    _mp_real_div_fraction(t, D, S, n);

    /* log(1 + u) = +- 2 atanh(t) */
    _mp_real_atanh_rs(s, NULL, t, n);
    mpn_lshift(s, s, n, 1);

    /* log(1 + x) = log(1 + u) - sum_(d_i = 1) a_i + sum_(d_i = -1) b_i */
    flint_mpn_zero(Z, PAD4(F + 1) + 8);
    _vx_apply(Z, F, ta, tb, T->stride, plus, np, minus, nm);
    cz = _vx_normalize(Z, F);
    _vx_to_limbs(res, n, Z, F, VB * F - FLINT_BITS * n);

    hi = (ulong) cz;
    if (neg)
        hi -= mpn_sub_n(res, res, s, n);
    else
        hi += mpn_add_n(res, res, s, n);

    /* below zero only by the rounding errors of a result near 0 */
    if (hi != 0)
    {
        FLINT_ASSERT((slong) hi < 0);
        flint_mpn_zero(res, n);
    }

    TMP_END;
}

/* atan ***********************************************************************/

/*
    atan(x), 0 <= x < 1, by signed vectoring: (X, Y) = (1, x) is rotated
    by factors 1 - d_i i 2^-i, i = 1..r, d_i = +-1 when
    +-Y >= 0.75 2^-i X.  On the angle theta of (X, Y) this is the
    signed-digit reduction with the steps A_i = atan(2^-i) on both sides:
    |theta| < 2 A_i before step i holds (A_i < 2 A_(i+1)) for any threshold
    in [A_i - 2 A_(i+1), 2 A_(i+1)] (the lower end is negative), so the
    decisions tolerate about 0.23 2^-i, and they run on the models
    XM = X 2^62 (X < 1.65: the vector grows by at most
    sqrt(2) prod sqrt(1 + 4^-i)) and YM = Y 2^P, P = i0 + 61, per
    window, the heads of X and Y (two's complement limbs) rotated between
    windows.  Afterwards |theta| < 2^-r and X lies in [1, 1.65]: one
    division gives the residual |Y| / X for the atan series, and

        atan(x) = sum_(d_i = 1) A_i - sum_(d_i = -1) A_i +- atan(|Y| / X).

    Errors as for log1p, within the 4 r + 64 of _mp_real_atan_bitwise_rs.
*/

/* (S, hn) = (H, hn) >> i, arithmetic (two's complement) */
static void
_vx_head_sra(nn_ptr S, nn_srcptr H, slong hn, slong i)
{
    slong q = i / FLINT_BITS, m = hn - q, j;
    int b = (int) (i % FLINT_BITS);
    ulong ext = ((slong) H[hn - 1] < 0) ? ~UWORD(0) : 0;

    if (m <= 0)
    {
        for (j = 0; j < hn; j++)
            S[j] = ext;
        return;
    }
    if (b != 0)
    {
        mpn_rshift(S, H + q, m, b);
        S[m - 1] |= ext << (FLINT_BITS - b);
    }
    else
        flint_mpn_copyi(S, H + q, m);
    for (j = m; j < hn; j++)
        S[j] = ext;
}

/* the digits rotating (X, Y) = (hx, hy) (h + 1 limbs each, two's
   complement, the top one integral; initially (1, x)) onto the real
   axis by factors 1 - d_i i 2^-i, i = 1..r: d_i = 1 to plus, -1 to minus
   (indices increasing); the heads are consumed, s1 and s2 scratch of
   h + 1 limbs */
static void
_vx_atan_digits(nn_ptr hx, nn_ptr hy, slong h, slong r, slong * plus,
    slong * np_, slong * minus, slong * nm_, nn_ptr s1, nn_ptr s2)
{
    slong np = 0, nm = 0, c, i, j, hn = h + 1;

    for (c = 0; VW * c <= r; c++)
    {
        slong i0 = VW * c, i1 = FLINT_MIN(VW * (c + 1), r + 1);
        slong P = i0 + 61, sp = np, sm = nm;
        /* XM = X 2^62 (X < 2), YM = Y 2^P */
        slong XM = _vx_head_bits(hx, FLINT_BITS * h - 62);
        slong YM = _vx_head_bits(hy, FLINT_BITS * h - P);

        for (i = FLINT_MAX(i0, 1); i < i1; i++)
        {
            /* Y >= 0.75 X 2^-i resp. <= -0.75 X 2^-i */
            slong e = 63 + i - P;
            slong th = (XM >> e) + (XM >> (e + 1));
            slong vy = XM >> (e - 1);
            slong vx = YM >> FLINT_MIN(P + i - 62, FLINT_BITS - 1);
            slong dp = (YM >= th), dm = (YM <= -th);
            YM -= (vy & -dp) - (vy & -dm);
            XM += (vx & -dp) - (vx & -dm);
            plus[np] = i;
            np += dp;
            minus[nm] = i;
            nm += dm;
        }

        /* the heads, for the next window */
        if (i1 <= r)
        {
            for (j = sp; j < np; j++)
            {
                _vx_head_sra(s1, hx, hn, plus[j]);
                _vx_head_sra(s2, hy, hn, plus[j]);
                mpn_add_n(hx, hx, s2, hn);
                mpn_sub_n(hy, hy, s1, hn);
            }
            for (j = sm; j < nm; j++)
            {
                _vx_head_sra(s1, hx, hn, minus[j]);
                _vx_head_sra(s2, hy, hn, minus[j]);
                mpn_sub_n(hx, hx, s2, hn);
                mpn_add_n(hy, hy, s1, hn);
            }
        }
    }

    *np_ = np;
    *nm_ = nm;
}

void
_mp_real_atan_bitwise_vec(nn_ptr res, nn_srcptr x, slong n, int r)
{
    slong F, len, h, np, nm, q, wn, cz, cy;
    slong * plus, * minus;
    ulong * XR, * XI, * Z;
    nn_ptr Xl, u, D, t, y, hx, hy, s1, s2, tmp;
    const _vx_tab_struct * T;
    const ulong * ta;
    ulong hi;
    int neg;
    TMP_INIT;

    FLINT_ASSERT(n >= 1);
    FLINT_ASSERT(r >= 32 && r <= FLINT_BITS * n - 16);

    F = (FLINT_BITS * n + VB - 1) / VB + 1;
    len = F + 1;
    h = (r + 127 + FLINT_BITS - 1) / FLINT_BITS;
    q = r / VB + 1;
    wn = n + 1;

    T = _vx_tab_ensure(1, F, r, 1);
    ta = T->a + (T->F - F);

    TMP_START;
    plus = TMP_ALLOC(2 * (r + 2) * sizeof(slong));
    minus = plus + (r + 2);
    XR = TMP_ALLOC((PAD4(len + 1) + q + 16 + 4) * sizeof(ulong));
    XI = TMP_ALLOC((PAD4(len + 1) + q + 16 + 4) * sizeof(ulong));
    Z = TMP_ALLOC((PAD4(F + 1) + 8 + 4) * sizeof(ulong));
    XR = ALIGN32(XR);
    XI = ALIGN32(XI);
    Z = ALIGN32(Z);
    Xl = TMP_ALLOC((wn + wn + n + n + n + 4 * (h + 2)) * sizeof(ulong));
    u = Xl + wn;
    D = u + wn;
    t = D + n;
    y = t + n;
    hx = y + n;
    hy = hx + (h + 2);
    s1 = hy + (h + 2);
    s2 = s1 + (h + 2);
    tmp = TMP_ALLOC(((VB * len + 128) / FLINT_BITS + n + 8) * sizeof(ulong));

    /* the digits, from heads of (1, x) */
    flint_mpn_zero(hx, h + 2);
    hx[h] = 1;
    if (h <= n)
        flint_mpn_copyi(hy, x + n - h, h);
    else
    {
        flint_mpn_zero(hy, h - n);
        flint_mpn_copyi(hy + h - n, x, n);
    }
    hy[h] = 0;
    hy[h + 1] = 0;
    _vx_atan_digits(hx, hy, h, r, plus, &np, minus, &nm, s1, s2);

    /* (X, Y) = (1, x) prod (1 - d_i i 2^-i), the angle of (X, Y) below
       2^-r in magnitude: F fraction digits and an integral one each */
    flint_mpn_zero(XR, PAD4(len + 1) + q + 16);
    XR[F] = 1;
    _vx_from_limbs(XI, F, x, n, FLINT_BITS * n - VB * F, tmp);
    flint_mpn_zero(XI + F, PAD4(len + 1) + q + 16 - F);
    _vx_rotate(XR, XI, len, plus, np, minus, nm, 1);
    _vx_normalize(XR, len);
    cy = _vx_normalize(XI, len);
    _vx_to_limbs(Xl, wn, XR, len, VB * F - FLINT_BITS * n);
    _vx_to_limbs(u, n, XI, len, VB * F - FLINT_BITS * n);

    /* D = |Y|; X in [1, 1.65], since each factor adds |Y| 2^-i to X, up
       to the radix passes' truncations (units of the guard digit), which
       could leave it a hair below 1 if all digits are deep */
    neg = (cy < 0);
    if (neg)
        mpn_neg(D, u, n);
    else
        flint_mpn_copyi(D, u, n);
    if (Xl[n] == 0)
    {
        flint_mpn_zero(Xl, n);
        Xl[n] = 1;
    }
    FLINT_ASSERT(Xl[n] == 1);

    /* t = |Y| / X < 2^-r */
    FLINT_ASSERT(flint_mpn_zero_p(D + n - (r - 1) / FLINT_BITS, (r - 1) / FLINT_BITS));
    _mp_real_div_fraction(t, D, Xl, n);

    /* the angle of (X, Y): +- atan(t) */
    _mp_real_atan_rs(y, NULL, t, n);

    /* atan(x) = sum_(d_i = 1) A_i - sum_(d_i = -1) A_i +- atan(t) */
    flint_mpn_zero(Z, PAD4(F + 1) + 8);
    _vx_apply(Z, F, ta, ta, T->stride, minus, nm, plus, np);
    cz = _vx_normalize(Z, F);
    _vx_to_limbs(res, n, Z, F, VB * F - FLINT_BITS * n);

    hi = (ulong) cz;
    if (neg)
        hi -= mpn_sub_n(res, res, y, n);
    else
        hi += mpn_add_n(res, res, y, n);

    /* below zero only by the rounding errors of a result near 0 */
    if (hi != 0)
    {
        FLINT_ASSERT((slong) hi < 0);
        flint_mpn_zero(res, n);
    }

    TMP_END;
}

/* sin, cos, tan **************************************************************/

/*
    The reduction of _mp_real_tan_halfangle_mid with signed digits:
    x / 2 (exact in radix 2^56) is reduced as the exponential's argument,
    with the angle tables in both directions (threshold 0.75 2^-i,
    tolerance 0.23 2^-i as for atan; A_0 = pi/4 is never taken, x/2 < 1/2),
    to x / 2 = sum d_i A_i + t', |t'| < 2 A_(r+2); a negative t' receives
    one more A_r (A_(r+1) would not do: A_i - 2 A_(i+1) < 0), so that
    0 <= t' < A_r < 2^-r.  The residual vector, (1, tan t') or
    (cos t', sin t'), is rotated by the factors 1 + d_i i 2^-i to a
    positive multiple (DE, TT) of exp(i x / 2): DE in (0.87, 1.17) and
    TT < 0.56, the ranges _mp_real_tan_halfangle_finish expects (the
    rotation's growth is at most prod sqrt(1 + 4^-i) < 1.1645, besides
    1 / cos t' for the tangent).
*/

void
_mp_real_tan_halfangle_vec(nn_ptr ysin, nn_ptr ycos, nn_ptr ytan,
    nn_srcptr x, slong n, int r)
{
    slong F, len, rt, h, np, nm, q, wn;
    slong * plus, * minus;
    ulong * X, * XR, * XI;
    nn_ptr xh, t, TT, DE, ss, cc, head, tmp;
    const _vx_tab_struct * T;
    TMP_INIT;

    FLINT_ASSERT(n >= 1);
    FLINT_ASSERT(r >= 32 && r <= FLINT_BITS * n - 16);

    rt = r + 1;
    F = (FLINT_BITS * n + VB - 1) / VB + 1;
    len = F + 1;
    h = (rt + 127 + FLINT_BITS - 1) / FLINT_BITS;
    q = rt / VB + 1;
    wn = n + 1;

    T = _vx_tab_ensure(1, F, rt, h);

    TMP_START;
    plus = TMP_ALLOC(2 * (rt + 2) * sizeof(slong));
    minus = plus + (rt + 2);
    X = TMP_ALLOC((PAD4(F + 1) + 8 + 4) * sizeof(ulong));
    XR = TMP_ALLOC((PAD4(len + 1) + q + 16 + 4) * sizeof(ulong));
    XI = TMP_ALLOC((PAD4(len + 1) + q + 16 + 4) * sizeof(ulong));
    X = ALIGN32(X);
    XR = ALIGN32(XR);
    XI = ALIGN32(XI);
    xh = TMP_ALLOC((n + n + wn + wn + (n + 2) + (n + 2) + h + 2)
        * sizeof(ulong));
    t = xh + n;
    TT = t + n;
    DE = TT + wn;
    ss = DE + wn;
    cc = ss + (n + 2);
    head = cc + (n + 2);
    tmp = TMP_ALLOC(((VB * len + 128) / FLINT_BITS + n + 8) * sizeof(ulong));

    /* x / 2 -> F digits (exact), and to n limbs for the head */
    _vx_from_limbs(X, F, x, n, FLINT_BITS * n - VB * F + 1, tmp);
    flint_mpn_zero(X + F, PAD4(F + 1) + 8 - F);
    mpn_rshift(xh, x, n, 1);

    /* x / 2 = sum d_i A_i + t', t' in [0, A_r) (a negative remainder,
       above -2 A_(r+2), gets A_r more) */
    _vx_reduce(X, xh, n, F, rt, r, T, plus, &np, minus, &nm, head, h);
    _vx_to_limbs(t, n, X, F, VB * F - FLINT_BITS * n);

    /* (DE, TT): a positive multiple of exp(i t') */
    if (n >= MP_REAL_TAN_TAPERED_MIN && (n <= MP_REAL_TAN_TAPERED_MAX
            || (n <= MP_REAL_TAN_RS_MAX
                && _mp_real_series_rs_tan_ok(n, (flint_bitcnt_t) r))))
    {
        /* (1, tan t') */
        _mp_real_series_rs_tan(TT, t, n, (flint_bitcnt_t) r);
        TT[n] = 0;
        flint_mpn_zero(DE, n);
        DE[n] = 1;
    }
    else
    {
        /* (1 - g, sin t'), g = 1 - cos t' */
        _mp_real_sin_cos_reduced(ss, cc, NULL, t, n, (flint_bitcnt_t) r, 0);
        flint_mpn_copyi(TT, ss, n);
        TT[n] = 0;
        if (flint_mpn_zero_p(cc, n))
        {
            flint_mpn_zero(DE, n);
            DE[n] = 1;
        }
        else
        {
            mpn_neg(DE, cc, n);
            DE[n] = 0;
        }
    }

    /* times W = prod (1 + d_i i 2^-i) */
    _vx_from_limbs(XR, len, DE, wn, FLINT_BITS * n - VB * F, tmp);
    flint_mpn_zero(XR + len, PAD4(len + 1) + q + 16 - len);
    _vx_from_limbs(XI, len, TT, wn, FLINT_BITS * n - VB * F, tmp);
    flint_mpn_zero(XI + len, PAD4(len + 1) + q + 16 - len);
    _vx_rotate(XR, XI, len, plus, np, minus, nm, 0);
    _vx_normalize(XR, len);
    _vx_normalize(XI, len);
    _vx_to_limbs(DE, wn, XR, len, VB * F - FLINT_BITS * n);
    _vx_to_limbs(TT, wn, XI, len, VB * F - FLINT_BITS * n);
    FLINT_ASSERT(TT[n] == 0);

    _mp_real_tan_halfangle_finish(ysin, ycos, ytan, TT, DE, n);

    TMP_END;
}

#endif
