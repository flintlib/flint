/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "longlong.h"
#include "mpn_extras.h"
#include "mp_real.h"
#include "impl.h"

#if FLINT_BITS == 64

/*
    The series of the hardcoded tiers (exp_rs.c, trig_rs.c): exp, sin and
    cos, sinh and cosh, atan, atanh of a fraction x < 2^-64 resp. 2^-32.
    Every one is a sum of terms c_k v^k / D over a variable v (x for exp,
    z = x^2 for the others) with one-word coefficients c_k over a shared
    one-word denominator D (20! for the factorial families, 3 5 ... 33 for
    atan and atanh), and the powers v^j have zb j / 64 leading zero limbs
    exactly (zb = 64 resp. 32 for exp, 128 resp. 64 for z).  One engine
    evaluates them all: rectangular splitting with block length m (the
    powers v^2 .. v^m by high products of the shrinking lengths), the
    blocks Horner-combined from the top by boundary products with v^m,
    every term of a block one mpn_addmul_1 (mpn_submul_1 for the negative
    terms of the alternating series, over an accumulator that may wrap
    negative, with a branch-free running sign extension) over the window
    of limbs its power occupies, and one division by D (mpn_divrem_1, or
    a high product by a truncated inverse of D, whose top n limbs equal
    floor(B^n / D) for every n).

    This is the algorithm of the per-size straight-line routines of
    dev/gen_mp_real_exp_rs_hard.py and dev/gen_mp_real_trig_rs_hard.py,
    executing the same calls in the same order (the results are
    identical, and so are the error bounds documented there: exp at
    most 8 ulps, 11 with m = 2 from N = 5, the odd series 12, the even
    10) from the same tuned parameters (m, divrem), kept in the tables
    below, with the windows computed at run time: a schedule of the
    calls first, then a loop over it.  That costs a few hundred
    instructions per call (from twice the time of the straight-line
    routine at four limbs to 10 or 15 percent at twenty), so the
    generators emit straight-line routines only at the sizes on the
    hot paths of the elementary functions: the 2^-64 tier, to which the
    dense kernels and the larger bitwise kernels reduce (exp from 9
    limbs on, atan and atanh from 11, sin and cos from 11), and sin and
    cos at 8 to 10 limbs at
    2^-32 (the sin_cos_opt_<n> kernels, pinned to r = 32); about
    95 KB of code.  This engine (about 12 KB) serves the rest: the
    small sizes, which the per-size kernels cover themselves, sinh and
    cosh (the hyperbolic functions go through exp), and the 2^-32 tier
    of exp, atan and atanh, which the elementary functions reach only
    in rare cases of the dense reductions (the bitwise kernels called
    directly reach it with their r = 32 below 12, 18 and 19 limbs).
*/

/* exp: c_k = 20!/k!, k <= 20 */
static const ulong _rs_exp_cs[21] = {
    UWORD(0x21c3677c82b40000), UWORD(0x21c3677c82b40000), UWORD(0x10e1b3be415a0000),
    UWORD(0x05a0913f6b1e0000), UWORD(0x0168244fdac78000), UWORD(0x004807432bc18000),
    UWORD(0x000c0135dca04000), UWORD(0x0001b707b1cdc000), UWORD(0x000036e0f639b800),
    UWORD(0x00000618fee9f800), UWORD(0x0000009c197dcc00), UWORD(0x0000000e30dce400),
    UWORD(0x000000012ebd1300), UWORD(0x0000000017499f00), UWORD(0x0000000001a9d480),
    UWORD(0x00000000001c6380), UWORD(0x000000000001c638), UWORD(0x0000000000001ab8),
    UWORD(0x000000000000017c), UWORD(0x0000000000000014), UWORD(0x0000000000000001),
};

/* sin, sinh: c_k = 20!/(2k+1)!, k < 10 */
static const ulong _rs_sin_cs[10] = {
    UWORD(0x21c3677c82b40000), UWORD(0x05a0913f6b1e0000), UWORD(0x004807432bc18000),
    UWORD(0x0001b707b1cdc000), UWORD(0x00000618fee9f800), UWORD(0x0000000e30dce400),
    UWORD(0x0000000017499f00), UWORD(0x00000000001c6380), UWORD(0x0000000000001ab8),
    UWORD(0x0000000000000014),
};

/* cos, cosh: c_k = 20!/(2k)!, k <= 10 */
static const ulong _rs_cos_cs[11] = {
    UWORD(0x21c3677c82b40000), UWORD(0x10e1b3be415a0000), UWORD(0x0168244fdac78000),
    UWORD(0x000c0135dca04000), UWORD(0x000036e0f639b800), UWORD(0x0000009c197dcc00),
    UWORD(0x000000012ebd1300), UWORD(0x0000000001a9d480), UWORD(0x000000000001c638),
    UWORD(0x000000000000017c), UWORD(0x0000000000000001),
};

/* atan, atanh: c_k = D_A/(2k+1), D_A = 3 5 ... 33, k < 17 */
static const ulong _rs_atan_cs[17] = {
    UWORD(0x57e22099c030d941), UWORD(0x1d4b60334010486b), UWORD(0x1193a01ec009c50d),
    UWORD(0x0c8e04a84006fa77), UWORD(0x09c3cabbc0056d79), UWORD(0x07fd48c828be9f63),
    UWORD(0x06c2a00bd3b4fd05), UWORD(0x05dbe00a400341af), UWORD(0x052b6b5456997631),
    UWORD(0x04a01ca9c6bf33db), UWORD(0x042f56e2c002537d), UWORD(0x03d22df06c87b067),
    UWORD(0x0383ecd2f3352769), UWORD(0x03414393eaac79d3), UWORD(0x0307cc289cb23c75),
    UWORD(0x02d5befcb39e7a9f), UWORD(0x02a9c2ed62ea3521),
};

/* floor(B^23 / 20!): its top n limbs are floor(B^n / 20!) */
static const ulong _rs_fac20_inv[23] = {
    UWORD(0x0d5bbd7da5ae6428), UWORD(0x1bad98a6674c198c), UWORD(0xbc517e5f66d3a87f),
    UWORD(0x8c4034f580f18563), UWORD(0x8a29220579be7513), UWORD(0x65d073df6b79975b),
    UWORD(0x136a63b6555a5d9d), UWORD(0x2a30831afad3a3cf), UWORD(0x29e468b9e7817789),
    UWORD(0x787fc17d48fd72dd), UWORD(0x89f5a798a1fa6762), UWORD(0x6904b18f03db7880),
    UWORD(0xd5bda8a1d7ad266b), UWORD(0xa101729d2a4e36b9), UWORD(0x37b755d40a8a2944),
    UWORD(0x99e5bcf9a6e1356e), UWORD(0x5f835cb6c5bac38c), UWORD(0x685d9de972e229fe),
    UWORD(0xe06eea9ef3f3f761), UWORD(0xc80a68dcfd232c95), UWORD(0x72b4afe3c2eaeff7),
    UWORD(0x950ae900808941ea), UWORD(0x0000000000000007),
};

/* floor(B^36 / D_A): its top n limbs are floor(B^n / D_A) */
static const ulong _rs_atan_inv[36] = {
    UWORD(0xa712d9fb491679f7), UWORD(0x1555114573323f55), UWORD(0x638f51e23acaf934),
    UWORD(0x91bdd7a112df5d4a), UWORD(0xe74abe8741ad49bb), UWORD(0x6262abac0f3682b8),
    UWORD(0x1e6917fd5f2332f8), UWORD(0x6ad782597a5afa07), UWORD(0xbaf811435147449a),
    UWORD(0x997ec8b1e3a24a7f), UWORD(0x7b763f05629e51d5), UWORD(0x1879f351c84fa92a),
    UWORD(0x8400629b77e30324), UWORD(0x4a1d2ce74ce4a512), UWORD(0xec5ce71d882b3228),
    UWORD(0x6d09e3d09dbb114e), UWORD(0xb24d787c406eb993), UWORD(0x41a8a546abc3ae2a),
    UWORD(0xaa089e81fba3fb8e), UWORD(0x009567a3baa0800f), UWORD(0x2b71702a85c43ade),
    UWORD(0xe000cb20b4cbe41d), UWORD(0x8ba038d7d1c6f17d), UWORD(0xc6a39d3bfd3cdfad),
    UWORD(0x323ba627b8b33fcf), UWORD(0x2e7c85a8745a2bd0), UWORD(0x179d7ce8b1a7b2be),
    UWORD(0x342d511dcb4d462e), UWORD(0x9923359a3873b5d9), UWORD(0xcbf14d32f8f81b4f),
    UWORD(0x5f73d47adb80d262), UWORD(0xc33da2bddf6d18d4), UWORD(0x740f5ae9256bf70a),
    UWORD(0x703747bb94ab3028), UWORD(0xe9b752f717ca44ab), UWORD(0x0000000000000002),
};

/* exp at 2^-64, by N = xn: (m, divrem) */
static const unsigned char _rs_exp64_par[22][2] = {
    {0, 0}, {0, 0}, {0, 0}, {0, 0}, {2, 0}, {3, 0}, {3, 0}, {5, 0},
    {3, 0}, {3, 0}, {4, 0}, {4, 0}, {4, 0}, {4, 1}, {4, 1}, {4, 1},
    {3, 1}, {4, 1}, {4, 1}, {4, 1}, {4, 1}, {4, 1},
};

/* exp at 2^-32, by N (xn = ceil(N / 2)): (m, divrem) */
static const unsigned char _rs_exp32_par[22][2] = {
    {0, 0}, {0, 0}, {0, 0}, {0, 0}, {1, 0}, {1, 0}, {4, 0}, {4, 0},
    {6, 0}, {6, 0}, {4, 0}, {4, 0}, {4, 0}, {4, 0}, {4, 0}, {4, 0},
    {4, 0}, {4, 0}, {4, 0}, {4, 0}, {4, 0}, {4, 1},
};

/* atan at 2^-64, by xn: (m, divrem) */
static const unsigned char _rs_atan64_par[36][2] = {
    {0, 0}, {0, 0}, {0, 0}, {0, 0}, {2, 1}, {2, 1}, {2, 1}, {2, 1},
    {2, 1}, {2, 0}, {4, 0}, {4, 0}, {4, 0}, {4, 0}, {4, 0}, {6, 0},
    {4, 0}, {2, 0}, {6, 0}, {2, 0}, {6, 0}, {4, 0}, {4, 0}, {4, 0},
    {4, 0}, {4, 0}, {4, 0}, {4, 0}, {6, 0}, {4, 0}, {4, 0}, {4, 0},
    {4, 0}, {4, 0}, {4, 0}, {4, 0},
};

/* atan at 2^-32, by xn: (m, divrem) */
static const unsigned char _rs_atan32_par[18][2] = {
    {0, 0}, {0, 0}, {2, 1}, {2, 1}, {2, 1}, {4, 1}, {4, 1}, {6, 1},
    {4, 1}, {4, 1}, {6, 0}, {4, 0}, {4, 0}, {4, 0}, {4, 0}, {4, 0},
    {4, 0}, {4, 0},
};

/* atanh at 2^-64, by xn: (m, divrem) */
static const unsigned char _rs_atanh64_par[36][2] = {
    {0, 0}, {0, 0}, {0, 0}, {0, 0}, {2, 1}, {2, 1}, {2, 1}, {2, 1},
    {2, 1}, {2, 1}, {4, 1}, {4, 1}, {4, 0}, {4, 0}, {4, 0}, {6, 0},
    {4, 0}, {6, 0}, {6, 0}, {6, 0}, {6, 0}, {4, 0}, {4, 0}, {4, 0},
    {4, 0}, {4, 0}, {4, 0}, {4, 0}, {4, 0}, {4, 0}, {4, 0}, {4, 0},
    {4, 0}, {4, 0}, {4, 0}, {4, 0},
};

/* atanh at 2^-32, by xn: (m, divrem) */
static const unsigned char _rs_atanh32_par[18][2] = {
    {0, 0}, {0, 0}, {2, 1}, {2, 1}, {2, 1}, {4, 1}, {4, 1}, {6, 1},
    {4, 1}, {4, 1}, {4, 0}, {6, 0}, {4, 0}, {6, 0}, {4, 0}, {4, 0},
    {4, 0}, {4, 0},
};

/* sincos at 2^-64, by xn: (m, divrem) */
static const unsigned char _rs_sincos64_par[23][2] = {
    {0, 0}, {0, 0}, {0, 0}, {0, 0}, {2, 1}, {2, 1}, {2, 1}, {2, 1},
    {2, 1}, {4, 1}, {4, 1}, {4, 0}, {4, 0}, {6, 0}, {6, 0}, {6, 0},
    {6, 0}, {6, 0}, {8, 0}, {8, 0}, {6, 0}, {8, 0}, {8, 0},
};

/* sincos at 2^-32, by xn: (m, divrem) */
static const unsigned char _rs_sincos32_par[11][2] = {
    {0, 0}, {0, 0}, {2, 1}, {2, 1}, {2, 1}, {4, 1}, {4, 1}, {6, 1},
    {6, 1}, {8, 1}, {8, 0},
};

/* sinhcosh at 2^-64, by xn: (m, divrem) */
static const unsigned char _rs_sinhcosh64_par[23][2] = {
    {0, 0}, {0, 0}, {0, 0}, {0, 0}, {2, 1}, {2, 1}, {2, 1}, {2, 1},
    {2, 1}, {4, 1}, {4, 1}, {4, 1}, {4, 1}, {6, 0}, {6, 0}, {6, 0},
    {6, 0}, {8, 0}, {8, 0}, {8, 0}, {8, 0}, {8, 0}, {8, 0},
};

/* sinhcosh at 2^-32, by xn: (m, divrem) */
static const unsigned char _rs_sinhcosh32_par[11][2] = {
    {0, 0}, {0, 0}, {2, 1}, {2, 1}, {2, 1}, {4, 1}, {4, 1}, {6, 1},
    {6, 1}, {8, 1}, {8, 1},
};

/* the engine ***************************************************************/

#define RS_MMAX 8
#define RS_XNMAX 36
#define RS_OPMAX 40

typedef struct
{
    slong xn;           /* fraction limbs */
    ulong zb;           /* leading zero bits per power of v */
    slong m;            /* block length */
    nn_ptr pw[RS_MMAX + 1];  /* v^j at its significant limbs
                                pw[j][0 .. xn - lz(j) - 1]; v^m readable
                                one limb below (zero) */
}
_rs_struct;

/* lz(j) = the leading zero limbs of v^j */
#define RS_LZ(j) ((slong) (((ulong) (j) * zb) / FLINT_BITS))

/* the factor v^a of v^j = v^a v^(j - a), the most balanced pair
   minimizing the alignment lz(j) - lz(a) - lz(j - a): zero for zb a
   multiple of 64; at zb = 32 one for v^2 (formed at one more limb, its
   top limb zero), zero for the others by these choices */
static const unsigned char _rs_pair_a[2][RS_MMAX + 1] = {
    { 0, 0, 1, 1, 2, 2, 3, 3, 4 },
    { 0, 0, 1, 1, 2, 2, 2, 3, 4 },
};

/* the powers v^2 .. v^m that the terms k_lo .. k_hi (and the
   boundaries) need, and their factors; v^1 = pw[1] is given; buf has
   room for (m + 1) (xn + 2) limbs */
static void
_rs_powers(_rs_struct * R, slong k_lo, slong k_hi, nn_ptr buf)
{
    slong m = R->m, xn = R->xn, j, k, r, a, b, n;
    ulong zb = R->zb;
    const unsigned char * pa = _rs_pair_a[zb % FLINT_BITS != 0];
    unsigned int need;

    if (k_hi - k_lo + 1 >= m)
        need = (1u << m) - 1;
    else
        for (need = 0, k = k_lo, r = k_lo % m; k <= k_hi; k++)
        {
            need |= 1u << r;
            r = (r + 1 == m) ? 0 : r + 1;
        }
    if (k_hi >= m)
        need |= 1u << m;
    for (j = m; j >= 3; j--)
        if ((need >> j) & 1)
            need |= (1u << pa[j]) | (1u << (j - pa[j]));

    for (j = 2; j <= m; j++)
    {
        if (!((need >> j) & 1))
            continue;
        a = pa[j];
        b = j - a;
        n = xn - RS_LZ(a) - RS_LZ(b);
        if (j == m && k_hi >= m)
            *buf++ = 0;         /* the zero limb below v^m */
        R->pw[j] = buf;
        if (a == b)
            flint_mpn_sqrhigh_n(buf, R->pw[a] + (xn - RS_LZ(a) - n), n);
        else
            flint_mpn_mulhigh_n(buf, R->pw[a] + (xn - RS_LZ(a) - n),
                R->pw[b] + (xn - RS_LZ(b) - n), n);
        buf += n + 1;
    }
}

/* the operations of a sum, scheduled first (so that the loop doing the
   arithmetic carries little state across its calls) */
#define RS_OP_ADD 0     /* acc += c v^p over the window of g limbs */
#define RS_OP_SUB 1     /* acc -= c v^p */
#define RS_OP_MUL 2     /* acc = c v^p (the first term) */
#define RS_OP_NEG 3     /* acc = -c v^p */
#define RS_OP_UNITS 4   /* acc += c at limb g */
#define RS_OP_BOUNDARY 5 /* acc = acc v^m over g limbs */

typedef struct
{
    nn_srcptr src;
    ulong c;
    int op;
    int g;
}
_rs_op_struct;

/* sum_{k = k_lo}^{k_hi} (-1)^(alt k) cs[k] v^k, scaled by B^xn, by
   descending rectangular splitting into one of the buffers sa, ta
   (xn + 1 limbs each); returns the buffer, sets *occ to its top
   occupied limb.

   Every term is one mpn_addmul_1 (mpn_submul_1 for the negative terms
   of the alternating series) over the window of the g limbs its power
   occupies, with the carry or borrow into limb g.  An alternating
   accumulator may wrap negative: fl is its running sign-extension limb
   (0 or -1), copied into the limbs a wider window newly covers.  For
   the alternating series k_lo = 0 and m is even, so that every block
   ends on its dominant nonnegative units term, after which fl = 0.
   The blocks are Horner-combined from the top by high products with
   v^m over one limb more than the accumulator's (v^m padded with its
   zero limb below where that exceeds its own length). */
static nn_ptr
_rs_sum(const _rs_struct * R, const ulong * cs, slong k_lo,
    slong k_hi, int alt, nn_ptr sa, nn_ptr ta, slong * pocc)
{
    _rs_op_struct ops[RS_OPMAX], * o, * end;
    slong m = R->m, xn = R->xn, b, p, p_lo, k, g, occ, n;
    ulong zb = R->zb, fl, cy;
    nn_ptr cur, oth, t;

    /* schedule */
    o = ops;
    occ = -1;
    for (b = k_hi - k_hi % m; b + m > k_lo; b -= m)
    {
        p_lo = FLINT_MAX(k_lo - b, 0);
        for (p = FLINT_MIN(m - 1, k_hi - b); p >= FLINT_MAX(p_lo, 1); p--)
        {
            k = b + p;
            o->src = R->pw[p] + (RS_LZ(k) - RS_LZ(p));
            o->c = cs[k];
            o->op = ((occ < 0) ? RS_OP_MUL : RS_OP_ADD) + (alt & k);
            o->g = occ = xn - RS_LZ(k);
            o++;
        }
        if (p_lo == 0)
        {
            FLINT_ASSERT(!(alt & b));
            o->c = cs[b];
            o->op = RS_OP_UNITS;
            o->g = occ = xn - RS_LZ(b);
            o++;
        }
        if (b > 0)
        {
            n = occ + 1;
            o->src = (n <= xn - RS_LZ(m)) ? R->pw[m] + (xn - RS_LZ(m) - n)
                                          : R->pw[m] - 1;
            o->op = RS_OP_BOUNDARY;
            o->g = n;
            o++;
            occ = n - 1;
        }
    }
    FLINT_ASSERT(o - ops <= RS_OPMAX);
    *pocc = occ;

    /* execute */
    cur = sa;
    oth = ta;
    occ = -1;
    fl = 0;
    for (end = o, o = ops; o < end; o++)
    {
        g = o->g;
        switch (o->op)
        {
            case RS_OP_ADD:
                while (occ < g)
                    cur[++occ] = fl;
                cy = mpn_addmul_1(cur, o->src, g, o->c);
                cur[g] += cy;
                fl += (cur[g] < cy);
                break;
            case RS_OP_SUB:
                while (occ < g)
                    cur[++occ] = fl;
                cy = mpn_submul_1(cur, o->src, g, o->c);
                fl -= (cur[g] < cy);
                cur[g] -= cy;
                break;
            case RS_OP_MUL:
                cur[g] = mpn_mul_1(cur, o->src, g, o->c);
                occ = g;
                break;
            case RS_OP_NEG:
                cur[g] = mpn_mul_1(cur, o->src, g, o->c);
                fl = -(ulong) mpn_neg(cur, cur, g + 1);
                occ = g;
                break;
            case RS_OP_UNITS:
                if (occ == g)
                    cur[g] += o->c;     /* over a carry limb (zb = 32) */
                else
                {
                    while (occ < g - 1)
                        cur[++occ] = fl;
                    cur[g] = o->c + fl;
                    occ = g;
                }
                fl = 0;
                break;
            default:
                FLINT_ASSERT(fl == 0);
                flint_mpn_mulhigh_n(oth, cur, o->src, g);
                t = cur;
                cur = oth;
                oth = t;
                occ = g - 1;
        }
    }

    return cur;
}

/* floor(acc / D) of (acc, dlen): in place by mpn_divrem_1 (after
   lowering acc by one if sub1), else into oth by the truncated inverse
   (inv, ilen limbs); returns the quotient */
static nn_ptr
_rs_div(nn_ptr cur, nn_ptr oth, slong dlen, int sub1, int divrem,
    ulong D, const ulong * inv, slong ilen)
{
    if (divrem)
    {
        if (sub1)
            mpn_sub_1(cur, cur, dlen, 1);
        mpn_divrem_1(cur, 0, cur, dlen, D);
        return cur;
    }
    flint_mpn_mulhigh_n(oth, cur, inv + ilen - dlen, dlen);
    return oth;
}

/* exp **********************************************************************/

/* the zb = 32 hand routine for N = 4 (xn = 2): 1 + x + x^2/2 + x^3/6 */
FLINT_FORCE_INLINE void
_rs_exp32_4(nn_ptr res, nn_srcptr x)
{
    ulong h, l, s, h3, l3, t;

    s = x[1] * x[1];                    /* exact: x[1] < 2^32 */
    umul_ppmm(h3, l3, s, x[1]);         /* x^3 B^2 ~= h3 < 2^32 */
    (void) l3;
    umul_ppmm(h, l, x[1], x[0]);
    t = s + ((h << 1) | (l >> 63));     /* x^2 B^2 */
    t = (t >> 1) + h3 / 6;              /* x^2/2 + x^3/6 */
    add_ssaaaa(res[1], res[0], x[1], x[0], 0, t);
    res[2] = 1;
}

/* (res, xn + 1) = sum_{k < N} x^k / k! for x < 2^-zb (zb = 64: xn = N,
   x's top limb zero; zb = 32: xn = ceil(N / 2)) */
void
_mp_real_exp_series_hard(nn_ptr res, nn_srcptr x, slong N, int zb)
{
    _rs_struct R;
    ulong buf[(RS_MMAX + 1) * (RS_XNMAX + 2)];
    ulong sa[RS_XNMAX + 2], ta[RS_XNMAX + 2];
    nn_ptr acc;
    slong xn, xs, g2, i, occ, dlen, n;
    int m, dr;

    FLINT_ASSERT(N >= 1 && N <= 21);
    xn = (zb == FLINT_BITS) ? N : (N + 1) / 2;
    R.xn = xn;
    R.zb = zb;
    R.pw[1] = (nn_ptr) x;
    xs = xn - zb / FLINT_BITS;

    if (zb == 32 && N == 4)
    {
        _rs_exp32_4(res, x);
        return;
    }

    if (N <= 2)
    {
        for (i = 0; i < xn; i++)
            res[i] = (N == 2 && i < xs) ? x[i] : 0;
        res[xn] = 1;
        return;
    }

    m = (zb == 64) ? _rs_exp64_par[N][0] : _rs_exp32_par[N][0];
    dr = (zb == 64) ? _rs_exp64_par[N][1] : _rs_exp32_par[N][1];
    R.m = FLINT_MAX(m, 2);

    g2 = xn - (2 * zb) / FLINT_BITS;

    if (N == 3)
    {
        /* x^2 (over one limb more than its own at zb = 32) */
        n = xn - 2 * (zb / FLINT_BITS);
        flint_mpn_sqrhigh_n(buf, x + (xs - n), n);
        mpn_rshift(ta, buf, g2, 1);
        for (i = 0; i < xn; i++)
            res[i] = (i < g2) ? ta[i] : 0;
    }
    else
    {
        _rs_powers(&R, 3, N - 1, buf);

        /* the terms k = 3 .. N - 1 over 20!, then + x^2/2 */
        acc = _rs_sum(&R, _rs_exp_cs, 3, N - 1, 0, ta, sa, &occ);
        dlen = occ + 1;
        if (dr)
            mpn_divrem_1(res, 0, acc, dlen, _rs_exp_cs[0]);
        else
            flint_mpn_mulhigh_n(res, acc, _rs_fac20_inv + 23 - dlen, dlen);

        acc = (acc == ta) ? sa : ta;
        mpn_rshift(acc, R.pw[2], g2, 1);
        if (dlen == g2)
            res[g2] = mpn_add_n(res, res, acc, g2);
        else
            res[g2] += mpn_add_n(res, res, acc, g2);
    }

    /* + x, + 1 */
    if (xs + 1 <= xn)
    {
        res[xs] = mpn_add_n(res, res, x, xs);
        res[xn] = 1;
    }
    else
    {
        res[xn] = 1 + mpn_add_n(res, res, x, xs);
    }
}

/* the odd and even series (sin, cos, sinh, cosh, atan, atanh) **************/

#define RS_ATAN_NMAX 17
#define RS_SC_NMAX 10

/* ycos = 1 -+ z/2 (+ q, qlen limbs) */
static void
_rs_head_even(nn_ptr ycos, const _rs_struct * R, nn_srcptr q, slong qlen,
    int alt, nn_ptr th)
{
    slong xn = R->xn, gz = xn - R->zb / FLINT_BITS, i;
    ulong cy;

    mpn_rshift(th, R->pw[1], gz, 1);
    if (alt)
    {
        cy = mpn_neg(ycos, th, gz);
        for (i = gz; i < xn; i++)
            ycos[i] = -cy;
        ycos[xn] = 1 - cy;
    }
    else
    {
        flint_mpn_copyi(ycos, th, gz);
        for (i = gz; i < xn; i++)
            ycos[i] = 0;
        ycos[xn] = 1;
    }
    if (q != NULL)
    {
        cy = mpn_add_n(ycos, ycos, q, qlen);
        mpn_add_1(ycos + qlen, ycos + qlen, gz - qlen, cy);
    }
}

/* out = x * (odd series), from the accumulated sum */
static void
_rs_finish_odd(nn_ptr out, const _rs_struct * R, nn_srcptr x, int xzb,
    nn_ptr cur, nn_ptr oth, slong occ, int alt, int dr, ulong D,
    const ulong * inv, slong ilen)
{
    slong xn = R->xn, xsig = (xzb == 64) ? xn - 1 : xn, n2;
    nn_ptr q;
    ulong cy;

    if (alt)
    {
        /* the sum includes the units term; q has xn significant limbs */
        q = _rs_div(cur, oth, xn + 1, 1, dr, D, inv, ilen);
        if (xzb == 64)
        {
            flint_mpn_mulhigh_n(out, q + 1, x, xn - 1);
            out[xn - 1] = 0;
        }
        else
            flint_mpn_mulhigh_n(out, q, x, xn);
    }
    else
    {
        /* x + x (S - 1): q = (S - 1) B^xn has n2 significant limbs */
        q = _rs_div(cur, oth, occ + 1, 0, dr, D, inv, ilen);
        n2 = occ + xsig - xn;
        flint_mpn_mulhigh_n(out, q + (occ - n2), x + (xsig - n2), n2);
        cy = mpn_add_n(out, out, x, n2);
        if (xzb == 64)
            out[xn - 1] = mpn_add_1(out + n2, x + n2, xsig - n2, cy);
        else
            mpn_add_1(out + n2, x + n2, xsig - n2, cy);
    }
}

/* z = x^2 at its significant limbs */
static void
_rs_square(nn_ptr z, nn_srcptr x, slong xn, int xzb)
{
    if (xzb == 64)
        flint_mpn_sqrhigh_n(z, x + 1, xn - 2);
    else
        flint_mpn_sqrhigh_n(z, x, xn);
}

/* (res, xn) = sum_k (-+1)^k x^(2k+1) / (2k+1) for x < 2^-xzb, xzb = 64
   (x's top limb zero) or 32 */
void
_mp_real_atan_series_hard(nn_ptr res, nn_srcptr x, slong xn, int xzb, int alt)
{
    _rs_struct R;
    ulong buf[(RS_MMAX + 1) * (RS_XNMAX + 2)];
    ulong z[RS_XNMAX + 1], sa[RS_XNMAX + 2], ta[RS_XNMAX + 2];
    slong NA, i, occ;
    const unsigned char (*par)[2];
    nn_ptr cur;

    FLINT_ASSERT(xn >= 1 && xn <= ((xzb == 64) ? 35 : 17));

    NA = FLINT_MIN(FLINT_MAX(1, (xzb == FLINT_BITS) ? xn / 2 : xn), RS_ATAN_NMAX);
    if (NA == 1)
    {
        for (i = 0; i < xn; i++)
            res[i] = (xzb == 64 && i == xn - 1) ? 0 : x[i];
        return;
    }

    par = (xzb == 64) ? (alt ? _rs_atan64_par : _rs_atanh64_par)
                      : (alt ? _rs_atan32_par : _rs_atanh32_par);
    R.xn = xn;
    R.zb = 2 * xzb;
    R.m = par[xn][0];
    R.pw[1] = z;

    _rs_square(z, x, xn, xzb);
    _rs_powers(&R, alt ? 0 : 1, NA - 1, buf);
    cur = _rs_sum(&R, _rs_atan_cs, alt ? 0 : 1, NA - 1, alt, sa, ta, &occ);
    _rs_finish_odd(res, &R, x, xzb, cur, (cur == sa) ? ta : sa, occ, alt,
        !par[xn][1], _rs_atan_cs[0], _rs_atan_inv, 36);
}

/* (ysin, xn + 1), (ycos, xn + 1) = the odd and even factorial series
   (sin and cos, or sinh and cosh), either may be NULL */
void
_mp_real_sin_cos_series_hard(nn_ptr ysin, nn_ptr ycos, nn_srcptr x,
    slong xn, int xzb, int alt)
{
    _rs_struct R;
    ulong buf[(RS_MMAX + 1) * (RS_XNMAX + 2)];
    ulong z[RS_XNMAX + 1], sa[RS_XNMAX + 2], ta[RS_XNMAX + 2], th[RS_XNMAX];
    slong NS, NC, i, occ, kmax;
    const unsigned char (*par)[2];
    nn_ptr cur, q;
    int dr;

    FLINT_ASSERT(xn >= 1 && xn <= ((xzb == 64) ? 22 : 10));

    NS = FLINT_MIN(FLINT_MAX(1, (xzb == FLINT_BITS) ? xn / 2 : xn), RS_SC_NMAX);
    NC = FLINT_MIN((xzb == FLINT_BITS) ? (xn + 1) / 2 : xn, RS_SC_NMAX + 1);

    par = (xzb == 64) ? (alt ? _rs_sincos64_par : _rs_sinhcosh64_par)
                      : (alt ? _rs_sincos32_par : _rs_sinhcosh32_par);
    R.xn = xn;
    R.zb = 2 * xzb;
    R.m = FLINT_MAX(par[xn][0], 2);
    R.pw[1] = z;
    dr = !par[xn][1];

    if (NC >= 2)
        _rs_square(z, x, xn, xzb);

    /* the powers both series need */
    kmax = FLINT_MAX((NS >= 2) ? NS - 1 : 0, (NC >= 3) ? NC - 1 : 0);
    if (kmax >= 2)
        _rs_powers(&R, 0, kmax, buf);

    if (ysin != NULL)
    {
        if (NS >= 2)
        {
            cur = _rs_sum(&R, _rs_sin_cs, alt ? 0 : 1, NS - 1, alt, sa, ta, &occ);
            _rs_finish_odd(ysin, &R, x, xzb, cur, (cur == sa) ? ta : sa,
                occ, alt, dr, _rs_sin_cs[0], _rs_fac20_inv, 23);
        }
        else
        {
            for (i = 0; i < xn; i++)
                ysin[i] = (xzb == 64 && i == xn - 1) ? 0 : x[i];
        }
        ysin[xn] = 0;
    }

    if (ycos != NULL)
    {
        if (NC >= 3)
        {
            cur = _rs_sum(&R, _rs_cos_cs, 2, NC - 1, alt, sa, ta, &occ);
            q = _rs_div(cur, (cur == sa) ? ta : sa, occ + 1, 0, dr,
                _rs_cos_cs[0], _rs_fac20_inv, 23);
            _rs_head_even(ycos, &R, q, occ, alt, th);
        }
        else if (NC == 2)
        {
            _rs_head_even(ycos, &R, NULL, 0, alt, th);
        }
        else
        {
            for (i = 0; i < xn; i++)
                ycos[i] = 0;
            ycos[xn] = 1;
        }
    }
}

#endif
