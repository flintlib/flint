/*
    Copyright (C) 2026 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdint.h>
#include <string.h>
#include "longlong.h"
#include "mpn_extras.h"
#include "nmod.h"
#include "nmod_vec.h"
#include "nmod_poly.h"
#include "impl.h"

/*
    Karatsuba with AVX2 (see impl.h for when it is used): coefficients nlo
    to nhi - 1 of the product, for any window.

    Windows. Full and low products (nlo = 0) are full products of the
    operands truncated to nhi. Other windows are middle products, computed
    by the transposed Karatsuba (see _mid_rec).

    Blocks. The block length is the shorter length (min(len2, m) for a
    middle product of m outputs) padded to a multiple of 16. A product cuts
    the longer operand into blocks; a middle product cuts its outputs, or
    with few outputs the shorter operand, into blocks. A short last block
    is computed with shorter blocks over the other dimension when that is
    cheaper (_split_last_block). The recursion halves the number of blocks
    of 16, an odd number b into ceil(b/2) and floor(b/2).

    Leaves (16 x 16 coefficients). A product leaf broadcasts one
    coefficient of b and multiplies it by windows of 4 coefficients of a,
    into 8 accumulators that hold the 31 output coefficients; the windows
    that leave the block read zeros, as the operands of products are copied
    into blocks of 16 with 4 zeros on each side ("gapped"). A middle-product
    leaf has 4 accumulators. Nodes of 32 coefficients combine their three
    leaves in registers.

    Arithmetic. In one word, all sums and products are taken modulo 2^64
    without reduction (vec4n_mul, vec4n_mul_signed for the middle products,
    whose operands are signed), and each output coefficient is reduced once
    at the end: this is exact under the conditions of
    _nmod_poly_karatsuba_one_word. Otherwise, doubles for n < 2^50 (see
    below). Both use machine_vectors.h.
*/

#if NMOD_POLY_HAVE_KARATSUBA

#include "machine_vectors.h"

#define KARA_ALIGN64(p) ((void *) (((uintptr_t) (p) + 63) & ~(uintptr_t) 63))

/* x[0..3] + y[0..3] and x[0..3] - y[0..3] in one word */
FLINT_FORCE_INLINE vec4n
_sum4(const ulong * x, const ulong * y)
{
    return vec4n_add(vec4n_load_unaligned(x), vec4n_load_unaligned(y));
}

FLINT_FORCE_INLINE vec4n
_diff4(const ulong * x, const ulong * y)
{
    return vec4n_sub(vec4n_load_unaligned(x), vec4n_load_unaligned(y));
}

/* without this, gcc computes all the products of a leaf first and spills
   them, or keeps the loaded windows in stack temporaries */
#if defined(__GNUC__)
# define KARA_BARRIER4(a0, a1, a2, a3) \
    __asm__ volatile("" : "+x"(a0), "+x"(a1), "+x"(a2), "+x"(a3) : : "memory")
#else
# define KARA_BARRIER4(a0, a1, a2, a3)
#endif

/* blocks of 16 coefficients of an operand of a product, at a stride of 24,
   with 4 zeros before and after each block */
#define KARA_GS 24
#define KARA_GLEN(n) (KARA_GS * ((n) / 16) + 8)
#define KARA_PAD16(n) (((n) + 15) & ~WORD(15))

/* the number of leaves of the recursion for b blocks of 16 */
static slong
_leaves(slong b)
{
    return (b <= 1) ? 1 : 2 * _leaves((b + 1) / 2) + _leaves(b / 2);
}

/*
    A product cuts the longer operand into blocks of n = PAD16(len2)
    coefficients, and a middle product cuts its outputs (or the shorter
    operand) into blocks of n: a last block of r < n coefficients costs as
    much as a full one. Whether to compute it instead with blocks of
    PAD16(r) over the other dimension, of length l, which costs
    ceil(l / PAD16(r)) recursions of size PAD16(r) instead of one of size n
    (counted in leaves).
*/
static int
_split_last_block(slong r, slong n, slong l)
{
    slong rp = KARA_PAD16(r);

    if (r <= 0 || rp >= n)
        return 0;
    return ((l + rp - 1) / rp) * _leaves(rp / 16) < _leaves(n / 16);
}

/* ---- one word, products ------------------------------------------------ */

/* acc[t] = coefficients 4t to 4t + 3 of a[0..15] b[0..15], for t < 8; a is
   a block of a gapped operand */
#define KARA_MUL_T(t, j, acc) \
    if (4 * (t) - (j) >= -3 && 4 * (t) - (j) <= 15) \
        acc = vec4n_add(acc, \
                vec4n_mul(vec4n_load_unaligned(a + 4 * (t) - (j)), bj));

#define KARA_MUL_J(j) \
    { \
        const vec4n bj = vec4n_set_n(b[j]); \
        KARA_MUL_T(0, j, a0) KARA_MUL_T(1, j, a1) \
        KARA_MUL_T(2, j, a2) KARA_MUL_T(3, j, a3) \
        KARA_MUL_T(4, j, a4) KARA_MUL_T(5, j, a5) \
        KARA_MUL_T(6, j, a6) KARA_MUL_T(7, j, a7) \
        KARA_BARRIER4(a0, a1, a2, a3); \
        KARA_BARRIER4(a4, a5, a6, a7); \
    }

FLINT_FORCE_INLINE void
_leaf_mul_acc(vec4n * acc, const ulong * a, const ulong * b)
{
    vec4n a0, a1, a2, a3, a4, a5, a6, a7;

    a0 = a1 = a2 = a3 = a4 = a5 = a6 = a7 = vec4n_zero();

    KARA_MUL_J(0) KARA_MUL_J(1) KARA_MUL_J(2) KARA_MUL_J(3)
    KARA_MUL_J(4) KARA_MUL_J(5) KARA_MUL_J(6) KARA_MUL_J(7)
    KARA_MUL_J(8) KARA_MUL_J(9) KARA_MUL_J(10) KARA_MUL_J(11)
    KARA_MUL_J(12) KARA_MUL_J(13) KARA_MUL_J(14) KARA_MUL_J(15)

    acc[0] = a0; acc[1] = a1; acc[2] = a2; acc[3] = a3;
    acc[4] = a4; acc[5] = a5; acc[6] = a6; acc[7] = a7;
}

/* c[0..31] = a[0..15] b[0..15] (c[31] = 0), c 32-byte aligned */
FLINT_FORCE_INLINE void
_leaf_mul(ulong * c, const ulong * a, const ulong * b)
{
    vec4n acc[8];
    int t;

    _leaf_mul_acc(acc, a, b);
    for (t = 0; t < 8; t++)
        vec4n_store_aligned(c + 4 * t, acc[t]);
}

/* gapped block s = x + y (x, y blocks of gapped operands) */
FLINT_FORCE_INLINE void
_gap_add16(ulong * s, const ulong * x, const ulong * y)
{
    const vec4n z = vec4n_zero();

    vec4n_store_aligned(s - 4, z);
    vec4n_store_aligned(s, _sum4(x, y));
    vec4n_store_aligned(s + 4, _sum4(x + 4, y + 4));
    vec4n_store_aligned(s + 8, _sum4(x + 8, y + 8));
    vec4n_store_aligned(s + 12, _sum4(x + 12, y + 12));
    vec4n_store_aligned(s + 16, z);
}

FLINT_FORCE_INLINE void
_gap_copy16(ulong * s, const ulong * x)
{
    const vec4n z = vec4n_zero();

    vec4n_store_aligned(s - 4, z);
    vec4n_store_aligned(s, vec4n_load_aligned(x));
    vec4n_store_aligned(s + 4, vec4n_load_aligned(x + 4));
    vec4n_store_aligned(s + 8, vec4n_load_aligned(x + 8));
    vec4n_store_aligned(s + 12, vec4n_load_aligned(x + 12));
    vec4n_store_aligned(s + 16, z);
}

/* c[0..63] = a b for 32 coefficients, a and b gapped; the third product
   is combined in registers */
static void
_node32_mul(ulong * c, const ulong * a, const ulong * b)
{
    ulong buf[64 + 8];
    ulong * sa = (ulong *) KARA_ALIGN64(buf) + 4;
    ulong * sb = sa + KARA_GS;
    vec4n p[8];
    int j;

    /* sums first, so that their stores have retired when the leaf reads
       windows across them */
    _gap_add16(sa, a, a + KARA_GS);
    _gap_add16(sb, b, b + KARA_GS);
    _leaf_mul(c, a, b);
    _leaf_mul(c + 32, a + KARA_GS, b + KARA_GS);
    _leaf_mul_acc(p, sa, sb);

    for (j = 0; j < 4; j++)
    {
        vec4n u = vec4n_load_aligned(c + 4 * j);
        vec4n x = vec4n_load_aligned(c + 16 + 4 * j);
        vec4n y = vec4n_load_aligned(c + 32 + 4 * j);
        vec4n v = vec4n_load_aligned(c + 48 + 4 * j);

        vec4n_store_aligned(c + 16 + 4 * j,
                vec4n_sub(vec4n_add(x, p[j]), vec4n_add(u, y)));
        vec4n_store_aligned(c + 32 + 4 * j,
                vec4n_sub(vec4n_add(y, p[j + 4]), vec4n_add(x, v)));
    }
}

static slong
_mul_rec_scratch(slong n)
{
    slong s = 0, h;

    while (n > 32)
    {
        h = 16 * ((n / 16 + 1) / 2);
        s += 2 * KARA_GLEN(h) + 2 * h;
        n = h;
    }

    return s;
}

/* gapped copy of x (len coefficients, zero-padded to n, a multiple of 16)
   at g, 32-byte aligned; returns the address of block 0 */
static ulong *
_gap(ulong * g, const ulong * x, slong len, slong n)
{
    const vec4n z = vec4n_zero();
    ulong * g0 = g + 4, * blk;
    slong k, j;

    for (k = 0; k < n / 16; k++)
    {
        blk = g0 + KARA_GS * k;
        vec4n_store_aligned(blk - 4, z);
        if (16 * k + 16 <= len)
        {
            const ulong * xk = x + 16 * k;
            vec4n_store_aligned(blk, vec4n_load_unaligned(xk));
            vec4n_store_aligned(blk + 4, vec4n_load_unaligned(xk + 4));
            vec4n_store_aligned(blk + 8, vec4n_load_unaligned(xk + 8));
            vec4n_store_aligned(blk + 12, vec4n_load_unaligned(xk + 12));
        }
        else
        {
            for (j = 0; j < 16; j++)
                blk[j] = (16 * k + j < len) ? x[16 * k + j] : 0;
        }
        vec4n_store_aligned(blk + 16, z);
    }

    return g0;
}

/* c[0 .. 2n - 1] = a b (c[2n - 1] = 0), n a multiple of 16, a and b
   gapped, c and scr 32-byte aligned */
static void
_mul_rec(ulong * c, const ulong * a, const ulong * b, slong n, ulong * scr)
{
    slong m, hb, lb, h, l, i, k, e;
    ulong * sa, * sb, * p1, * rest;

    if (n == 16)
    {
        _leaf_mul(c, a, b);
        return;
    }

    if (n == 32)
    {
        _node32_mul(c, a, b);
        return;
    }

    m = n / 16;
    hb = (m + 1) / 2;
    lb = m - hb;
    h = 16 * hb;
    l = 16 * lb;
    sa = scr + 4;
    sb = scr + KARA_GLEN(h) + 4;
    p1 = scr + 2 * KARA_GLEN(h);
    rest = p1 + 2 * h;

    for (k = 0; k < hb; k++)
    {
        if (k < lb)
        {
            _gap_add16(sa + KARA_GS * k, a + KARA_GS * k,
                       a + KARA_GS * (k + hb));
            _gap_add16(sb + KARA_GS * k, b + KARA_GS * k,
                       b + KARA_GS * (k + hb));
        }
        else
        {
            _gap_copy16(sa + KARA_GS * k, a + KARA_GS * k);
            _gap_copy16(sb + KARA_GS * k, b + KARA_GS * k);
        }
    }

    _mul_rec(c, a, b, h, rest);
    _mul_rec(c + 2 * h, a + KARA_GS * hb, b + KARA_GS * hb, l, rest);
    _mul_rec(p1, sa, sb, h, rest);

    /* c[h + i] += p1[i] - P0[i] - P2[i], both halves in one pass (P2 has
       2l coefficients: its high half stops at e = 2l - h) */
    e = 2 * l - h;
    for (i = 0; i < h; i += 4)
    {
        vec4n u = vec4n_load_aligned(c + i);
        vec4n x = vec4n_load_aligned(c + h + i);
        vec4n y = vec4n_load_aligned(c + 2 * h + i);
        vec4n v = (i < e) ? vec4n_load_aligned(c + 3 * h + i) : vec4n_zero();
        vec4n q0 = vec4n_load_aligned(p1 + i);
        vec4n q1 = vec4n_load_aligned(p1 + h + i);

        vec4n_store_aligned(c + h + i,
                            vec4n_sub(vec4n_add(x, q0), vec4n_add(u, y)));
        vec4n_store_aligned(c + 2 * h + i,
                            vec4n_sub(vec4n_add(y, q1), vec4n_add(x, v)));
    }
}

/* out[0 .. len1 + len2 - 2] = a b, exact one-word values, for len1 >= len2
   >= 1; out 32-byte aligned, with room for len1 + 2 PAD16(len2) words */
static void
_mul_raw(ulong * out, const ulong * a, slong len1, const ulong * b, slong len2)
{
    slong n = KARA_PAD16(len2), k0, i;
    ulong * w, * ga, * gb, * p, * scr, * pa, * pb;
    TMP_INIT;

    TMP_START;
    w = KARA_ALIGN64(TMP_ALLOC((2 * KARA_GLEN(n) + 2 * n
                                + _mul_rec_scratch(n) + 8) * sizeof(ulong)));
    ga = w;
    gb = ga + KARA_GLEN(n);
    p = gb + KARA_GLEN(n);
    scr = p + 2 * n;

    pb = _gap(gb, b, len2, n);
    _mul_rec(out, _gap(ga, a, FLINT_MIN(n, len1), n), pb, n, scr);

    for (k0 = n; k0 < len1; k0 += n)
    {
        pa = _gap(ga, a + k0, FLINT_MIN(n, len1 - k0), n);
        _mul_rec(p, pa, pb, n, scr);

        /* the first n overlap the previous block, whose last word is 0;
           the whole block is written, its last word 0 too */
        for (i = 0; i < n; i += 4)
            vec4n_store_aligned(out + k0 + i, _sum4(out + k0 + i, p + i));
        for ( ; i < 2 * n; i += 4)
            vec4n_store_aligned(out + k0 + i, vec4n_load_aligned(p + i));
    }

    TMP_END;
}

/* ---- one word, middle products (transposed Karatsuba) ------------------ */

/* MP(a, b) for n coefficients: r[k] = sum_{i < n} a[k + n - 1 - i] b[i] for
   k < n, a having 2n - 1 coefficients */

/* acc[t] += coefficients 4t to 4t + 3 of MP(a, b) for 16 coefficients,
   t < 4; the operands are signed, below 2^31 in absolute value, stored in
   words in two's complement (vec4n_mul_signed reads their low 32 bits) */
#define KARA_MID_TERM(k) vec4n_mul_signed(vec4n_load_unaligned(a + (k)), bi)

#define KARA_MID_I(i) \
    { \
        const vec4n bi = vec4n_set_n(b[i]); \
        a0 = vec4n_add(a0, KARA_MID_TERM(15 - (i))); \
        a1 = vec4n_add(a1, KARA_MID_TERM(19 - (i))); \
        a2 = vec4n_add(a2, KARA_MID_TERM(23 - (i))); \
        a3 = vec4n_add(a3, KARA_MID_TERM(27 - (i))); \
        KARA_BARRIER4(a0, a1, a2, a3); \
    }

FLINT_FORCE_INLINE void
_leaf_mid_acc(vec4n * acc, const ulong * a, const ulong * b)
{
    vec4n a0 = acc[0], a1 = acc[1], a2 = acc[2], a3 = acc[3];

    KARA_MID_I(0) KARA_MID_I(1) KARA_MID_I(2) KARA_MID_I(3)
    KARA_MID_I(4) KARA_MID_I(5) KARA_MID_I(6) KARA_MID_I(7)
    KARA_MID_I(8) KARA_MID_I(9) KARA_MID_I(10) KARA_MID_I(11)
    KARA_MID_I(12) KARA_MID_I(13) KARA_MID_I(14) KARA_MID_I(15)

    acc[0] = a0; acc[1] = a1; acc[2] = a2; acc[3] = a3;
}

FLINT_FORCE_INLINE void
_leaf_mid(ulong * r, const ulong * a, const ulong * b)
{
    vec4n acc[4];

    acc[0] = acc[1] = acc[2] = acc[3] = vec4n_zero();
    _leaf_mid_acc(acc, a, b);
    vec4n_store_unaligned(r, acc[0]);
    vec4n_store_unaligned(r + 4, acc[1]);
    vec4n_store_unaligned(r + 8, acc[2]);
    vec4n_store_unaligned(r + 12, acc[3]);
}

/* s[0..30] = x[0..30] + y[0..30], s 32-byte aligned */
FLINT_FORCE_INLINE void
_add31(ulong * s, const ulong * x, const ulong * y)
{
    int i;

    for (i = 0; i < 28; i += 4)
        vec4n_store_aligned(s + i, _sum4(x + i, y + i));
    s[28] = x[28] + y[28];
    s[29] = x[29] + y[29];
    s[30] = x[30] + y[30];
}

/* MP for 32 coefficients, the outer leaves starting from +- beta */
static void
_node32_mid(ulong * r, const ulong * a, const ulong * b)
{
    ulong buf[96 + 8];
    ulong * db = (ulong *) KARA_ALIGN64(buf), * s1 = db + 32, * s2 = db + 64;
    vec4n beta[4], acc[4];
    int t;

    for (t = 0; t < 16; t += 4)
        vec4n_store_aligned(db + t, _diff4(b + t, b + 16 + t));
    _add31(s1, a, a + 16);
    _add31(s2, a + 16, a + 32);

    beta[0] = beta[1] = beta[2] = beta[3] = vec4n_zero();
    _leaf_mid_acc(beta, a + 16, db);

    for (t = 0; t < 4; t++)
        acc[t] = beta[t];
    _leaf_mid_acc(acc, s1, b + 16);
    for (t = 0; t < 4; t++)
        vec4n_store_unaligned(r + 4 * t, acc[t]);

    for (t = 0; t < 4; t++)
        acc[t] = vec4n_sub(vec4n_zero(), beta[t]);
    _leaf_mid_acc(acc, s2, b);
    for (t = 0; t < 4; t++)
        vec4n_store_unaligned(r + 16 + 4 * t, acc[t]);
}

static slong
_mid_rec_scratch(slong n)
{
    slong s = 0, h;

    while (n > 32)
    {
        h = 16 * ((n / 16 + 1) / 2);
        s += 6 * h;
        n = h;
    }

    return s;
}

/*
    r[0 .. n - 1] = MP(a, b), n a multiple of 16, scr 32-byte aligned.

    With b = b0 + x^h b1, b0 of h and b1 of l coefficients, and the windows
    Y0 = a[l - h ..], Y1 = a[l ..], Y2 = a[n ..] of 2h - 1 coefficients (zero
    outside a), r[0 .. h - 1] = MP(Y1, b0) + MP(Y0, b1') and r[h .. n - 1] =
    MP(Y2, b0) + MP(Y1, b1') up to l outputs, b1' being b1 padded with
    h - l zeros; so with beta = MP(Y1, b0 - b1'),

        r[0 .. h - 1] = MP(Y0 + Y1, b1') + beta,
        r[h .. n - 1] = MP(Y1 + Y2, b0) - beta    (first l outputs).

    For h = l (an even number of blocks), these are three middle products
    of size h. For h = l + 16 (an odd number), MP(Y0 + Y1, b1') is a middle
    product of size l for its first l outputs (on (Y0 + Y1)[16 ..], which
    does not read Y0 outside a) and a row of l / 16 leaves for the last 16;
    the first l outputs of MP(Y1 + Y2, b0) are a middle product of size l
    with the first l coefficients of b0 and a row of l / 16 leaves with its
    last 16.
*/
static void
_mid_rec(ulong * r, const ulong * a, const ulong * b, slong n, ulong * scr)
{
    slong m, h, l, i, q;

    if (n == 16)
    {
        _leaf_mid(r, a, b);
        return;
    }

    if (n == 32)
    {
        _node32_mid(r, a, b);
        return;
    }

    m = n / 16;
    h = 16 * ((m + 1) / 2);
    l = n - h;

    if (h == l)
    {
        ulong * s1 = scr, * s2 = scr + 2 * h, * db = scr + 4 * h;
        ulong * beta = db + h, * rest = beta + h;
        const ulong * A0 = a, * A1 = a + h, * A2 = a + 2 * h;
        const ulong * b0 = b, * b1 = b + h;

        for (i = 0; i < h; i += 4)
            vec4n_store_aligned(db + i, _diff4(b0 + i, b1 + i));
        for (i = 0; i < 2 * h - 4; i += 4)
        {
            vec4n_store_aligned(s1 + i, _sum4(A0 + i, A1 + i));
            vec4n_store_aligned(s2 + i, _sum4(A1 + i, A2 + i));
        }
        for ( ; i < 2 * h - 1; i++)
        {
            s1[i] = A0[i] + A1[i];
            s2[i] = A1[i] + A2[i];
        }

        _mid_rec(beta, A1, db, h, rest);
        _mid_rec(r, s1, b1, h, rest);
        _mid_rec(r + h, s2, b0, h, rest);

        for (i = 0; i < h; i += 4)
        {
            vec4n bt = vec4n_load_aligned(beta + i);
            vec4n_store_unaligned(r + i,
                                  vec4n_add(vec4n_load_unaligned(r + i), bt));
            vec4n_store_unaligned(r + h + i,
                    vec4n_sub(vec4n_load_unaligned(r + h + i), bt));
        }
    }
    else
    {
        /* t1 = (Y0 + Y1)[16 ..]: t1[v] = a[v] + a[v + h] for v < 2l + 15;
           s2 = Y1 + Y2: s2[u] = a[l + u] + a[n + u] for u < n - 1 */
        ulong * db = scr, * beta = db + h, * t1 = beta + h, * s2 = t1 + 2 * h;
        ulong * rest = s2 + 2 * h;
        const ulong * b1 = b + h;
        vec4n acc[4];
        int t;

        for (i = 0; i < l; i += 4)
            vec4n_store_aligned(db + i, _diff4(b + i, b1 + i));
        for ( ; i < h; i += 4)
            vec4n_store_aligned(db + i, vec4n_load_unaligned(b + i));
        for (i = 0; i + 4 <= 2 * l + 15; i += 4)
            vec4n_store_aligned(t1 + i, _sum4(a + i, a + h + i));
        for ( ; i < 2 * l + 15; i++)
            t1[i] = a[i] + a[h + i];
        for (i = 0; i + 4 <= n - 1; i += 4)
            vec4n_store_aligned(s2 + i, _sum4(a + l + i, a + n + i));
        for ( ; i < n - 1; i++)
            s2[i] = a[l + i] + a[n + i];

        _mid_rec(beta, a + l, db, h, rest);

        /* r[0 .. l - 1] and r[l .. h - 1] */
        _mid_rec(r, t1, b1, l, rest);
        for (t = 0; t < 4; t++)
            acc[t] = vec4n_load_aligned(beta + l + 4 * t);
        for (q = 0; q < l / 16; q++)
            _leaf_mid_acc(acc, t1 + n - 32 - 16 * q, b1 + 16 * q);
        for (t = 0; t < 4; t++)
            vec4n_store_unaligned(r + l + 4 * t, acc[t]);
        for (i = 0; i < l; i += 4)
            vec4n_store_unaligned(r + i, _sum4(r + i, beta + i));

        /* r[h .. n - 1] */
        _mid_rec(r + h, s2 + 16, b, l, rest);
        for (q = 0; q < l / 16; q++)
        {
            for (t = 0; t < 4; t++)
                acc[t] = _diff4(r + h + 16 * q + 4 * t, beta + 16 * q + 4 * t);
            _leaf_mid_acc(acc, s2 + 16 * q, b + l);
            for (t = 0; t < 4; t++)
                vec4n_store_unaligned(r + h + 16 * q + 4 * t, acc[t]);
        }
    }
}

/* x[0 .. count - 1] = a[start .. start + count - 1], zero outside [0, len) */
static void
_gather(ulong * x, const ulong * a, slong len, slong start, slong count)
{
    slong lo = FLINT_MAX(start, 0), hi = FLINT_MIN(start + count, len);

    if (lo >= hi)
    {
        flint_mpn_zero(x, count);
        return;
    }

    if (lo > start)
        flint_mpn_zero(x, lo - start);
    flint_mpn_copyi(x + (lo - start), a + lo, hi - lo);
    if (hi - start < count)
        flint_mpn_zero(x + (hi - start), count - (hi - start));
}

/* out[k] = sum_{i < len2} a[nlo + k - i] b[i] for k < m (a zero outside
   [0, len1)), exact one-word values; out is not aligned */
static void
_mid_raw(ulong * out, const ulong * a, slong len1, const ulong * b, slong len2,
         slong nlo, slong m)
{
    slong N, k0, s, i, st;
    const ulong * xp, * yp;
    ulong * w, * x, * y, * r, * scr;
    TMP_INIT;

    TMP_START;

    N = KARA_PAD16(FLINT_MIN(len2, m));
    w = KARA_ALIGN64(TMP_ALLOC((4 * N + _mid_rec_scratch(N) + 8)
                               * sizeof(ulong)));
    x = w;
    y = x + 2 * N;
    r = y + N;
    scr = r + N;

    if (m >= len2)
    {
        /* blocks of N outputs */
        if (len2 == N)
            yp = b;
        else
        {
            _gather(y, b, len2, 0, N);
            yp = y;
        }

        for (k0 = 0; k0 < m; k0 += N)
        {
            st = nlo + k0 - N + 1;
            if (st >= 0 && st + 2 * N - 1 <= len1)
                xp = a + st;
            else
            {
                _gather(x, a, len1, st, 2 * N - 1);
                xp = x;
            }

            if (m - k0 >= N)
                _mid_rec(out + k0, xp, yp, N, scr);
            else if (k0 > 0 && _split_last_block(m - k0, N, len2))
            {
                _mid_raw(out + k0, a, len1, b, len2, nlo + k0, m - k0);
                break;
            }
            else
            {
                _mid_rec(r, xp, yp, N, scr);
                flint_mpn_copyi(out + k0, r, m - k0);
            }
        }
    }
    else
    {
        /* blocks of N coefficients of b */
        for (s = 0; s * N < len2; s++)
        {
            if ((s + 1) * N <= len2)
                yp = b + s * N;
            else if (s > 0 && _split_last_block(len2 - s * N, N, m))
            {
                ulong * t = TMP_ALLOC(m * sizeof(ulong));
                _mid_raw(t, a, len1, b + s * N, len2 - s * N, nlo - s * N, m);
                for (i = 0; i < m; i++)
                    out[i] += t[i];
                break;
            }
            else
            {
                _gather(y, b, len2, s * N, N);
                yp = y;
            }

            st = nlo - s * N - N + 1;
            if (st >= 0 && st + 2 * N - 1 <= len1)
                xp = a + st;
            else
            {
                _gather(x, a, len1, st, 2 * N - 1);
                xp = x;
            }

            _mid_rec(r, xp, yp, N, scr);
            if (s == 0)
                flint_mpn_copyi(out, r, m);
            else
                for (i = 0; i < m; i++)
                    out[i] += r[i];
        }
    }

    TMP_END;
}

/* ---- doubles (n < 2^50) ------------------------------------------------ */

#if NMOD_POLY_HAVE_KARATSUBA_FP50

/*
    Doubles, with the arithmetic of machine_vectors.h. The inputs are
    residues in [0, n); the sums made by the recursion are reduced as
    v - n round(v/n) (fp50_red: at most n/2 + 2^-52 |v| in absolute value
    for |v| < 2^53, so within [-n/2 - 2, n/2 + 2]).

    A product w b is split exactly as h + l, with h = fl(w b) and
    l = fma(w, b, -h), and h is reduced at once as e = h - q n with
    q = round(w (b/n)): b/n is computed once per coefficient of b, so that
    q does not wait for h. For n < 2^50, |w b/n - q| <= 3/4 (q is rounded
    with the constant 3 2^51, in one fused operation) and
    |l| < n^2 2^-53 < n/8, so |e| < 7n/8.

    A leaf sums the e and the l of each output separately and reduces the
    e after their first 8 terms, so |E| < n/2 + 2 + 7n. Since the l of the
    last 8 terms cancel in E + L, |E + L| < n/2 + 2 + 8 (3n/4) + 8 (n/8),
    below 2^53; E + L is reduced once more at the end of the leaf.
*/
typedef struct
{
    vec4d P, PINV, C;
    double p, pinv;
}
fp50_ctx;

static void
fp50_init(fp50_ctx * c, ulong n)
{
    c->p = (double) n;
    c->pinv = 1.0 / c->p;
    c->P = vec4d_set_d(c->p);
    c->PINV = vec4d_set_d(c->pinv);
    c->C = vec4d_set_d(6755399441055744.0);     /* 3 2^51 */
}

/* v - n round(v/n), as vec4d_reduce_to_pm1n, but always rounding with the
   constant 3 2^51, as in the leaves (vec4d_reduce_to_pm1n uses vroundpd
   where FLINT_FFT_SMALL_ROUND_USE_NATIVE is set): 2 uops on every target,
   against 3 for mul and vroundpd on Intel; |v/n| <= 2^51 holds here */
FLINT_FORCE_INLINE vec4d
fp50_red(vec4d v, const fp50_ctx * c)
{
    vec4d q = vec4d_sub(vec4d_fmadd_fused(v, c->PINV, c->C), c->C);

    return vec4d_fnmadd_fused(q, c->P, v);
}

FLINT_FORCE_INLINE double
fp50_red1(double v, const fp50_ctx * c)
{
    const double C = 6755399441055744.0;    /* 3 2^51 */
    double q = vec1d_sub(vec1d_fmadd(v, c->pinv, C), C);

    return vec1d_fnmadd(q, c->p, v);
}

/* x[0..3] + y[0..3] and x[0..3] - y[0..3], reduced */
FLINT_FORCE_INLINE vec4d
fp50_sum4(const double * x, const double * y, const fp50_ctx * c)
{
    return fp50_red(vec4d_add(vec4d_load_unaligned(x),
                              vec4d_load_unaligned(y)), c);
}

FLINT_FORCE_INLINE vec4d
fp50_diff4(const double * x, const double * y, const fp50_ctx * c)
{
    return fp50_red(vec4d_sub(vec4d_load_unaligned(x),
                              vec4d_load_unaligned(y)), c);
}

#if defined(__GNUC__)
# define KARA_BARRIER8D(a0, a1, a2, a3, a4, a5, a6, a7) \
    __asm__ volatile("" : "+x"(a0), "+x"(a1), "+x"(a2), "+x"(a3), \
                          "+x"(a4), "+x"(a5), "+x"(a6), "+x"(a7) : : "memory")
#else
# define KARA_BARRIER8D(a0, a1, a2, a3, a4, a5, a6, a7)
#endif

/* one product w b into E (reduced part) and L (low part), bq = b/n; the
   three fused operations must be fused for the bounds above */
#define FP50_STEP(w, bj, bq, E, L) \
    do { \
        vec4d _h = vec4d_mul(w, bj); \
        vec4d _q = vec4d_sub(vec4d_fmadd_fused(w, bq, cx->C), cx->C); \
        vec4d _l = vec4d_fmsub_fused(w, bj, _h); \
        E = vec4d_add(E, vec4d_fnmadd_fused(_q, cx->P, _h)); \
        L = vec4d_add(L, _l); \
    } while (0)

/* product leaf, outputs 16 half to 16 half + 15 (half = 0, 1) into
   acc[0..3], reduced; a gapped */
#define FP50_MUL_T(t, j, E, L) \
    if (4 * (t) - (j) >= -3 && 4 * (t) - (j) <= 15) \
    { \
        vec4d _w = vec4d_load_unaligned(a + 4 * (t) - (j)); \
        FP50_STEP(_w, bj, bq, E, L); \
    }

#define FP50_MUL_J(j) \
    if ((j) < 8 || (j) < jn) \
    { \
        const vec4d bj = vec4d_set_d(b[j]); \
        const vec4d bq = vec4d_mul(bj, cx->PINV); \
        if (half == 0) \
        { \
            FP50_MUL_T(0, j, E0, L0) FP50_MUL_T(1, j, E1, L1) \
            FP50_MUL_T(2, j, E2, L2) FP50_MUL_T(3, j, E3, L3) \
        } \
        else \
        { \
            FP50_MUL_T(4, j, E0, L0) FP50_MUL_T(5, j, E1, L1) \
            FP50_MUL_T(6, j, E2, L2) FP50_MUL_T(7, j, E3, L3) \
        } \
        if ((j) == 7) \
        { \
            E0 = fp50_red(E0, cx); E1 = fp50_red(E1, cx); \
            E2 = fp50_red(E2, cx); E3 = fp50_red(E3, cx); \
        } \
        KARA_BARRIER8D(E0, E1, E2, E3, L0, L1, L2, L3); \
    }

/* the terms of b[0 .. jn - 1] only (jn >= 8) */
FLINT_FORCE_INLINE void
fp50_leaf_mul_half(vec4d * acc, const double * a, const double * b, int half,
                   slong jn, const fp50_ctx * cx)
{
    const vec4d z = vec4d_zero();
    vec4d E0 = z, E1 = z, E2 = z, E3 = z, L0 = z, L1 = z, L2 = z, L3 = z;

    FP50_MUL_J(0) FP50_MUL_J(1) FP50_MUL_J(2) FP50_MUL_J(3)
    FP50_MUL_J(4) FP50_MUL_J(5) FP50_MUL_J(6) FP50_MUL_J(7)
    FP50_MUL_J(8) FP50_MUL_J(9) FP50_MUL_J(10) FP50_MUL_J(11)
    FP50_MUL_J(12) FP50_MUL_J(13) FP50_MUL_J(14) FP50_MUL_J(15)

    acc[0] = fp50_red(vec4d_add(E0, L0), cx);
    acc[1] = fp50_red(vec4d_add(E1, L1), cx);
    acc[2] = fp50_red(vec4d_add(E2, L2), cx);
    acc[3] = fp50_red(vec4d_add(E3, L3), cx);
}

FLINT_FORCE_INLINE void
fp50_leaf_mul(double * c, const double * a, const double * b,
              const fp50_ctx * cx)
{
    vec4d acc[4];
    int t;

    fp50_leaf_mul_half(acc, a, b, 0, 16, cx);
    for (t = 0; t < 4; t++)
        vec4d_store_aligned(c + 4 * t, acc[t]);
    fp50_leaf_mul_half(acc, a, b, 1, 16, cx);
    for (t = 0; t < 4; t++)
        vec4d_store_aligned(c + 16 + 4 * t, acc[t]);
}

/* the leaf for b of jn <= 16 coefficients (shorter operands of at most 16
   coefficients): the rows from max(jn, 8) on are skipped, and outputs 16 to
   31 are zero unless high */
static void
fp50_leaf_mul_short(double * c, const double * a, const double * b, slong jn,
                    int high, const fp50_ctx * cx)
{
    vec4d acc[4];
    int t;

    fp50_leaf_mul_half(acc, a, b, 0, jn, cx);
    for (t = 0; t < 4; t++)
        vec4d_store_aligned(c + 4 * t, acc[t]);
    if (high)
        fp50_leaf_mul_half(acc, a, b, 1, jn, cx);
    else
        for (t = 0; t < 4; t++)
            acc[t] = vec4d_zero();
    for (t = 0; t < 4; t++)
        vec4d_store_aligned(c + 16 + 4 * t, acc[t]);
}

FLINT_FORCE_INLINE void
fp50_gap_add16(double * s, const double * x, const double * y,
               const fp50_ctx * cx)
{
    const vec4d z = vec4d_zero();
    int i;

    vec4d_store_aligned(s - 4, z);
    for (i = 0; i < 16; i += 4)
        vec4d_store_aligned(s + i, fp50_sum4(x + i, y + i, cx));
    vec4d_store_aligned(s + 16, z);
}

FLINT_FORCE_INLINE void
fp50_gap_copy16(double * s, const double * x)
{
    const vec4d z = vec4d_zero();
    int i;

    vec4d_store_aligned(s - 4, z);
    for (i = 0; i < 16; i += 4)
        vec4d_store_aligned(s + i, vec4d_load_aligned(x + i));
    vec4d_store_aligned(s + 16, z);
}

static slong
fp50_mul_rec_scratch(slong n)
{
    slong s = 0, h;

    while (n > 16)
    {
        h = 16 * ((n / 16 + 1) / 2);
        s += 2 * KARA_GLEN(h) + 2 * h;
        n = h;
    }

    return s;
}

/* c[0 .. 2n - 1] = a b, reduced; as _mul_rec */
static void
fp50_mul_rec(double * c, const double * a, const double * b, slong n,
             double * scr, const fp50_ctx * cx)
{
    slong m, hb, lb, h, l, i, k, e;
    double * sa, * sb, * p1, * rest;

    if (n == 16)
    {
        fp50_leaf_mul(c, a, b, cx);
        return;
    }

    m = n / 16;
    hb = (m + 1) / 2;
    lb = m - hb;
    h = 16 * hb;
    l = 16 * lb;
    sa = scr + 4;
    sb = scr + KARA_GLEN(h) + 4;
    p1 = scr + 2 * KARA_GLEN(h);
    rest = p1 + 2 * h;

    for (k = 0; k < hb; k++)
    {
        if (k < lb)
        {
            fp50_gap_add16(sa + KARA_GS * k, a + KARA_GS * k,
                           a + KARA_GS * (k + hb), cx);
            fp50_gap_add16(sb + KARA_GS * k, b + KARA_GS * k,
                           b + KARA_GS * (k + hb), cx);
        }
        else
        {
            fp50_gap_copy16(sa + KARA_GS * k, a + KARA_GS * k);
            fp50_gap_copy16(sb + KARA_GS * k, b + KARA_GS * k);
        }
    }

    fp50_mul_rec(c, a, b, h, rest, cx);
    fp50_mul_rec(c + 2 * h, a + KARA_GS * hb, b + KARA_GS * hb, l, rest, cx);
    fp50_mul_rec(p1, sa, sb, h, rest, cx);

    e = 2 * l - h;
    for (i = 0; i < h; i += 4)
    {
        vec4d u = vec4d_load_aligned(c + i);
        vec4d x = vec4d_load_aligned(c + h + i);
        vec4d y = vec4d_load_aligned(c + 2 * h + i);
        vec4d v = (i < e) ? vec4d_load_aligned(c + 3 * h + i) : vec4d_zero();
        vec4d q0 = vec4d_load_aligned(p1 + i);
        vec4d q1 = vec4d_load_aligned(p1 + h + i);

        vec4d_store_aligned(c + h + i, fp50_red(vec4d_sub(vec4d_add(x, q0),
                                                vec4d_add(u, y)), cx));
        vec4d_store_aligned(c + 2 * h + i, fp50_red(vec4d_sub(vec4d_add(y, q1),
                                                    vec4d_add(x, v)), cx));
    }
}

/* middle-product leaf: acc[t] += 4 outputs of MP(a, b) for 16 coefficients,
   then reduced; acc enters reduced */
#define FP50_MID_I(i) \
    { \
        const vec4d bi = vec4d_set_d(b[i]); \
        const vec4d bq = vec4d_mul(bi, cx->PINV); \
        vec4d w; \
        w = vec4d_load_unaligned(a + 15 - (i)); FP50_STEP(w, bi, bq, E0, L0); \
        w = vec4d_load_unaligned(a + 19 - (i)); FP50_STEP(w, bi, bq, E1, L1); \
        w = vec4d_load_unaligned(a + 23 - (i)); FP50_STEP(w, bi, bq, E2, L2); \
        w = vec4d_load_unaligned(a + 27 - (i)); FP50_STEP(w, bi, bq, E3, L3); \
        if ((i) == 7) \
        { \
            E0 = fp50_red(E0, cx); E1 = fp50_red(E1, cx); \
            E2 = fp50_red(E2, cx); E3 = fp50_red(E3, cx); \
        } \
        KARA_BARRIER8D(E0, E1, E2, E3, L0, L1, L2, L3); \
    }

FLINT_FORCE_INLINE void
fp50_leaf_mid_acc(vec4d * acc, const double * a, const double * b,
                  const fp50_ctx * cx)
{
    const vec4d z = vec4d_zero();
    vec4d E0 = acc[0], E1 = acc[1], E2 = acc[2], E3 = acc[3];
    vec4d L0 = z, L1 = z, L2 = z, L3 = z;

    FP50_MID_I(0) FP50_MID_I(1) FP50_MID_I(2) FP50_MID_I(3)
    FP50_MID_I(4) FP50_MID_I(5) FP50_MID_I(6) FP50_MID_I(7)
    FP50_MID_I(8) FP50_MID_I(9) FP50_MID_I(10) FP50_MID_I(11)
    FP50_MID_I(12) FP50_MID_I(13) FP50_MID_I(14) FP50_MID_I(15)

    acc[0] = fp50_red(vec4d_add(E0, L0), cx);
    acc[1] = fp50_red(vec4d_add(E1, L1), cx);
    acc[2] = fp50_red(vec4d_add(E2, L2), cx);
    acc[3] = fp50_red(vec4d_add(E3, L3), cx);
}

/* d = x + y reduced, len doubles, d aligned */
static void
fp50_add(double * d, const double * x, const double * y, slong len,
         const fp50_ctx * cx)
{
    slong i;

    for (i = 0; i + 4 <= len; i += 4)
        vec4d_store_aligned(d + i, fp50_sum4(x + i, y + i, cx));
    for ( ; i < len; i++)
        d[i] = fp50_red1(x[i] + y[i], cx);
}

/* d = x - y reduced, len a multiple of 4 (blocks of b), d aligned */
static void
fp50_sub(double * d, const double * x, const double * y, slong len,
         const fp50_ctx * cx)
{
    slong i;

    FLINT_ASSERT(len % 4 == 0);

    for (i = 0; i < len; i += 4)
        vec4d_store_aligned(d + i, fp50_diff4(x + i, y + i, cx));
}

static slong
fp50_mid_rec_scratch(slong n)
{
    slong s = 0, h;

    while (n > 16)
    {
        h = 16 * ((n / 16 + 1) / 2);
        s += 6 * h;
        n = h;
    }

    return s;
}

/* r[0 .. n - 1] = MP(a, b) reduced, as _mid_rec (same splits) */
static void
fp50_mid_rec(double * r, const double * a, const double * b, slong n,
             double * scr, const fp50_ctx * cx)
{
    slong m, h, l, i, q;
    vec4d acc[4];
    int t;

    if (n == 16)
    {
        acc[0] = acc[1] = acc[2] = acc[3] = vec4d_zero();
        fp50_leaf_mid_acc(acc, a, b, cx);
        for (t = 0; t < 4; t++)
            vec4d_store_unaligned(r + 4 * t, acc[t]);
        return;
    }

    m = n / 16;
    h = 16 * ((m + 1) / 2);
    l = n - h;

    if (h == l)
    {
        double * s1 = scr, * s2 = scr + 2 * h, * db = scr + 4 * h;
        double * beta = db + h, * rest = beta + h;

        fp50_sub(db, b, b + h, h, cx);
        fp50_add(s1, a, a + h, 2 * h - 1, cx);
        fp50_add(s2, a + h, a + 2 * h, 2 * h - 1, cx);

        if (h == 16)
        {
            /* node of 32: the outer leaves start from +- beta */
            vec4d bt[4];

            bt[0] = bt[1] = bt[2] = bt[3] = vec4d_zero();
            fp50_leaf_mid_acc(bt, a + h, db, cx);
            for (t = 0; t < 4; t++)
                acc[t] = bt[t];
            fp50_leaf_mid_acc(acc, s1, b + h, cx);
            for (t = 0; t < 4; t++)
                vec4d_store_unaligned(r + 4 * t, acc[t]);
            for (t = 0; t < 4; t++)
                acc[t] = vec4d_sub(vec4d_zero(), bt[t]);
            fp50_leaf_mid_acc(acc, s2, b, cx);
            for (t = 0; t < 4; t++)
                vec4d_store_unaligned(r + 16 + 4 * t, acc[t]);
            return;
        }

        fp50_mid_rec(beta, a + h, db, h, rest, cx);
        fp50_mid_rec(r, s1, b + h, h, rest, cx);
        fp50_mid_rec(r + h, s2, b, h, rest, cx);

        for (i = 0; i < h; i += 4)
        {
            vec4d_store_unaligned(r + i, fp50_sum4(r + i, beta + i, cx));
            vec4d_store_unaligned(r + h + i,
                                  fp50_diff4(r + h + i, beta + i, cx));
        }
    }
    else
    {
        double * db = scr, * beta = db + h, * t1 = beta + h, * s2 = t1 + 2 * h;
        double * rest = s2 + 2 * h;
        const double * b1 = b + h;

        fp50_sub(db, b, b1, l, cx);
        for (i = l; i < h; i += 4)
            vec4d_store_aligned(db + i, vec4d_load_unaligned(b + i));
        fp50_add(t1, a, a + h, 2 * l + 15, cx);
        fp50_add(s2, a + l, a + n, n - 1, cx);

        fp50_mid_rec(beta, a + l, db, h, rest, cx);

        fp50_mid_rec(r, t1, b1, l, rest, cx);
        for (t = 0; t < 4; t++)
            acc[t] = vec4d_load_aligned(beta + l + 4 * t);
        for (q = 0; q < l / 16; q++)
            fp50_leaf_mid_acc(acc, t1 + n - 32 - 16 * q, b1 + 16 * q, cx);
        for (t = 0; t < 4; t++)
            vec4d_store_unaligned(r + l + 4 * t, acc[t]);
        for (i = 0; i < l; i += 4)
            vec4d_store_unaligned(r + i, fp50_sum4(r + i, beta + i, cx));

        fp50_mid_rec(r + h, s2 + 16, b, l, rest, cx);
        for (q = 0; q < l / 16; q++)
        {
            for (t = 0; t < 4; t++)
                acc[t] = fp50_diff4(r + h + 16 * q + 4 * t,
                                    beta + 16 * q + 4 * t, cx);
            fp50_leaf_mid_acc(acc, s2 + 16 * q, b + l, cx);
            for (t = 0; t < 4; t++)
                vec4d_store_unaligned(r + h + 16 * q + 4 * t, acc[t]);
        }
    }
}

/* words below 2^52 to doubles, and reduced doubles back to residues (the
   correction of the sign tests v < 0, so that -0.0 would give 0, not n) */
static void
fp50_to_double(double * d, nn_srcptr x, slong len)
{
    slong i;

    for (i = 0; i + 4 <= len; i += 4)
        vec4d_store_unaligned(d + i,
                vec4n_convert_limited_vec4d(vec4n_load_unaligned(x + i)));
    for ( ; i < len; i++)
        d[i] = (double) x[i];
}

static void
fp50_from_double(nn_ptr x, const double * d, slong len, const fp50_ctx * cx)
{
    const vec4d z = vec4d_zero();
    slong i;

    for (i = 0; i + 4 <= len; i += 4)
    {
        vec4d v = fp50_red(vec4d_load_unaligned(d + i), cx);
        v = vec4d_add(v, vec4d_bit_and(vec4d_cmp_lt(v, z), cx->P));
        vec4n_store_unaligned(x + i, vec4d_convert_limited_vec4n(v));
    }
    for ( ; i < len; i++)
    {
        double v = fp50_red1(d[i], cx);
        if (v < 0)
            v += cx->p;
        x[i] = (ulong) v;
    }
}

/* res[0 .. nout - 1] = the first nout coefficients of a b mod n, for
   len1 >= len2 >= 1 and nout <= len1 + len2 - 1 */
static void
fp50_mul(nn_ptr res, nn_srcptr a, slong len1, nn_srcptr b, slong len2,
         slong nout, const nmod_t * mod)
{
    fp50_ctx cx;
    slong n = KARA_PAD16(len2), k0, i;
    double * w, * ad, * bd, * ga, * gb, * p, * out, * scr, * pa, * pb;
    TMP_INIT;

    fp50_init(&cx, mod->n);

    TMP_START;
    w = KARA_ALIGN64(TMP_ALLOC((len1 + len2 + 2 * KARA_GLEN(n) + 4 * n
                                + KARA_PAD16(len1) + fp50_mul_rec_scratch(n)
                                + 32) * sizeof(double)));
    ga = w;
    gb = ga + KARA_GLEN(n);
    p = gb + KARA_GLEN(n);
    out = p + 2 * n;
    scr = out + KARA_PAD16(len1) + 2 * n;
    ad = scr + fp50_mul_rec_scratch(n);
    bd = ad + len1;

    fp50_to_double(ad, a, len1);
    fp50_to_double(bd, b, len2);

    pb = (double *) _gap((ulong *) gb, (const ulong *) bd, len2, n);
    for (k0 = 0; k0 < len1; k0 += n)
    {
        slong la = FLINT_MIN(n, len1 - k0);
        double * d = (k0 == 0) ? out : p;

        pa = (double *) _gap((ulong *) ga, (const ulong *) ad + k0, la, n);
        if (n == 16)
            fp50_leaf_mul_short(d, pa, pb, len2, la + len2 > 17, &cx);
        else
            fp50_mul_rec(d, pa, pb, n, scr, &cx);
        if (k0 != 0)
        {
            for (i = 0; i < n; i += 4)
                vec4d_store_unaligned(out + k0 + i,
                        vec4d_add(vec4d_load_unaligned(out + k0 + i),
                                  vec4d_load_aligned(p + i)));
            for ( ; i < 2 * n; i += 4)
                vec4d_store_unaligned(out + k0 + i, vec4d_load_aligned(p + i));
        }
    }

    fp50_from_double(res, out, nout, &cx);

    TMP_END;
}

/* res[k] = sum_{i < len2} a[nlo + k - i] b[i] mod n for k < m, a zero
   outside [0, len1) */
static void
fp50_mulmid(nn_ptr res, nn_srcptr a, slong len1, nn_srcptr b, slong len2,
            slong nlo, slong m, const nmod_t * mod)
{
    fp50_ctx cx;
    slong N, k0, s, i, st, mdone = m;
    double * w, * ad, * bd, * x, * y, * r, * out, * scr;
    nn_ptr extra = NULL;
    TMP_INIT;

    fp50_init(&cx, mod->n);
    N = KARA_PAD16(FLINT_MIN(len2, m));

    TMP_START;
    w = KARA_ALIGN64(TMP_ALLOC((len1 + len2 + 4 * N + m
                                + fp50_mid_rec_scratch(N) + 16)
                           * sizeof(double)));
    x = w;
    y = x + 2 * N;
    r = y + N;
    scr = r + N;
    out = scr + fp50_mid_rec_scratch(N);
    ad = out + m;
    bd = ad + len1;

    fp50_to_double(ad, a, len1);
    fp50_to_double(bd, b, len2);

    if (m >= len2)
    {
        _gather((ulong *) y, (const ulong *) bd, len2, 0, N);
        for (k0 = 0; k0 < m; k0 += N)
        {
            if (k0 > 0 && m - k0 < N && _split_last_block(m - k0, N, len2))
            {
                mdone = k0;
                break;
            }
            st = nlo + k0 - N + 1;
            _gather((ulong *) x, (const ulong *) ad, len1, st, 2 * N - 1);
            fp50_mid_rec(r, x, y, N, scr, &cx);
            memcpy(out + k0, r, FLINT_MIN(N, m - k0) * sizeof(double));
        }
    }
    else
    {
        for (s = 0; s * N < len2; s++)
        {
            if (s > 0 && (s + 1) * N > len2
                    && _split_last_block(len2 - s * N, N, m))
            {
                extra = TMP_ALLOC(m * sizeof(ulong));
                fp50_mulmid(extra, a, len1, b + s * N, len2 - s * N,
                            nlo - s * N, m, mod);
                break;
            }
            _gather((ulong *) y, (const ulong *) bd, len2, s * N, N);
            _gather((ulong *) x, (const ulong *) ad, len1,
                    nlo - s * N - N + 1, 2 * N - 1);
            fp50_mid_rec(r, x, y, N, scr, &cx);
            if (s == 0)
                memcpy(out, r, m * sizeof(double));
            else
                for (i = 0; i < m; i++)
                    out[i] = fp50_red1(out[i] + r[i], &cx);
        }
    }

    fp50_from_double(res, out, mdone, &cx);
    if (mdone < m)
        fp50_mulmid(res + mdone, a, len1, b, len2, nlo + mdone, m - mdone,
                    mod);
    if (extra != NULL)
        _nmod_vec_add(res, res, extra, m, *mod);

    TMP_END;
}

#endif

/* ---- the rule ---------------------------------------------------------- */

/*
    Whether to take the Karatsuba for a product of the given kind (see
    impl.h for the shapes: s, l, unbalanced). The rule in one word is in
    impl.h (_nmod_poly_want_karatsuba_one_word), as the dispatchers inline
    it; the rule of the doubles is below. _nmod_poly_karatsuba_in_range
    (impl.h) checks necessary conditions first: keep it weaker than both.
*/

#if NMOD_POLY_HAVE_KARATSUBA_FP50

/*
    The rule of the doubles (n < 2^50, where one word does not hold: in
    practice from 29 bits, or for very long operands). The block length s
    must lie in [lo, hi], outside some lengths that pad badly.

    lo, the shortest blocks, by range of n (rows) and shape (columns), in
    two tables: the AVX2 rules, and the AVX-512 rules
    (FLINT_NMOD_POLY_KARA_AVX512_RULES), where the rows and the classical
    code are faster. The Karatsuba starts later where the classical code is
    fast: up to n = 1515531528 (about 2^30.5) its dot products split into
    two words; below 2^30, balanced low products where the classical code
    gives way to KS (10 + b^2/10). Middle products of few outputs
    (len2 >= 4m) start late where the classical code is fast, and so do
    unbalanced shapes, where the longer operand is cut into blocks.

    Exceptions: lengths of balanced squares and low products that cost as
    much as longer ones; with the AVX2 rules, products at 29-32 bits up to
    128 take nearly empty last blocks, but not 33-36 and 65-68; up to
    2^30.5, products and middle products of 12 blocks only from 187; with
    IFMA, which leaves the doubles what its kernel does not take, middle
    products of a much longer operand (len2 >= 8m) not up to 2^30.5 and
    from 64 above 32 bits. Below 2^28, middle products of few outputs
    (4m <= len2 < 8m) never and balanced ones only up to 128: main's are
    faster there than at 2^28-2^29, where the row was measured. Up to
    2^30.5, middle products and balanced low products below 128 not when a
    short last block makes the Karatsuba cost more than 1.2 (l/s) times the
    s x s block (main's dot products grow with l).

    hi, the longest blocks: where fft_small takes over, from the number of
    primes it needs: 192 with two (squares 128), 512 with three, from about
    46 bits (squares 256, low products 384), two thirds of that for
    unbalanced shapes, and with two primes also for 9s/8 < l <= 3s/2
    (fft_small's transform stays the same there while the Karatsuba grows);
    where it transforms n directly, 128, and for unbalanced shapes 64 if the
    longer side is less than 8 times the block or for low products, else 32;
    48 (squares 32, unbalanced 16) for its own primes. At most 192 up to
    2^30.5 (squares 128), against the classical code in two words, and with
    the AVX-512 rules squares at most 64 up to 32 bits.

    Measured in Nov 2026 on Cascade Lake (AVX-512 rules), Comet Lake (AVX2
    rules), Ice Lake, Zen 4 and a Sapphire Rapids VM.
*/

/* the ranges of n of the rule (rows of the tables) */
#define FP_N29 0    /* n < 2^29 */
#define FP_N30 1    /* 2^29 <= n < 2^30 */
#define FP_N30H 2   /* 2^30 <= n <= 1515531528 (about 2^30.5) */
#define FP_N31 3    /* 1515531528 < n < 2^31 */
#define FP_N32 4    /* 2^31 <= n < 2^32 */
#define FP_N50 5    /* 2^32 <= n < 2^50 */

/* the shapes (columns): balanced full products, squares, low and middle
   products (the kinds); unbalanced full and low products, l < 8s and
   l >= 8s; unbalanced middle products, l < 8s and l >= 8s; middle products
   of few outputs, 4m <= len2 < 8m and len2 >= 8m */
#define FP_MUL _NMOD_POLY_KIND_MUL
#define FP_SQR _NMOD_POLY_KIND_SQR
#define FP_LOW _NMOD_POLY_KIND_LOW
#define FP_MID _NMOD_POLY_KIND_MID
#define FP_LM_U8 4
#define FP_LM_UF 5
#define FP_MID_U8 6
#define FP_MID_UF 7
#define FP_FEW 8
#define FP_FEW8 9

#define FP_NEVER 10000

static slong
_fp_lo(int avx512, int row, int shape)
{
    /* X: never; 94 and 100: 10 + b^2/10 at 29 and 30 bits (below 29 bits
       one word holds up to the longest blocks of the doubles) */
    enum { X = FP_NEVER };
    static const short lo[2][6][10] = {
    {   /* AVX2 rules */
        /*              mul  sqr  low  mid   u8   uf  mu8  muf  few few8 */
        /* < 2^29    */ { 28,   X,  94, 128,  32,  32, 128, 128, 128, 128 },
        /* < 2^30    */ { 28,  28, 100,  28,  32,  32,  28,  28, 128, 128 },
        /* <= 2^30.5 */ {  0,   0,   0,   0,  16,  24,  16,  24, 128, 128 },
        /* < 2^31    */ {  0,   0,   0,   0,  16,  20,  16,  20,  16,  16 },
        /* < 2^32    */ {  0,   0,   0,   0,  16,  20,  16,  20,  16,  16 },
        /* < 2^50    */ {  0,   0,  16,   0,   0,   0,   0,   0,  16,  16 }
    },
    {   /* AVX-512 rules */
        /*              mul  sqr  low  mid   u8   uf  mu8  muf  few few8 */
        /* < 2^29    */ { 28,   X,  94, 128,  32,  32, 128, 128, 128, 128 },
        /* < 2^30    */ { 28,  28, 100,  28,  32,  32,  28,  28, 128, 128 },
        /* <= 2^30.5 */ {  0,   0,   0,   0,  24,  24,  24,  24, 128, 128 },
        /* < 2^31    */ {  0,   0,   0,   0,  24,  24,  24,  24,  24,  24 },
        /* < 2^32    */ {  0,   0,   0,   0,  24,  24,  24,  24,  24,  24 },
        /* < 2^50    */ {  0,   0,  32,   0,  12,  16,  12,  16,  16,  32 }
    }};

    return lo[avx512][row][shape];
}

/* whether s is one of the lengths of balanced squares and low products
   left out, in [a, b] (0: none) */
static int
_fp_left_out(int avx512, int row, int kind, slong s)
{
    static const short ex[2][6][2][4] = {
    {   /* AVX2 rules:      squares            low products */
        /* < 2^29    */ {{ 0,  0,  0,  0}, { 0,   0,  0,   0}},
        /* < 2^30    */ {{33, 40,  0,  0}, { 0,   0,  0,   0}},
        /* <= 2^30.5 */ {{33, 40,  0,  0}, {17,  31, 33, 105}},
        /* < 2^31    */ {{33, 40,  0,  0}, {17,  27,  0,   0}},
        /* < 2^32    */ {{33, 40,  0,  0}, {17,  27,  0,   0}},
        /* < 2^50    */ {{17, 23,  0,  0}, {17,  27, 33,  40}}
    },
    {   /* AVX-512 rules */
        /* < 2^29    */ {{ 0,  0,  0,  0}, { 0,   0,  0,   0}},
        /* < 2^30    */ {{33, 40,  0,  0}, { 0,   0,  0,   0}},
        /* <= 2^30.5 */ {{21, 23, 33, 40}, {17,  31, 33, 105}},
        /* < 2^31    */ {{21, 23, 33, 40}, {17, 112,  0,   0}},
        /* < 2^32    */ {{33, 40,  0,  0}, {17, 112,  0,   0}},
        /* < 2^50    */ {{17, 23,  0,  0}, {33,  40,  0,   0}}
    }};
    const short * e = ex[avx512][row][kind == _NMOD_POLY_KIND_LOW];

    return (s >= e[0] && s <= e[1]) || (s >= e[2] && s <= e[3]);
}

/* the padding check (_nmod_poly_karatsuba_blocks_ok) of the doubles */
static int
_fp_pad(int kind, slong bits, ulong n, slong s)
{
    if (kind == _NMOD_POLY_KIND_LOW)
        return (n <= UWORD(1515531528) || s >= 113)
            ? _NMOD_POLY_KARATSUBA_PAD_LOOSE : _NMOD_POLY_KARATSUBA_PAD_STRICT;
    if (kind == _NMOD_POLY_KIND_SQR)
        return (bits <= 35) ? _NMOD_POLY_KARATSUBA_PAD_STRICT
             : (bits <= 44) ? _NMOD_POLY_KARATSUBA_PAD_STRICT40
             : _NMOD_POLY_KARATSUBA_PAD_LOOSE;
    if (kind == _NMOD_POLY_KIND_MUL)
        return (bits <= 32) ? _NMOD_POLY_KARATSUBA_PAD_STRICT
                            : _NMOD_POLY_KARATSUBA_PAD_LOOSE;
    return (n <= UWORD(1515531528)) ? _NMOD_POLY_KARATSUBA_PAD_STRICT
                                    : _NMOD_POLY_KARATSUBA_PAD_LOOSE;
}

/* the leaves of a product or middle product whose longer dimension l is cut
   into blocks of PAD16(s), a short last block split (_split_last_block) */
static slong
_kara_leaves(slong l, slong s)
{
    slong n = KARA_PAD16(s), q = l / n, r = l - q * n, rp = KARA_PAD16(r);

    if (r == 0)
        return q * _leaves(n / 16);
    if (q > 0 && _split_last_block(r, n, s))
        return q * _leaves(n / 16) + ((s + rp - 1) / rp) * _leaves(rp / 16);
    return (q + 1) * _leaves(n / 16);
}

static int
_want_fp50(slong len1, slong len2, slong m, int kind,
           slong bits, slong s, slong l, ulong n)
{
    const int avx512 = FLINT_NMOD_POLY_KARA_AVX512_RULES;
    int unbalanced = (l >= 2 * s), near = (l < 8 * s);
    int few = (kind == _NMOD_POLY_KIND_MID && len2 >= 4 * m);
    int up_to_30h = (n <= UWORD(1515531528));
    int one_prime = (2 * bits + FLINT_BIT_COUNT(len2) <= 99);
    int row, shape, pad;
    slong lo, hi;

    row = (bits < 30) ? FP_N29 : (bits == 30) ? FP_N30
        : (bits == 31) ? (up_to_30h ? FP_N30H : FP_N31)
        : (bits == 32) ? FP_N32 : FP_N50;

    if (few)
        shape = near ? FP_FEW : FP_FEW8;
    else if (unbalanced)
        shape = (kind == _NMOD_POLY_KIND_MID)
            ? (near ? FP_MID_U8 : FP_MID_UF) : (near ? FP_LM_U8 : FP_LM_UF);
    else
        shape = kind;   /* FP_SQR, FP_LOW, FP_MUL and FP_MID are the kinds */

    lo = _fp_lo(avx512, row, shape);

#if NMOD_POLY_HAVE_IFMA
    if (shape == FP_FEW8)
    {
        if (up_to_30h)
            return 0;
        if (bits > 32)
            lo = FLINT_MAX(lo, 64);
    }
#endif

    if (s < lo)
        return 0;

    if (!unbalanced && (kind == _NMOD_POLY_KIND_SQR
                        || kind == _NMOD_POLY_KIND_LOW)
            && _fp_left_out(avx512, row, kind, s))
        return 0;

    /* below 2^28, where one word fails only for long operands, main's middle
       products are faster than at 2^28-2^29, where the row was measured */
    if (n < (UWORD(1) << 28) && kind == _NMOD_POLY_KIND_MID
            && (shape == FP_FEW || (!unbalanced && s > 128)))
        return 0;

    /* up to 2^30.5, main's dot products grow with l, while a short last
       block costs the Karatsuba a padded block: not more than 1.2 (l/s)
       times the leaves of the s x s block, for middle products and balanced
       low products below 128 */
    if (row == FP_N30H && s < 128 && (kind == _NMOD_POLY_KIND_MID
                || (kind == _NMOD_POLY_KIND_LOW && !unbalanced))
            && 5 * s * _kara_leaves(l, s)
               > 6 * l * _leaves(KARA_PAD16(s) / 16))
        return 0;

    pad = _fp_pad(kind, bits, n, s);

    if (!avx512 && kind == _NMOD_POLY_KIND_MUL && bits >= 29
            && bits <= 32 && s <= 128)
    {
        if ((s > 32 && s <= 36) || (s > 64 && s <= 68))
            return 0;
        pad = _NMOD_POLY_KARATSUBA_PAD_LOOSE;
    }

    if (!_nmod_poly_karatsuba_blocks_ok(s, kind, pad))
        return 0;

    if (up_to_30h && s > 176 && s < 187
            && (kind == _NMOD_POLY_KIND_MUL
                || kind == _NMOD_POLY_KIND_MID))
        return 0;

    /* both branches allow s <= 32: ask fft_small only beyond */
    if (s > 32 && _nmod_poly_fft_small_direct(len1, len2, n))
    {
        hi = !unbalanced ? 128
           : (near || kind == _NMOD_POLY_KIND_LOW) ? 64 : 32;
    }
    else
    {
        if (one_prime)
            hi = (kind == _NMOD_POLY_KIND_SQR) ? 128 : 192;
        else
            hi = (kind == _NMOD_POLY_KIND_SQR) ? 256
               : (kind == _NMOD_POLY_KIND_LOW) ? 384 : 512;
        if (unbalanced || (one_prime && 8 * l > 9 * s && 2 * l <= 3 * s))
            hi = 2 * hi / 3;
    }

    if (bits == 50 && s > 16 && _nmod_poly_is_fft_small_prime(n))
        hi = unbalanced ? 16 : (kind == _NMOD_POLY_KIND_SQR) ? 32 : 48;

    if (up_to_30h)
        hi = FLINT_MIN(hi, (kind == _NMOD_POLY_KIND_SQR) ? 128 : 192);

    if (avx512 && kind == _NMOD_POLY_KIND_SQR && bits <= 32)
        hi = FLINT_MIN(hi, 64);

    return s <= hi;
}

#endif

int
_nmod_poly_want_karatsuba(slong len1, slong len2, slong m, int kind, ulong n)
{
    slong bits;
    int mid = (kind == _NMOD_POLY_KIND_MID);

    bits = FLINT_BIT_COUNT(n);
    if (bits > 50 || ((len2 < 12 || m < 12)
            && !_nmod_poly_karatsuba_from_8(len2, m, kind, bits)))
        return 0;

    if (bits <= 30 && _nmod_poly_karatsuba_one_word(n,
                                mid ? FLINT_MIN(len2, m) : len2, len2, mid))
        return _nmod_poly_want_karatsuba_one_word(len1, len2, m, kind, bits, n);

#if NMOD_POLY_HAVE_KARATSUBA_FP50
    return _want_fp50(len1, len2, m, kind, bits,
                      mid ? FLINT_MIN(len2, m) : len2,
                      mid ? FLINT_MAX(len2, m)
                          : (kind == _NMOD_POLY_KIND_LOW) ? m : len1, n);
#else
    return 0;
#endif
}

/* ---- entry point ------------------------------------------------------- */

/* the one-word code where it is exact (unless force_fp50), else the
   doubles for n < 2^50, else the classical code */
static void
_mulmid_karatsuba(nn_ptr res, nn_srcptr poly1, slong len1, nn_srcptr poly2,
                  slong len2, slong nlo, slong nhi, const nmod_t * mod,
                  int force_fp50)
{
    TMP_INIT;

    FLINT_ASSERT(len1 >= 1 && len2 >= 1);
    FLINT_ASSERT(0 <= nlo && nlo < nhi && nhi <= len1 + len2 - 1);

    if (nlo == 0)
    {
        /* the full product of the operands truncated to nhi */
        slong l1 = FLINT_MIN(len1, nhi), l2 = FLINT_MIN(len2, nhi);
        ulong * c;

        if (l1 < l2)
        {
            FLINT_SWAP(nn_srcptr, poly1, poly2);
            FLINT_SWAP(slong, l1, l2);
        }

        if (!force_fp50 && _nmod_poly_karatsuba_one_word(mod->n, l2, l2, 0))
        {
            /* q in whole blocks of nb, then r; a division only for more
               than one block */
            slong nb = KARA_PAD16(l2), i;
            slong q = (l1 > nb) ? (l1 / nb) * nb : l1, r = l1 - q;

            TMP_START;
            c = KARA_ALIGN64(TMP_ALLOC((l1 + 2 * KARA_PAD16(l2) + 8)
                                       * sizeof(ulong)));
            if (r > 0 && _split_last_block(r, nb, l2))
            {
                ulong * c2 = KARA_ALIGN64(TMP_ALLOC((l2 + 2 * KARA_PAD16(r)
                                                     + 8) * sizeof(ulong)));
                _mul_raw(c, poly1, q, poly2, l2);
                _mul_raw(c2, poly2, l2, poly1 + q, r);
                for (i = 0; i < r + l2 - 1; i++)
                    c[q + i] = (i < nb) ? c[q + i] + c2[i] : c2[i];
            }
            else
                _mul_raw(c, poly1, l1, poly2, l2);
            _nmod_vec_reduce(res, c, nhi, *mod);
            TMP_END;
        }
#if NMOD_POLY_HAVE_KARATSUBA_FP50
        else if (mod->n < (UWORD(1) << 50))
        {
            slong nb = KARA_PAD16(l2), k;
            slong q = (l1 > nb) ? (l1 / nb) * nb : l1, r = l1 - q;

            if (r > 0 && _split_last_block(r, nb, l2))
            {
                slong n1 = FLINT_MIN(nhi, q + l2 - 1);

                fp50_mul(res, poly1, q, poly2, l2, n1, mod);
                if (nhi > q)
                {
                    slong n2 = FLINT_MIN(nhi - q, r + l2 - 1);
                    nn_ptr t;

                    TMP_START;
                    t = TMP_ALLOC(n2 * sizeof(ulong));
                    fp50_mul(t, poly2, l2, poly1 + q, r, n2, mod);
                    for (k = q; k < q + n2; k++)
                        res[k] = (k < n1) ? nmod_add(res[k], t[k - q], *mod)
                                          : t[k - q];
                    TMP_END;
                }
            }
            else
                fp50_mul(res, poly1, l1, poly2, l2, nhi, mod);
        }
#endif
        else
            _nmod_poly_mulmid_classical(res, poly1, l1, poly2, l2, 0, nhi, *mod);
    }
    else
    {
        slong m;
        ulong * c;

        _nmod_poly_mulmid_trim(&poly1, &len1, &poly2, &len2, &nlo, &nhi);
        m = nhi - nlo;

        if (!force_fp50 && _nmod_poly_karatsuba_one_word(mod->n,
                                            FLINT_MIN(len2, m), len2, 1))
        {
            TMP_START;
            c = TMP_ALLOC(m * sizeof(ulong));
            _mid_raw(c, poly1, len1, poly2, len2, nlo, m);
            _nmod_vec_reduce(res, c, m, *mod);
            TMP_END;
        }
#if NMOD_POLY_HAVE_KARATSUBA_FP50
        else if (mod->n < (UWORD(1) << 50))
            fp50_mulmid(res, poly1, len1, poly2, len2, nlo, m, mod);
#endif
        else
            _nmod_poly_mulmid_classical(res, poly1, len1, poly2, len2, nlo,
                                        nlo + m, *mod);
    }
}

/* mod by value, as in the dispatchers, which can then jump here with
   their own arguments (tail calls, see impl.h) */
void
_nmod_poly_mulmid_karatsuba(nn_ptr res, nn_srcptr poly1, slong len1,
                       nn_srcptr poly2, slong len2, slong nlo, slong nhi, nmod_t mod)
{
    _mulmid_karatsuba(res, poly1, len1, poly2, len2, nlo, nhi, &mod, 0);
}

void
_nmod_poly_mulmid_karatsuba_fp50(nn_ptr res, nn_srcptr poly1, slong len1,
                       nn_srcptr poly2, slong len2, slong nlo, slong nhi, nmod_t mod)
{
    _mulmid_karatsuba(res, poly1, len1, poly2, len2, nlo, nhi, &mod, 1);
}

#else

int
_nmod_poly_want_karatsuba(slong len1, slong len2, slong m, int kind, ulong n)
{
    (void) len1;
    (void) len2;
    (void) m;
    (void) kind;
    (void) n;
    return 0;
}

/* without AVX2 (or a 32-bit build): the classical code */
void
_nmod_poly_mulmid_karatsuba(nn_ptr res, nn_srcptr poly1, slong len1,
                       nn_srcptr poly2, slong len2, slong nlo, slong nhi, nmod_t mod)
{
    _nmod_poly_mulmid_classical(res, poly1, len1, poly2, len2, nlo, nhi, mod);
}

void
_nmod_poly_mulmid_karatsuba_fp50(nn_ptr res, nn_srcptr poly1, slong len1,
                       nn_srcptr poly2, slong len2, slong nlo, slong nhi, nmod_t mod)
{
    _nmod_poly_mulmid_classical(res, poly1, len1, poly2, len2, nlo, nhi, mod);
}

#endif
