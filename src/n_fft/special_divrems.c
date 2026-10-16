/*
    Copyright (C) 2025 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "n_fft/impl.h"

/*---------------------------------------*/
/* "c-circulant" division with remainder */
/*---------------------------------------*/

/* division by x**d - c, with precomputation c_precomp for mulmod by c */
/* assumes len > 0 and d > 0 */
/* lazy_4_4: in [0, 4*n) | out [0, 4*n) | max < 4n */
void _n_fft_divrem_circulant_lazy_4_4(nn_ptr p, slong len, ulong d, ulong c, ulong c_precomp, ulong n, ulong n2)
{
    ulong i, val0, val1;

    for (i = len - 1; i >= d; i--)
    {
        /* p[i-d] = p[i-d] + c * p[i] */
        val0 = p[i-d];
        val1 = p[i];
        if (val0 >= n2)
            val0 -= n2;                                      /* [0, 2n) */
        N_MULMOD_PRECOMP_LAZY(val1, c, val1, c_precomp, n);  /* [0, 2n) */
        p[i-d] = val0 + val1;                                /* [0, 4n) */
    }
}

/* same, lazy_2_2: in [0, 2*n) | out [0, 2*n) | max < 4n */
void _n_fft_divrem_circulant_lazy_2_2(nn_ptr p, slong len, ulong d, ulong c, ulong c_precomp, ulong n, ulong n2)
{
    ulong i, val0, val1;

    for (i = len - 1; i >= d; i--)
    {
        /* p[i-d] = p[i-d] + c * p[i] */
        val0 = p[i-d];
        val1 = p[i];
        N_MULMOD_PRECOMP_LAZY(val1, c, val1, c_precomp, n);  /* [0, 2n) */
        val1 += val0;                                        /* [0, 4n) */
        if (val1 >= n2)
            val1 -= n2;                                      /* [0, 2n) */
        p[i-d] = val1;
    }
}

/* division by x**d - 1 (not lazy: [0, n) -> [0, n)) */
/* assumes len >= d > 0 */
void _n_fft_divrem_circulant1(nn_ptr p, slong len, ulong d, ulong n)
{
    slong i;
    ulong j, r;

    r = len % d;
    i = len - r - d;  /* multiple of d, >= 0 by assumption */

    for (j = 0; j < r; j++)
        p[i+j] = n_addmod(p[i+j], p[d+i+j], n);

    i -= d;
    while (i >= 0)
    {
        for (j = 0; j < d; j++)
            p[i+j] = n_addmod(p[i+j], p[d+i+j], n);
        i -= d;
    }
}

/* division by x**d - 1, lazy_2_4: in [0, 2*n) | out [0, 4*n) | max < 4n */
/* assumes 0 < d < len <= 2*d (so that each p[i-d] receives a single addition) */
static void _n_fft_divrem_circulant1_lazy_2_4(nn_ptr p, slong len, ulong d)
{
    FLINT_ASSERT((ulong) len <= 2*d);
    for (ulong i = len - 1; i >= d; i--)
        p[i-d] = p[i-d] + p[i];
}

/* division by x**d - 1, lazy_2_2: in [0, 2*n) | out [0, 2*n) | max < 4n */
/* assumes len > 0 and d > 0 */
static void _n_fft_divrem_circulant1_lazy_2_2(nn_ptr p, slong len, ulong d, ulong n2)
{
    for (ulong i = len - 1; i >= d; i--)
    {
        p[i-d] = p[i-d] + p[i];
        if (p[i-d] >= n2)
            p[i-d] -= n2;
    }
}

/*---------------------------------------------------*/
/* "c-circulant" division with remainder, transposed */
/*   -> expand sequence mod x**d - c                 */
/*---------------------------------------------------*/

/* transposed version of above function:
 *    Input:
 *        :p: vector of length >= len
 *        :len: target length of expansion
 *        :d: positive integer
 *        :c: element of base field, with precomputation c_precomp
 *    Effect:
 *    find coefficients d...len-1 of p by unrolling the recurrence with
 *    charpoly x**d - c
 *    -> explicitly: for all d*i+j with j < d such that d*i+j < len,
 *        p_{d*i+j} == p_j * c**i
 *   (in particular, first d entries are unchanged)
 *    Laziness: the input entries p[0:d] may be any ulong values (they are
 *    only used as the second operand of a lazy mulmod), they are left
 *    unchanged, and the new entries p[d:len] are in [0, 2n)
 **/
void _n_fft_divrem_circulant_lazy_4_2_t(nn_ptr p, ulong len, ulong d, ulong c, ulong c_precomp, ulong n)
{
    ulong i;

    for (i = 0; i+d < len; i++)
        N_MULMOD_PRECOMP_LAZY(p[i+d], c, p[i], c_precomp, n);  /* [0, 2n) */
}

/* non-lazy, just copies */
void _n_fft_divrem_circulant1_t(nn_ptr p, ulong len, ulong d)
{
    ulong i;
    for (i = 0; i+d < len; i++)
        p[i+d] = p[i];
}


/*------------------------------------------*/
/* remainder mod product of x-w             */
/*------------------------------------------*/

/** Remainder in division modulo
 * product of x - w for specific set w of roots of unity.
 *     Input:
 *         :p: array, representing a polynomial of degree < len
 *         :len: positive integer (length of p)
 *         :d: positive integer
 *         :depth: nonnegative integer, current depth in root tree
 *         :node: nonnegative integer, current node in root tree
 *         :F: n_fft_args_t
 *     Requirements: d <= 2**depth, d even, 2**depth * node + d <= 2**F->depth
 *     Action:
 *         store in p[0:d] the coefficients of the remainder in the division
 *         of p by prod(x - w_k for k in range(d)), where the w_k are the first
 *         d evaluation points of the DFT at (depth, node) (see dft.c), namely
 *         w_{2j} = F->tab_w[2**depth * node + 2*j] and w_{2j+1} = -w_{2j}
 *         (note that the other input coefficients p[d:len] may be overwritten)
 * lazy_4_4: in [0, 4*n) | out [0, 4*n) | max < 4n
 */

/** Algorithm.
 *
 * - Base case: d is a power of 2
 * the modulus has the form x**d - w for some well-chosen w; this is just one
 * call to _n_fft_divrem_circulant_lazy_4_4
 *
 * Assume now d is not a power of 2, and suppose we have updated node and depth
 * so that 2**(depth-1) < d <= 2**depth.
 *
 * Write w = F->tab_w[2*node] and w2 = w**2, so that the evaluation points at
 * (depth, node) are the roots of x**(2**depth) - w2 (see dft.c); note that
 * w2 == F->tab_w[node] if node is even, and w2 == -F->tab_w[node-1] if node is
 * odd (this is how the code below retrieves w2 and its precomputation).
 * Write e = 2**depth / 2 (largest power of 2 strictly less than d).
 *
 * The computation has two main parts:
 *
 *  1/ If deg(p) >= 2*e, then start by reducing p mod x**(2*e) - w2
 * indeed, x**(2*e) - w2 == (x**e - w) * (x**e + w) is a multiple of the target modulus
 *     prod(x - w_k for k in range(d))
 * since the product of terms for 0 <= k < e is x**e - w,
 * and the product of terms for e <= k < d divides x**e + w
 *
 * 2/ Assuming now deg(p) < 2*e, from p == p_lo + x**e * p_hi where deg(p_lo) < e,
 * one deduces    p == (p_lo + w * p_hi) + (x**e - w) * p_hi
 * with deg(p_hi) < 2*e - e = e, so p_lo + w * p_hi is already reduced modulo x**e - w
 * So, we apply the algorithm recursively on p_hi with the d-e remaining roots
 * (node 2*node+1 at depth-1), yielding the remainder
 *     r = p_hi mod prod(x - w_k for k in range(e, d)),
 * it remains to combine: the sought remainder is
 *    (p_lo + w * p_hi) + (x**e - w) * r  ==  (p_lo + w * (p_hi - r)) + x**e * r
 */

/** Notes about possible efficiency improvements.
 *
 * - Using temporary storage would allow to compute p_hi - r and then multiply
 *   by w; instead of using two multiplications. Not clear if that will help,
 *   maybe an attempt should be considered.
 *
 * - When d is close to a power of 2, it may be more efficient to proceed
 *   otherwise. Say d == 2**k - 1, we want the remainder modulo a polynomial of
 *   the form (x**d - w2) / (x - w1). Instead of the procedure below, we may
 *   compute (x - w1) * p == q * (x**d - w2) + r (we do not need q) and then
 *   find the quotient of r by x - w1. Do we expect this to be faster than the
 *   method below? Does this generalize easily to other cases (d further from a
 *   power of 2)?
 */

void _n_fft_rem_prod_root1_node_lazy_4_4(nn_ptr p, ulong len, ulong d,
                                        ulong depth, ulong node, n_fft_args_t F)
{
    /* base case: if len <= d, polynomial is already reduced */
    if (len <= d)
        return;

    /* currently, 2 <= d <= 2**depth --> ensure 2**(depth-1) < d <= 2**depth */
    ulong depth_d = n_clog2_ge2(d);
    node = node << (depth - depth_d);
    depth = depth_d;
    ulong e = UWORD(1) << depth;

    /* if d is a power of 2, i.e. d == e, we are just reducing modulo */
    /* x**d - tab_w[2*node]**2, where tab_w[2*node]**2 == tab_w[node] if node */
    /* is even and == -tab_w[node-1] if node is odd (see the algorithm above) */
    if (d == e)
    {
        if (node % 2 == 0)
            _n_fft_divrem_circulant_lazy_4_4(p, len, d,
                                                 F->tab_w[node],
                                                 F->tab_w[node+1],
                                                 F->mod, F->mod2);
        else
            _n_fft_divrem_circulant_lazy_4_4(p, len, d,
                                                 F->mod - F->tab_w[node-1],
                                                 n_mulmod_precomp_shoup_negate(F->tab_w[node]),
                                                 F->mod, F->mod2);
        return;
    }

    e = e/2;
    /* from here on, 1 <= e == 2**(depth-1) < d <= 2**depth */

    if (len > 2*e)
    {
        if (node % 2 == 0)
            _n_fft_divrem_circulant_lazy_4_4(p, len, 2*e,
                                                 F->tab_w[node], F->tab_w[node+1],
                                                 F->mod, F->mod2);
        else
            _n_fft_divrem_circulant_lazy_4_4(p, len, 2*e,
                                                 F->mod - F->tab_w[node-1],
                                                 n_mulmod_precomp_shoup_negate(F->tab_w[node]),
                                                 F->mod, F->mod2);
        len = 2*e;
    }

    const ulong w = F->tab_w[2*node];
    const ulong w_pr = F->tab_w[2*node+1];
    ulong val0, val1;

    for (ulong i = 0; i < len - e; i++)
    {
        /* p[i] += w * p[e + i] */
        val0 = p[i];
        val1 = p[e+i];
        if (val0 >= F->mod2)
            val0 -= F->mod2;              /* [0, 2n) */
        N_MULMOD_PRECOMP_LAZY(val1, w, val1, w_pr, F->mod);  /* [0, 2n) */
        p[i] = val0 + val1;               /* [0, 4n) */
    }

    _n_fft_rem_prod_root1_node_lazy_4_4(p + e, len-e, d-e, depth-1, 2*node+1, F);

    for (ulong i = 0; i < d - e; i++)
    {
        /* p[i] -= w * p[e + i] */
        val0 = p[i];
        val1 = p[e+i];
        if (val0 >= F->mod2)
            val0 -= F->mod2;              /* [0, 2n) */
        N_MULMOD_PRECOMP_LAZY(val1, w, val1, w_pr, F->mod);  /* [0, 2n) */
        p[i] = val0 + F->mod2 - val1;     /* [0, 4n) */
    }
}

void _n_fft_rem_prod_root1_node_lazy_2_2(nn_ptr p, ulong len, ulong d,
                                        ulong depth, ulong node, n_fft_args_t F)
{
    /* base case: if len <= d, polynomial is already reduced */
    if (len <= d)
        return;

    /* currently, 2 <= d <= 2**depth --> ensure 2**(depth-1) < d <= 2**depth */
    ulong depth_d = n_clog2_ge2(d);
    node = node << (depth - depth_d);
    depth = depth_d;
    ulong e = UWORD(1) << depth;

    /* if d is a power of 2, i.e. d == e, we are just reducing modulo */
    /* x**d - tab_w[2*node]**2, where tab_w[2*node]**2 == tab_w[node] if node */
    /* is even and == -tab_w[node-1] if node is odd (see the algorithm above) */
    if (d == e)
    {
        if (node % 2 == 0)
            _n_fft_divrem_circulant_lazy_2_2(p, len, d,
                                                 F->tab_w[node],
                                                 F->tab_w[node+1],
                                                 F->mod, F->mod2);
        else
            _n_fft_divrem_circulant_lazy_2_2(p, len, d,
                                                 F->mod - F->tab_w[node-1],
                                                 n_mulmod_precomp_shoup_negate(F->tab_w[node]),
                                                 F->mod, F->mod2);
        return;
    }

    e = e/2;
    /* from here on, 1 <= e == 2**(depth-1) < d <= 2**depth */

    if (len > 2*e)
    {
        if (node % 2 == 0)
            _n_fft_divrem_circulant_lazy_2_2(p, len, 2*e,
                                                 F->tab_w[node], F->tab_w[node+1],
                                                 F->mod, F->mod2);
        else
            _n_fft_divrem_circulant_lazy_2_2(p, len, 2*e,
                                                 F->mod - F->tab_w[node-1],
                                                 n_mulmod_precomp_shoup_negate(F->tab_w[node]),
                                                 F->mod, F->mod2);
        len = 2*e;
    }

    const ulong w = F->tab_w[2*node];
    const ulong w_pr = F->tab_w[2*node+1];
    ulong val0, val1;

    for (ulong i = 0; i < len - e; i++)
    {
        /* p[i] += w * p[e + i] */
        val0 = p[i];
        val1 = p[e+i];
        N_MULMOD_PRECOMP_LAZY(val1, w, val1, w_pr, F->mod);  /* [0, 2n) */
        val1 += val0;                     /* [0, 4n) */
        if (val1 >= F->mod2)
            val1 -= F->mod2;              /* [0, 2n) */
        p[i] = val1;
    }

    _n_fft_rem_prod_root1_node_lazy_2_2(p + e, len-e, d-e, depth-1, 2*node+1, F);

    for (ulong i = 0; i < d - e; i++)
    {
        /* p[i] -= w * p[e + i] */
        val0 = p[i];
        val1 = p[e+i];
        N_MULMOD_PRECOMP_LAZY(val1, w, val1, w_pr, F->mod);  /* [0, 2n) */
        val1 = val0 + F->mod2 - val1;     /* [0, 4n) */
        if (val1 >= F->mod2)
            val1 -= F->mod2;              /* [0, 2n) */
        p[i] = val1;
    }
}

/* same as _n_fft_rem_prod_root1_node_lazy_2_2 for node == 0 (roots of x**len - 1), */
/* lazy_2_4: in [0, 2*n) | out [0, 4*n) | max < 4n */
void _n_fft_rem_prod_root1_lazy_2_4(nn_ptr p, ulong len, ulong d,
                                             ulong depth, n_fft_args_t F)
{
    /* base case: if len <= d, polynomial is already reduced */
    if (len <= d)
        return;

    /* currently, 2 <= d <= 2**depth --> ensure 2**(depth-1) < d <= 2**depth */
    ulong depth_d = n_clog2_ge2(d);
    depth = depth_d;
    ulong e = UWORD(1) << depth;

    /* if d is a power of 2, i.e. d == e, we are just reducing modulo */
    /* x**d - tab_w[2*node]**2, where tab_w[2*node]**2 == tab_w[node] if node */
    /* is even and == -tab_w[node-1] if node is odd (see the algorithm above) */
    if (d == e)
    {
        _n_fft_divrem_circulant1_lazy_2_4(p, len, d);
        return;
    }

    e = e/2;
    /* from here on, 1 <= e == 2**(depth-1) < d <= 2**depth */

    if (len > 2*e)
    {
        _n_fft_divrem_circulant1_lazy_2_2(p, len, 2*e, F->mod2);
        len = 2*e;
    }

    for (ulong i = 0; i < len - e; i++)
    {
        p[i] += p[e+i];  /* [0, 4n) */
        if (p[i] >= F->mod2)
            p[i] -= F->mod2;
    }

    _n_fft_rem_prod_root1_node_lazy_2_2(p + e, len-e, d-e, depth-1, 1, F);

    for (ulong i = 0; i < d - e; i++)
        p[i] += F->mod2 - p[e+i];  /* [0, 4n) */
}

/*---------------------------------------------*/
/* transposed remainder (recurrence extension) */
/*              mod product of x-w             */
/*---------------------------------------------*/

/** Transposed remainder (recurrence extension) modulo
 * product of x - w for specific set w of roots of unity.
 *     Input:
 *         :p: array
 *         :len: positive integer (length of p)
 *         :d: positive integer, multiple of 2
 *         :depth: nonnegative integer, current depth in root tree
 *         :node: nonnegative integer, current node in root tree
 *         :F: n_fft_args_t
 *     Requirements: d <= 2**depth, d even, 2**depth * node + d <= 2**F->depth
 *     Action:
 *         store in p[d:len] the coefficients obtained by unrolling the
 *         length-d recurrence provided by
 *            prod(x**2 - w_k**2) == prod((x-w_k)(x+w_k))
 *            for all w_k = F->tab_w[2**depth * node + 2*k], 0 <= k < d/2
 *         on the initial d coefficients p[:d]
 *         (note that the input coefficients p[d:len] are ignored and overwritten)
 *     Algorithm:
 *     one may derive explanations similar to those for _n_fft_rem_prod_root1_node_lazy_4_4
 *     or just see this as its direct transpose
 * lazy_4_4: in [0, 4*n) | out [0, 4*n) | max < 4n
 */

void _n_fft_rem_prod_root1_t_node_lazy_4_4(nn_ptr p, ulong len, ulong d,
                                          ulong depth, ulong node, n_fft_args_t F)
{
    /* base case: if len <= d, all requested terms are known */
    if (len <= d)
        return;

    /* currently, 2 <= d <= 2**depth --> ensure 2**(depth-1) < d <= 2**depth */
    ulong depth_d = n_clog2_ge2(d);
    node = node << (depth - depth_d);
    depth = depth_d;
    ulong e = UWORD(1) << depth;

    /* if d is a power of 2, i.e. d == e, we are just unrolling the recurrence */
    /* with charpoly x**d - tab_w[2*node]**2, where tab_w[2*node]**2 == tab_w[node] */
    /* if node is even and == -tab_w[node-1] if node is odd */
    if (d == e)
    {
        if (node % 2 == 0)
            _n_fft_divrem_circulant_lazy_4_2_t(p, len, d,
                                                   F->tab_w[node], F->tab_w[node+1], F->mod);
        else
            _n_fft_divrem_circulant_lazy_4_2_t(p, len, d,
                                                   F->mod - F->tab_w[node-1],
                                                   n_mulmod_precomp_shoup_negate(F->tab_w[node]),
                                                   F->mod);
        return;
    }

    ulong llen = FLINT_MIN(e, len);
    e = e/2;
    /* from here on, 1 <= e == 2**(depth-1) < d < llen <= 2**depth */

    const ulong w = F->tab_w[2*node];
    const ulong w_pr = F->tab_w[2*node+1];
    ulong val0, val1;

    for (ulong i = 0; i < d - e; i++)
    {
        /* p[e + i] -= w * p[i] */
        val0 = p[i];
        val1 = p[e+i];
        if (val1 >= F->mod2)
            val1 -= F->mod2;              /* [0, 2n) */
        N_MULMOD_PRECOMP_LAZY(val0, w, val0, w_pr, F->mod);  /* [0, 2n) */
        p[e+i] = val1 + F->mod2 - val0;   /* [0, 4n) */
    }

    _n_fft_rem_prod_root1_t_node_lazy_4_4(p + e, llen-e, d-e, depth-1, 2*node+1, F);

    for (ulong i = 0; i < llen - e; i++)
    {
        /* p[e + i] += w * p[i] */
        val0 = p[i];
        val1 = p[e+i];
        if (val1 >= F->mod2)
            val1 -= F->mod2;              /* [0, 2n) */
        N_MULMOD_PRECOMP_LAZY(val0, w, val0, w_pr, F->mod);  /* [0, 2n) */
        p[e+i] = val0 + val1;             /* [0, 4n) */
    }

    if (len > 2*e)
    {
        if (node % 2 == 0)
            _n_fft_divrem_circulant_lazy_4_2_t(p, len, 2*e,
                                                   F->tab_w[node], F->tab_w[node+1], F->mod);
        else
            _n_fft_divrem_circulant_lazy_4_2_t(p, len, 2*e,
                                                   F->mod - F->tab_w[node-1],
                                                   n_mulmod_precomp_shoup_negate(F->tab_w[node]),
                                                   F->mod);
    }
}

