/*
    Copyright (C) 2026 Fredrik Johansson
    Developed using Claude Fable 5

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h>
#include <stdio.h>
#include "longlong.h"
#include "mpn_extras.h"
#include "fmpz.h"
#include "fmpq.h"
#include "arf.h"
#include "arb.h"
#include "gr.h"
#include "gr_dft.h"

/* Fixed-point real and complex GR contexts for numerical DFTs, modeled
   on the arithmetic of nfloat/nfixed.c. An element of the real context
   consists of nlimbs + 1 contiguous limbs: a sign limb (0 or 1)
   followed by nlimbs fraction limbs representing an absolute value in
   [0, 1). An element of the complex context consists of a real and an
   imaginary part laid out contiguously (2 (nlimbs + 1) limbs), matching
   the layout assumption of the complex mode of the DFT plans.

   Additions and subtractions are exact as long as the result satisfies
   |res| < 1; as a safety net (used by the root table construction,
   where results can exceed 1 in magnitude by an ulp or two near the
   quarter points), arithmetic operations whose exact result would
   reach magnitude 1 return GR_UNABLE (leaving a wrapped value), so
   that a successful status at the end of a computation certifies that
   no overflow occurred and the error bounds hold unconditionally;
   the constants one and neg_one clamp to the largest representable
   value 1 - ulp instead of wrapping. Multiplications truncate the
   exact product toward zero, with error less than 1 ulp (bounded by
   2 ulp in all error analyses for safety). The values 1, -1, i, -i are
   not representable and are "rounded" to magnitude 1 - ulp; the DFT
   kernels never multiply by these values explicitly (they are handled
   by the free rotation classes of complex mode).

   These contexts implement only the operations required by the DFT
   machinery and basic GR functionality; they are not rings in the
   strict sense (there is no exact 1), and the predicates is_one and
   is_neg_one always return T_FALSE. */

typedef struct
{
    slong nlimbs;
}
_gr_dft_nfixed_ctx_struct;

#define NFIXED_CTX_NLIMBS(ctx) (((_gr_dft_nfixed_ctx_struct *)(ctx))->nlimbs)

/* Tuning: number of limbs from which the complex multiplication uses
   the 3-multiplication Karatsuba formula with carry-extended internal
   values (see _nfxc_mul_generic) instead of the 4-multiplication
   schoolbook formula. Provisional value: the measured crossover with
   a full-product stand-in for mulhigh is 12-16 limbs; the genuine
   flint_mpn_mulhigh_n shifts it upward since the saved multiplication
   is cheaper relative to the extra additions. */
#define GR_DFT_NFIXED_KARATSUBA_CUTOFF 22

/* Magnitude primitives *******************************************/

static void
_nfx_saturate(nn_ptr res, slong n)
{
    slong i;
    for (i = 0; i < n; i++)
        res[i] = ~UWORD(0);
}

/* res = |a| - |b| as a magnitude, returning 1 if the sign flips
   (that is, |a| < |b|); computed as a single subtraction with a
   conditional two's complement negation */
static int
_nfx_signed_sub(nn_ptr res, nn_srcptr a, nn_srcptr b, slong n)
{
    /* comparing first and subtracting in the right order beats an
       unconditional subtraction with conditional two's-complement
       negation at every generic size (see profile/p-gr_dft_nfixed):
       the comparison usually resolves within the top limbs, while
       the negation is a full second pass */
    if (mpn_cmp(a, b, n) >= 0)
    {
        mpn_sub_n(res, a, b, n);
        return 0;
    }
    else
    {
        mpn_sub_n(res, b, a, n);
        return 1;
    }
}

/* res = a + b; returns the carry (1 if the exact result has magnitude
   at least 1, in which case res holds the wrapped value). Carries are
   converted to GR_UNABLE by the operation entry points rather than
   saturating: a GR_SUCCESS status from a computation then certifies
   that all intermediate values stayed in range, making the error
   bounds of gr_dft_precomp_nfixed_bound unconditionally valid. */
static ulong
_nfx_add(nn_ptr res, nn_srcptr a, nn_srcptr b, slong n)
{
    ulong asgn = a[0], bsgn = b[0];

    if (asgn == bsgn)
    {
        res[0] = asgn;
        return mpn_add_n(res + 1, a + 1, b + 1, n);
    }
    else
    {
        res[0] = asgn ^ (ulong) _nfx_signed_sub(res + 1, a + 1, b + 1, n);
        return 0;
    }
}

static ulong
_nfx_sub(nn_ptr res, nn_srcptr a, nn_srcptr b, slong n)
{
    ulong asgn = a[0], bsgn = b[0];

    if (asgn != bsgn)
    {
        res[0] = asgn;
        return mpn_add_n(res + 1, a + 1, b + 1, n);
    }
    else
    {
        res[0] = asgn ^ (ulong) _nfx_signed_sub(res + 1, a + 1, b + 1, n);
        return 0;
    }
}

/* res = (sa, ma) + (sb, mb) for detached sign/magnitude operands,
   writing a complete sign + magnitude element; returns the carry */
static ulong
_nfx_signed_combine(nn_ptr res, ulong sa, nn_srcptr ma,
        ulong sb, nn_srcptr mb, slong n)
{
    if (sa == sb)
    {
        res[0] = sa;
        return mpn_add_n(res + 1, ma, mb, n);
    }
    else
    {
        res[0] = sa ^ (ulong) _nfx_signed_sub(res + 1, ma, mb, n);
        return 0;
    }
}

/* Fully inlined magnitude kernels for 1-4 limbs. All results may
   alias the inputs. The additions return a carry, the subtractions
   compute |a - b| with a sign-flip flag via conditional negation, and
   the multiplications return the high N limbs of the exact 2N-limb
   product (truncation toward zero, error < 1 ulp). */

/* The sums are performed by the carry chains from longlong.h (as in
   nfloat), extended by a zero high word so that the carry out of the
   magnitude lands in the top output; measured faster than portable-C
   carry recovery at every size (increasingly so with the length of
   the chain, e.g. 4.59 -> 2.82 ns at 4 limbs on x86-64). */

#define MAG_ADD_1(cy, r, a, b) \
    add_ssaaaa(cy, (r)[0], \
            UWORD(0), (a)[0], UWORD(0), (b)[0])

#define MAG_ADD_2(cy, r, a, b) \
    add_sssaaaaaa(cy, (r)[1], (r)[0], \
            UWORD(0), (a)[1], (a)[0], UWORD(0), (b)[1], (b)[0])

#define MAG_ADD_3(cy, r, a, b) \
    add_ssssaaaaaaaa(cy, (r)[2], (r)[1], (r)[0], \
            UWORD(0), (a)[2], (a)[1], (a)[0], \
            UWORD(0), (b)[2], (b)[1], (b)[0])

#define MAG_ADD_4(cy, r, a, b) \
    add_sssssaaaaaaaaaa(cy, (r)[3], (r)[2], (r)[1], (r)[0], \
            UWORD(0), (a)[3], (a)[2], (a)[1], (a)[0], \
            UWORD(0), (b)[3], (b)[2], (b)[1], (b)[0])

#define MAG_NEG_1(r) do { (r)[0] = -(r)[0]; } while (0)
#define MAG_NEG_2(r) \
    do { \
        (r)[0] = -(r)[0]; \
        (r)[1] = ~(r)[1] + ((r)[0] == 0); \
    } while (0)
#define MAG_NEG_3(r) \
    do { \
        ulong _z; \
        (r)[0] = -(r)[0]; _z = ((r)[0] == 0); \
        (r)[1] = ~(r)[1] + _z; _z &= ((r)[1] == 0); \
        (r)[2] = ~(r)[2] + _z; \
    } while (0)
#define MAG_NEG_4(r) \
    do { \
        ulong _z; \
        (r)[0] = -(r)[0]; _z = ((r)[0] == 0); \
        (r)[1] = ~(r)[1] + _z; _z &= ((r)[1] == 0); \
        (r)[2] = ~(r)[2] + _z; _z &= ((r)[2] == 0); \
        (r)[3] = ~(r)[3] + _z; \
    } while (0)

#define MAG_SSUB_1(bw, r, a, b) \
    do { \
        /* branchless ordered subtraction (measurably faster than
           subtract-then-negate at one limb; see
           profile/p-gr_dft_nfixed) */ \
        ulong _a0 = (a)[0], _b0 = (b)[0]; \
        (bw) = _a0 < _b0; \
        (r)[0] = (bw) ? _b0 - _a0 : _a0 - _b0; \
    } while (0)

/* At 2-4 limbs the differences use the borrow chains from
   longlong.h, extended by a zero high word: the top output is 0 or
   all ones according to the sign of a - b, which is normalized to a
   0/1 flip flag inside the (rare) negation branch. */

#define MAG_SSUB_2(bw, r, a, b) \
    do { \
        sub_dddmmmsss(bw, (r)[1], (r)[0], \
                UWORD(0), (a)[1], (a)[0], UWORD(0), (b)[1], (b)[0]); \
        if (bw) { (bw) = 1; MAG_NEG_2(r); } \
    } while (0)

#define MAG_SSUB_3(bw, r, a, b) \
    do { \
        sub_ddddmmmmssss(bw, (r)[2], (r)[1], (r)[0], \
                UWORD(0), (a)[2], (a)[1], (a)[0], \
                UWORD(0), (b)[2], (b)[1], (b)[0]); \
        if (bw) { (bw) = 1; MAG_NEG_3(r); } \
    } while (0)

#define MAG_SSUB_4(bw, r, a, b) \
    do { \
        sub_dddddmmmmmsssss(bw, (r)[3], (r)[2], (r)[1], (r)[0], \
                UWORD(0), (a)[3], (a)[2], (a)[1], (a)[0], \
                UWORD(0), (b)[3], (b)[2], (b)[1], (b)[0]); \
        if (bw) { (bw) = 1; MAG_NEG_4(r); } \
    } while (0)

/* r must not alias a, b */
#define MAG_MULHI_1(r, a, b) \
    do { \
        ulong _lo; \
        umul_ppmm((r)[0], _lo, (a)[0], (b)[0]); \
        (void) _lo; \
    } while (0)

#define MAG_MULHI_2(r, a, b) \
    do { \
        ulong _p1, _p0; \
        FLINT_MPN_MUL_2X2((r)[1], (r)[0], _p1, _p0, \
                (a)[1], (a)[0], (b)[1], (b)[0]); \
        (void) _p1; (void) _p0; \
    } while (0)

#define MAG_MULHI_3(r, a, b) flint_mpn_mulhigh_n((r), (a), (b), 3)
#define MAG_MULHI_4(r, a, b) flint_mpn_mulhigh_n((r), (a), (b), 4)

/* Real ring methods (generic parts shared by all limb counts) ****/

static int
_nfx_add_op(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
{
    /* branch-free overflow detection: GR_UNABLE == 2 */
    return (int) (_nfx_add(res, x, y, NFIXED_CTX_NLIMBS(ctx)) << 1);
}

static int
_nfx_sub_op(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
{
    return (int) (_nfx_sub(res, x, y, NFIXED_CTX_NLIMBS(ctx)) << 1);
}

/* generic multiplication: mulhigh with truncation toward the exact
   high product; error at most 2 ulp */
static int
_nfx_mul_op(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
{
    slong n = NFIXED_CTX_NLIMBS(ctx);
    nn_ptr r = res;
    nn_srcptr a = x, b = y;

    r[0] = a[0] ^ b[0];

    if (res != x && res != y)
    {
        flint_mpn_mulhigh_n(r + 1, a + 1, b + 1, n);
    }
    else
    {
        ulong t[GR_DFT_NFIXED_MAX_NLIMBS];
        slong i;
        flint_mpn_mulhigh_n(t, a + 1, b + 1, n);
        for (i = 0; i < n; i++)
            r[i + 1] = t[i];
    }

    return GR_SUCCESS;
}

/* division by an unsigned integer: floor division of the magnitude,
   hence truncation toward zero with error below 1 ulp */
static int
_nfx_div_ui_op(gr_ptr res, gr_srcptr x, ulong y, gr_ctx_t ctx)
{
    slong n = NFIXED_CTX_NLIMBS(ctx);
    nn_ptr r = res;
    nn_srcptr a = x;

    if (y == 0)
        return GR_DOMAIN;

    r[0] = a[0];
    mpn_divrem_1(r + 1, 0, a + 1, n, y);

    return GR_SUCCESS;
}

/* squaring: sqrhigh saves about half the limb products */
static int
_nfx_sqr_op(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
{
    slong n = NFIXED_CTX_NLIMBS(ctx);
    nn_ptr r = res;
    nn_srcptr a = x;

    r[0] = 0;

    if (res != x)
    {
        flint_mpn_sqrhigh(r + 1, a + 1, n);
    }
    else
    {
        ulong t[GR_DFT_NFIXED_MAX_NLIMBS];
        slong i;
        flint_mpn_sqrhigh(t, a + 1, n);
        for (i = 0; i < n; i++)
            r[i + 1] = t[i];
    }

    return GR_SUCCESS;
}

/* Sized real methods for 1-4 limbs *******************************/

#define DEF_NFX_SIZED(N) \
static int \
_nfx_add_op_##N(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx) \
{ \
    nn_ptr r = res; \
    nn_srcptr a = x, b = y; \
    if (a[0] == b[0]) \
    { \
        ulong cy; \
        r[0] = a[0]; \
        MAG_ADD_##N(cy, r + 1, a + 1, b + 1); \
        return (int) (cy << 1);   /* GR_UNABLE on overflow */ \
    } \
    else \
    { \
        ulong bw; \
        MAG_SSUB_##N(bw, r + 1, a + 1, b + 1); \
        r[0] = a[0] ^ bw; \
    } \
    return GR_SUCCESS; \
} \
static int \
_nfx_sub_op_##N(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx) \
{ \
    nn_ptr r = res; \
    nn_srcptr a = x, b = y; \
    if (a[0] != b[0]) \
    { \
        ulong cy; \
        r[0] = a[0]; \
        MAG_ADD_##N(cy, r + 1, a + 1, b + 1); \
        return (int) (cy << 1); \
    } \
    else \
    { \
        ulong bw; \
        MAG_SSUB_##N(bw, r + 1, a + 1, b + 1); \
        r[0] = a[0] ^ bw; \
    } \
    return GR_SUCCESS; \
} \
static int \
_nfx_mul_op_##N(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx) \
{ \
    nn_ptr r = res; \
    nn_srcptr a = x, b = y; \
    ulong t[N]; \
    slong i; \
    MAG_MULHI_##N(t, a + 1, b + 1); \
    r[0] = a[0] ^ b[0]; \
    for (i = 0; i < N; i++) \
        r[i + 1] = t[i]; \
    return GR_SUCCESS; \
}

DEF_NFX_SIZED(1)
DEF_NFX_SIZED(2)
DEF_NFX_SIZED(3)
DEF_NFX_SIZED(4)

/* Complex multiplication *****************************************/

#define NFXC_IM(x, n) ((nn_ptr) (x) + ((n) + 1))
#define NFXC_IM_SRC(x, n) ((nn_srcptr) (x) + ((n) + 1))

static int
_nfxc_div_ui_op(gr_ptr res, gr_srcptr x, ulong y, gr_ctx_t ctx)
{
    slong n = NFIXED_CTX_NLIMBS(ctx);
    int status = GR_SUCCESS;

    status |= _nfx_div_ui_op(res, x, y, ctx);
    status |= _nfx_div_ui_op(NFXC_IM(res, n), NFXC_IM_SRC(x, n), y, ctx);

    return status;
}


/* Sized schoolbook complex multiplication: all reads happen into
   stack temporaries before any write, so aliasing is safe with no
   copies of the result. Per-component error at most 2 ulp
   (two exact-truncation products). */
#define DEF_NFXC_MUL_SIZED(N) \
static int \
_nfxc_mul_op_##N(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx) \
{ \
    nn_srcptr xr = x, xi = NFXC_IM_SRC(x, N); \
    nn_srcptr yr = y, yi = NFXC_IM_SRC(y, N); \
    nn_ptr rr = res, ri = NFXC_IM(res, N); \
    ulong t1[N], t2[N], t3[N], t4[N]; \
    ulong s1, s2, s3, s4, f; \
    s1 = xr[0] ^ yr[0]; MAG_MULHI_##N(t1, xr + 1, yr + 1); \
    s2 = xi[0] ^ yi[0]; MAG_MULHI_##N(t2, xi + 1, yi + 1); \
    s3 = xr[0] ^ yi[0]; MAG_MULHI_##N(t3, xr + 1, yi + 1); \
    s4 = xi[0] ^ yr[0]; MAG_MULHI_##N(t4, xi + 1, yr + 1); \
    ulong ov = 0; \
    /* re = t1 - t2 */ \
    if (s1 != s2) \
    { \
        rr[0] = s1; \
        MAG_ADD_##N(ov, rr + 1, t1, t2); \
    } \
    else \
    { \
        MAG_SSUB_##N(f, rr + 1, t1, t2); \
        rr[0] = s1 ^ f; \
    } \
    /* im = t3 + t4 */ \
    if (s3 == s4) \
    { \
        ulong cy; \
        ri[0] = s3; \
        MAG_ADD_##N(cy, ri + 1, t3, t4); \
        ov |= cy; \
    } \
    else \
    { \
        MAG_SSUB_##N(f, ri + 1, t3, t4); \
        ri[0] = s3 ^ f; \
    } \
    return (int) (ov << 1); \
}

DEF_NFXC_MUL_SIZED(1)
DEF_NFXC_MUL_SIZED(2)
DEF_NFXC_MUL_SIZED(3)
DEF_NFXC_MUL_SIZED(4)

/* Generic complex multiplication. Below the Karatsuba cutoff, the
   schoolbook formula with 4 mulhighs (per-component error at most
   4 ulp: two products of at most 2 ulp each, combined exactly).

   From the cutoff onward, the complex Karatsuba formula

       re = t1 - t2,  im = (xr + xi)(yr + yi) - t1 - t2,
       t1 = xr yr, t2 = xi yi

   with 3 mulhighs. The sums s = xr + xi (and likewise for y) can
   reach 2 in magnitude and are kept as carry-extended sign-magnitude
   values with an extra high limb: s = h 2^(64 n) + f, h in {0, 1}.
   The middle product is then computed with a single mulhigh on the
   fraction parts plus conditional exact additions,

       s1 s2 / 2^(64 n) = h1 h2 2^(64 n) + h1 f2 + h2 f1
                          + mulhigh(f1, f2) + O(2 ulp),

   whose magnitude is below 4 and hence fits an (n+1)-limb extended
   value (top limb at most 3). The final combination im = m - t1 - t2
   runs in exact (n+1)-limb sign-magnitude arithmetic (an extra high
   limb for the internal operations), and the result is stored back
   to n limbs, saturating in case the caller violated the |t| < 1
   contract. Per-component error at most 6 ulp (three mulhighs of at
   most 2 ulp each, all other operations exact). */
static int
_nfxc_mul_impl(gr_ptr res, gr_srcptr x, gr_srcptr y, int use_karatsuba,
        gr_ctx_t ctx)
{
    slong n = NFIXED_CTX_NLIMBS(ctx);
    nn_srcptr xr = x, xi = NFXC_IM_SRC(x, n);
    nn_srcptr yr = y, yi = NFXC_IM_SRC(y, n);
    nn_ptr rr = res, ri = NFXC_IM(res, n);

    if (!use_karatsuba)
    {
        ulong t1[GR_DFT_NFIXED_MAX_NLIMBS];
        ulong t2[GR_DFT_NFIXED_MAX_NLIMBS];
        ulong t3[GR_DFT_NFIXED_MAX_NLIMBS];
        ulong t4[GR_DFT_NFIXED_MAX_NLIMBS];
        ulong s1, s2, s3, s4;

        s1 = xr[0] ^ yr[0]; flint_mpn_mulhigh_n(t1, xr + 1, yr + 1, n);
        s2 = xi[0] ^ yi[0]; flint_mpn_mulhigh_n(t2, xi + 1, yi + 1, n);
        s3 = xr[0] ^ yi[0]; flint_mpn_mulhigh_n(t3, xr + 1, yi + 1, n);
        s4 = xi[0] ^ yr[0]; flint_mpn_mulhigh_n(t4, xi + 1, yr + 1, n);

        ulong ov;

        ov = _nfx_signed_combine(rr, s1, t1, s2 ^ 1, t2, n);
        ov |= _nfx_signed_combine(ri, s3, t3, s4, t4, n);
        return (int) (ov << 1);
    }
    else
    {
        /* sign + (n+1)-limb extended magnitudes */
        ulong es1[GR_DFT_NFIXED_MAX_NLIMBS + 2];
        ulong es2[GR_DFT_NFIXED_MAX_NLIMBS + 2];
        ulong em[GR_DFT_NFIXED_MAX_NLIMBS + 2];
        ulong et[GR_DFT_NFIXED_MAX_NLIMBS + 2];
        ulong t1[GR_DFT_NFIXED_MAX_NLIMBS];
        ulong t2[GR_DFT_NFIXED_MAX_NLIMBS];
        ulong s1, s2, top;
        slong i;

        /* s1 = xr + xi, s2 = yr + yi */
        if (xr[0] == xi[0])
        {
            es1[0] = xr[0];
            es1[n + 1] = mpn_add_n(es1 + 1, xr + 1, xi + 1, n);
        }
        else
        {
            es1[0] = xr[0] ^ (ulong) _nfx_signed_sub(es1 + 1, xr + 1, xi + 1, n);
            es1[n + 1] = 0;
        }
        if (yr[0] == yi[0])
        {
            es2[0] = yr[0];
            es2[n + 1] = mpn_add_n(es2 + 1, yr + 1, yi + 1, n);
        }
        else
        {
            es2[0] = yr[0] ^ (ulong) _nfx_signed_sub(es2 + 1, yr + 1, yi + 1, n);
            es2[n + 1] = 0;
        }

        s1 = xr[0] ^ yr[0]; flint_mpn_mulhigh_n(t1, xr + 1, yr + 1, n);
        s2 = xi[0] ^ yi[0]; flint_mpn_mulhigh_n(t2, xi + 1, yi + 1, n);

        /* m = s1 s2 = h1 h2 + h1 f2 + h2 f1 + mulhigh(f1, f2) */
        em[0] = es1[0] ^ es2[0];
        flint_mpn_mulhigh_n(em + 1, es1 + 1, es2 + 1, n);
        top = es1[n + 1] & es2[n + 1];
        if (es1[n + 1])
            top += mpn_add_n(em + 1, em + 1, es2 + 1, n);
        if (es2[n + 1])
            top += mpn_add_n(em + 1, em + 1, es1 + 1, n);
        em[n + 1] = top;

        /* re = t1 - t2 */
        ulong ov = _nfx_signed_combine(rr, s1, t1, s2 ^ 1, t2, n);

        /* im = m - t1 - t2 in (n+1)-limb arithmetic */
        et[0] = s1 ^ 1;
        for (i = 0; i < n; i++)
            et[i + 1] = t1[i];
        et[n + 1] = 0;
        (void) _nfx_add(em, em, et, n + 1);   /* cannot carry */
        et[0] = s2 ^ 1;
        for (i = 0; i < n; i++)
            et[i + 1] = t2[i];
        (void) _nfx_add(em, em, et, n + 1);   /* cannot carry */

        ri[0] = em[0];
        ov |= (em[n + 1] != 0);
        for (i = 0; i < n; i++)
            ri[i + 1] = em[i + 1];

        return (int) (ov << 1);
    }
}

static int
_nfxc_mul_generic(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
{
    return _nfxc_mul_impl(res, x, y,
            NFIXED_CTX_NLIMBS(ctx) >= GR_DFT_NFIXED_KARATSUBA_CUTOFF, ctx);
}

/* forced-path entry points for tuning (see profile/p-gr_dft_nfixed) */
int
_gr_dft_nfixed_cmul_schoolbook(gr_ptr res, gr_srcptr x, gr_srcptr y,
        gr_ctx_t ctx)
{
    return _nfxc_mul_impl(res, x, y, 0, ctx);
}

int
_gr_dft_nfixed_cmul_karatsuba(gr_ptr res, gr_srcptr x, gr_srcptr y,
        gr_ctx_t ctx)
{
    return _nfxc_mul_impl(res, x, y, 1, ctx);
}

/* Complex squaring with 2 real multiplications:

       re = (a + b)(a - b),   im = 2 a b,

   where the sums a +- b (of magnitude up to sqrt(2)) are held as
   carry-extended sign-magnitude values with an extra high limb, and
   their product is computed with a single mulhigh on the fraction
   parts plus conditional exact additions selected by the carries, as
   in the Karatsuba multiplication above. The doubling of a b is an
   exact addition (saturating at the representation limit if the
   caller violated the |t| < 1 contract). Per-component errors: at
   most 2 ulp on the real part (one mulhigh), at most 4 ulp on the
   imaginary part (a doubled mulhigh); both within the error constant
   of the general complex multiplication. */
static int
_nfxc_sqr_generic(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
{
    slong n = NFIXED_CTX_NLIMBS(ctx);
    nn_srcptr a = x, b = NFXC_IM_SRC(x, n);
    nn_ptr rr = res, ri = NFXC_IM(res, n);
    ulong es1[GR_DFT_NFIXED_MAX_NLIMBS + 2];
    ulong es2[GR_DFT_NFIXED_MAX_NLIMBS + 2];
    ulong em[GR_DFT_NFIXED_MAX_NLIMBS + 2];
    ulong t[GR_DFT_NFIXED_MAX_NLIMBS];
    ulong sim, top;
    slong i;

    /* s1 = a + b, s2 = a - b, carry-extended */
    if (a[0] == b[0])
    {
        es1[0] = a[0];
        es1[n + 1] = mpn_add_n(es1 + 1, a + 1, b + 1, n);
        es2[0] = a[0] ^ (ulong) _nfx_signed_sub(es2 + 1, a + 1, b + 1, n);
        es2[n + 1] = 0;
    }
    else
    {
        es1[0] = a[0] ^ (ulong) _nfx_signed_sub(es1 + 1, a + 1, b + 1, n);
        es1[n + 1] = 0;
        es2[0] = a[0];
        es2[n + 1] = mpn_add_n(es2 + 1, a + 1, b + 1, n);
    }

    /* im = 2 a b, before res is written */
    sim = a[0] ^ b[0];
    flint_mpn_mulhigh_n(t, a + 1, b + 1, n);

    /* re = s1 s2 */
    rr[0] = es1[0] ^ es2[0];
    flint_mpn_mulhigh_n(em + 1, es1 + 1, es2 + 1, n);
    top = es1[n + 1] & es2[n + 1];
    if (es1[n + 1])
        top += mpn_add_n(em + 1, em + 1, es2 + 1, n);
    if (es2[n + 1])
        top += mpn_add_n(em + 1, em + 1, es1 + 1, n);
    for (i = 0; i < n; i++)
        rr[i + 1] = em[i + 1];

    ri[0] = sim;

    /* overflow of either component: GR_UNABLE, wrapped values */
    return (int) ((((ulong) (top != 0)) |
            mpn_add_n(ri + 1, t, t, n)) << 1);
}

/* per-component ulp error bound of one complex multiplication in the
   given context (used by the root table construction and the DFT
   bound computation) */
double
_gr_dft_nfixed_cmul_err_ulps(gr_ctx_t ctx)
{
    slong n = NFIXED_CTX_NLIMBS(ctx);

    if (n <= 2)
        return 2.0;   /* exact-truncation products */
    if (n < GR_DFT_NFIXED_KARATSUBA_CUTOFF)
        return 4.0;   /* schoolbook, 2 mulhighs per component */
    return 6.0;       /* Karatsuba, 3 mulhighs on the im part */
}

/* Shared scalar methods (all limb counts) ************************/

static int
_nfx_ctx_write(gr_stream_t out, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    status |= gr_stream_write(out, "Fixed-point numbers x, |x| < 1, with ");
    status |= gr_stream_write_si(out, NFIXED_CTX_NLIMBS(ctx));
    status |= gr_stream_write(out, "-limb precision (gr_dft_nfixed)");
    return status;
}

static void
_nfx_init(gr_ptr x, gr_ctx_t ctx)
{
    memset(x, 0, ctx->sizeof_elem);
}

static void
_nfx_clear(gr_ptr x, gr_ctx_t ctx)
{
}

static void
_nfx_swap(gr_ptr x, gr_ptr y, gr_ctx_t ctx)
{
    slong i, len = ctx->sizeof_elem / sizeof(ulong);
    nn_ptr xp = x, yp = y;
    for (i = 0; i < len; i++)
    {
        ulong t = xp[i];
        xp[i] = yp[i];
        yp[i] = t;
    }
}

static int
_nfx_set(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
{
    if (res != x)
        memcpy(res, x, ctx->sizeof_elem);
    return GR_SUCCESS;
}

static int
_nfx_zero(gr_ptr res, gr_ctx_t ctx)
{
    memset(res, 0, ctx->sizeof_elem);
    return GR_SUCCESS;
}

static int
_nfx_one(gr_ptr res, gr_ctx_t ctx)
{
    nn_ptr r = res;
    r[0] = 0;
    _nfx_saturate(r + 1, NFIXED_CTX_NLIMBS(ctx));
    return GR_SUCCESS;
}

static int
_nfx_neg_one(gr_ptr res, gr_ctx_t ctx)
{
    nn_ptr r = res;
    r[0] = 1;
    _nfx_saturate(r + 1, NFIXED_CTX_NLIMBS(ctx));
    return GR_SUCCESS;
}

static int
_nfx_set_si(gr_ptr res, slong x, gr_ctx_t ctx)
{
    if (x == 0)
        return _nfx_zero(res, ctx);
    if (x == 1)
        return _nfx_one(res, ctx);
    if (x == -1)
        return _nfx_neg_one(res, ctx);
    return GR_DOMAIN;
}

static int
_nfx_set_ui(gr_ptr res, ulong x, gr_ctx_t ctx)
{
    if (x == 0)
        return _nfx_zero(res, ctx);
    if (x == 1)
        return _nfx_one(res, ctx);
    return GR_DOMAIN;
}

static int
_nfx_neg(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
{
    if (res != x)
        memcpy(res, x, ctx->sizeof_elem);
    ((nn_ptr) res)[0] ^= 1;
    return GR_SUCCESS;
}

static int
_nfx_is_zero(nn_srcptr a, slong n)
{
    slong i;
    for (i = 1; i <= n; i++)
        if (a[i] != 0)
            return 0;
    return 1;
}

static truth_t
_nfx_is_zero_op(gr_srcptr x, gr_ctx_t ctx)
{
    return _nfx_is_zero(x, NFIXED_CTX_NLIMBS(ctx)) ? T_TRUE : T_FALSE;
}

static truth_t
_nfx_false_predicate(gr_srcptr x, gr_ctx_t ctx)
{
    return T_FALSE;
}

static truth_t
_nfx_equal(gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
{
    slong n = NFIXED_CTX_NLIMBS(ctx);
    nn_srcptr xp = x, yp = y;

    if (_nfx_is_zero(xp, n))
        return _nfx_is_zero(yp, n) ? T_TRUE : T_FALSE;

    if (xp[0] != yp[0])
        return T_FALSE;

    return (mpn_cmp(xp + 1, yp + 1, n) == 0) ? T_TRUE : T_FALSE;
}

/* writes exactly n + 1 limbs (one real part) */
static void
_nfx_randtest_limbs(nn_ptr r, flint_rand_t state, slong n)
{
    slong i;

    if (n_randint(state, 8) == 0)
    {
        for (i = 0; i <= n; i++)
            r[i] = 0;
        return;
    }

    r[0] = n_randint(state, 2);
    for (i = 1; i <= n; i++)
        r[i] = n_randtest(state);
    /* vary the magnitude */
    r[n] >>= n_randint(state, FLINT_BITS);
}

static int
_nfx_randtest(gr_ptr res, flint_rand_t state, gr_ctx_t ctx)
{
    _nfx_randtest_limbs(res, state, NFIXED_CTX_NLIMBS(ctx));
    return GR_SUCCESS;
}

void
gr_dft_nfixed_get_arf(arf_t res, gr_srcptr x, gr_ctx_t ctx)
{
    slong n = NFIXED_CTX_NLIMBS(ctx);
    nn_srcptr xp = x;
    fmpz_t t;

    fmpz_init(t);
    fmpz_set_ui_array(t, xp + 1, n);
    if (xp[0])
        fmpz_neg(t, t);
    arf_set_fmpz(res, t);
    arf_mul_2exp_si(res, res, -n * FLINT_BITS);
    fmpz_clear(t);
}

int
gr_dft_nfixed_set_arf(gr_ptr res, const arf_t x, gr_ctx_t ctx)
{
    slong n = NFIXED_CTX_NLIMBS(ctx);
    nn_ptr r = res;
    arf_t t;
    fmpz_t z;

    if (!arf_is_finite(x))
        return GR_DOMAIN;

    if (arf_cmpabs_2exp_si(x, 0) >= 0)
    {
        /* saturate values with |x| >= 1 */
        r[0] = (arf_sgn(x) < 0);
        _nfx_saturate(r + 1, n);
        return GR_SUCCESS;
    }

    arf_init(t);
    fmpz_init(z);

    arf_mul_2exp_si(t, x, n * FLINT_BITS);
    arf_get_fmpz(z, t, ARF_RND_DOWN);   /* toward zero */

    r[0] = (fmpz_sgn(z) < 0);
    fmpz_abs(z, z);
    fmpz_get_ui_array(r + 1, n, z);

    arf_clear(t);
    fmpz_clear(z);

    return GR_SUCCESS;
}

static int
_nfx_write(gr_stream_t out, gr_srcptr x, gr_ctx_t ctx)
{
    arf_t t;
    char buf[64];

    int status;

    arf_init(t);
    gr_dft_nfixed_get_arf(t, x, ctx);
    snprintf(buf, sizeof(buf), "%.17g", arf_get_d(t, ARF_RND_NEAR));
    status = gr_stream_write(out, buf);
    arf_clear(t);
    return status;
}

static truth_t
_nfx_ctx_predicate_true(gr_ctx_t ctx)
{
    return T_TRUE;
}

static int
_nfx_ctx_get_real_prec(slong * res, gr_ctx_t ctx)
{
    *res = NFIXED_CTX_NLIMBS(ctx) * FLINT_BITS;
    return GR_SUCCESS;
}

/* Complex shared methods *****************************************/

static int
_nfxc_ctx_write(gr_stream_t out, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    status |= gr_stream_write(out, "Complex fixed-point numbers x, |re(x)|, |im(x)| < 1, with ");
    status |= gr_stream_write_si(out, NFIXED_CTX_NLIMBS(ctx));
    status |= gr_stream_write(out, "-limb precision (gr_dft_nfixed_complex)");
    return status;
}

static int
_nfxc_one(gr_ptr res, gr_ctx_t ctx)
{
    slong n = NFIXED_CTX_NLIMBS(ctx);
    nn_ptr r = res;
    r[0] = 0;
    _nfx_saturate(r + 1, n);
    memset(NFXC_IM(res, n), 0, (n + 1) * sizeof(ulong));
    return GR_SUCCESS;
}

static int
_nfxc_neg_one(gr_ptr res, gr_ctx_t ctx)
{
    slong n = NFIXED_CTX_NLIMBS(ctx);
    nn_ptr r = res;
    r[0] = 1;
    _nfx_saturate(r + 1, n);
    memset(NFXC_IM(res, n), 0, (n + 1) * sizeof(ulong));
    return GR_SUCCESS;
}

static int
_nfxc_set_si(gr_ptr res, slong x, gr_ctx_t ctx)
{
    if (x == 0)
        return _nfx_zero(res, ctx);
    if (x == 1)
        return _nfxc_one(res, ctx);
    if (x == -1)
        return _nfxc_neg_one(res, ctx);
    return GR_DOMAIN;
}

static int
_nfxc_set_ui(gr_ptr res, ulong x, gr_ctx_t ctx)
{
    if (x == 0)
        return _nfx_zero(res, ctx);
    if (x == 1)
        return _nfxc_one(res, ctx);
    return GR_DOMAIN;
}

static int
_nfxc_neg(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
{
    slong n = NFIXED_CTX_NLIMBS(ctx);
    if (res != x)
        memcpy(res, x, ctx->sizeof_elem);
    ((nn_ptr) res)[0] ^= 1;
    NFXC_IM(res, n)[0] ^= 1;
    return GR_SUCCESS;
}

static int
_nfxc_add_op(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
{
    slong n = NFIXED_CTX_NLIMBS(ctx);
    ulong ov;
    ov = _nfx_add(res, x, y, n);
    ov |= _nfx_add(NFXC_IM(res, n), NFXC_IM_SRC(x, n), NFXC_IM_SRC(y, n), n);
    return (int) (ov << 1);
}

static int
_nfxc_sub_op(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
{
    slong n = NFIXED_CTX_NLIMBS(ctx);
    ulong ov;
    ov = _nfx_sub(res, x, y, n);
    ov |= _nfx_sub(NFXC_IM(res, n), NFXC_IM_SRC(x, n), NFXC_IM_SRC(y, n), n);
    return (int) (ov << 1);
}

/* sized componentwise complex add/sub via the sized real methods */
#define DEF_NFXC_ADDSUB_SIZED(N) \
static int \
_nfxc_add_op_##N(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx) \
{ \
    int status; \
    status = _nfx_add_op_##N(res, x, y, ctx); \
    status |= _nfx_add_op_##N(NFXC_IM(res, N), NFXC_IM_SRC(x, N), \
            NFXC_IM_SRC(y, N), ctx); \
    return status; \
} \
static int \
_nfxc_sub_op_##N(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx) \
{ \
    int status; \
    status = _nfx_sub_op_##N(res, x, y, ctx); \
    status |= _nfx_sub_op_##N(NFXC_IM(res, N), NFXC_IM_SRC(x, N), \
            NFXC_IM_SRC(y, N), ctx); \
    return status; \
}

DEF_NFXC_ADDSUB_SIZED(1)
DEF_NFXC_ADDSUB_SIZED(2)
DEF_NFXC_ADDSUB_SIZED(3)
DEF_NFXC_ADDSUB_SIZED(4)

static truth_t
_nfxc_is_zero_op(gr_srcptr x, gr_ctx_t ctx)
{
    slong n = NFIXED_CTX_NLIMBS(ctx);
    return (_nfx_is_zero(x, n) && _nfx_is_zero(NFXC_IM_SRC(x, n), n))
        ? T_TRUE : T_FALSE;
}

static truth_t
_nfxc_equal(gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
{
    truth_t r1, r2;
    slong n = NFIXED_CTX_NLIMBS(ctx);

    r1 = _nfx_equal(x, y, ctx);
    if (r1 == T_FALSE)
        return T_FALSE;
    r2 = _nfx_equal(NFXC_IM_SRC(x, n), NFXC_IM_SRC(y, n), ctx);
    return (r1 == T_TRUE && r2 == T_TRUE) ? T_TRUE :
           (r2 == T_FALSE ? T_FALSE : T_UNKNOWN);
}

static int
_nfxc_randtest(gr_ptr res, flint_rand_t state, gr_ctx_t ctx)
{
    slong n = NFIXED_CTX_NLIMBS(ctx);
    _nfx_randtest_limbs(res, state, n);
    _nfx_randtest_limbs(NFXC_IM(res, n), state, n);
    return GR_SUCCESS;
}

static int
_nfxc_write(gr_stream_t out, gr_srcptr x, gr_ctx_t ctx)
{
    slong n = NFIXED_CTX_NLIMBS(ctx);
    int status = GR_SUCCESS;
    status |= _nfx_write(out, x, ctx);
    status |= gr_stream_write(out, " + ");
    status |= _nfx_write(out, NFXC_IM_SRC(x, n), ctx);
    status |= gr_stream_write(out, "*I");
    return status;
}

/* Method tables: one per hardcoded limb count 1-4 plus generic ****/

#define NFX_TAB_COUNT 5

static gr_static_method_table _nfx_methods_tab[NFX_TAB_COUNT];
static gr_static_method_table _nfxc_methods_tab[NFX_TAB_COUNT];
static int _nfx_methods_initialized = 0;

#define NFX_COMMON_METHODS \
    {GR_METHOD_CTX_WRITE,           (gr_funcptr) _nfx_ctx_write}, \
    {GR_METHOD_CTX_IS_THREADSAFE,   (gr_funcptr) _nfx_ctx_predicate_true}, \
    {GR_METHOD_CTX_GET_REAL_PREC,   (gr_funcptr) _nfx_ctx_get_real_prec}, \
    {GR_METHOD_INIT,                (gr_funcptr) _nfx_init}, \
    {GR_METHOD_CLEAR,               (gr_funcptr) _nfx_clear}, \
    {GR_METHOD_SWAP,                (gr_funcptr) _nfx_swap}, \
    {GR_METHOD_RANDTEST,            (gr_funcptr) _nfx_randtest}, \
    {GR_METHOD_WRITE,               (gr_funcptr) _nfx_write}, \
    {GR_METHOD_ZERO,                (gr_funcptr) _nfx_zero}, \
    {GR_METHOD_ONE,                 (gr_funcptr) _nfx_one}, \
    {GR_METHOD_NEG_ONE,             (gr_funcptr) _nfx_neg_one}, \
    {GR_METHOD_IS_ZERO,             (gr_funcptr) _nfx_is_zero_op}, \
    {GR_METHOD_IS_ONE,              (gr_funcptr) _nfx_false_predicate}, \
    {GR_METHOD_IS_NEG_ONE,          (gr_funcptr) _nfx_false_predicate}, \
    {GR_METHOD_EQUAL,               (gr_funcptr) _nfx_equal}, \
    {GR_METHOD_SET,                 (gr_funcptr) _nfx_set}, \
    {GR_METHOD_SET_SI,              (gr_funcptr) _nfx_set_si}, \
    {GR_METHOD_SET_UI,              (gr_funcptr) _nfx_set_ui}, \
    {GR_METHOD_DIV_UI,              (gr_funcptr) _nfx_div_ui_op}, \
    {GR_METHOD_NEG,                 (gr_funcptr) _nfx_neg}

#define DEF_NFX_INPUT(N) \
static gr_method_tab_input _nfx_methods_input_##N[] = \
{ \
    NFX_COMMON_METHODS, \
    {GR_METHOD_ADD,                 (gr_funcptr) _nfx_add_op_##N}, \
    {GR_METHOD_SUB,                 (gr_funcptr) _nfx_sub_op_##N}, \
    {GR_METHOD_MUL,                 (gr_funcptr) _nfx_mul_op_##N}, \
    {0,                             (gr_funcptr) NULL} \
};

DEF_NFX_INPUT(1)
DEF_NFX_INPUT(2)
DEF_NFX_INPUT(3)
DEF_NFX_INPUT(4)

static gr_method_tab_input _nfx_methods_input_g[] =
{
    NFX_COMMON_METHODS,
    {GR_METHOD_ADD,                 (gr_funcptr) _nfx_add_op},
    {GR_METHOD_SUB,                 (gr_funcptr) _nfx_sub_op},
    {GR_METHOD_MUL,                 (gr_funcptr) _nfx_mul_op},
    {GR_METHOD_SQR,                 (gr_funcptr) _nfx_sqr_op},
    {0,                             (gr_funcptr) NULL}
};

#define NFXC_COMMON_METHODS \
    {GR_METHOD_CTX_WRITE,           (gr_funcptr) _nfxc_ctx_write}, \
    {GR_METHOD_CTX_IS_THREADSAFE,   (gr_funcptr) _nfx_ctx_predicate_true}, \
    {GR_METHOD_CTX_GET_REAL_PREC,   (gr_funcptr) _nfx_ctx_get_real_prec}, \
    {GR_METHOD_INIT,                (gr_funcptr) _nfx_init}, \
    {GR_METHOD_CLEAR,               (gr_funcptr) _nfx_clear}, \
    {GR_METHOD_SWAP,                (gr_funcptr) _nfx_swap}, \
    {GR_METHOD_RANDTEST,            (gr_funcptr) _nfxc_randtest}, \
    {GR_METHOD_WRITE,               (gr_funcptr) _nfxc_write}, \
    {GR_METHOD_ZERO,                (gr_funcptr) _nfx_zero}, \
    {GR_METHOD_ONE,                 (gr_funcptr) _nfxc_one}, \
    {GR_METHOD_NEG_ONE,             (gr_funcptr) _nfxc_neg_one}, \
    {GR_METHOD_IS_ZERO,             (gr_funcptr) _nfxc_is_zero_op}, \
    {GR_METHOD_IS_ONE,              (gr_funcptr) _nfx_false_predicate}, \
    {GR_METHOD_IS_NEG_ONE,          (gr_funcptr) _nfx_false_predicate}, \
    {GR_METHOD_EQUAL,               (gr_funcptr) _nfxc_equal}, \
    {GR_METHOD_SET,                 (gr_funcptr) _nfx_set}, \
    {GR_METHOD_SET_SI,              (gr_funcptr) _nfxc_set_si}, \
    {GR_METHOD_SET_UI,              (gr_funcptr) _nfxc_set_ui}, \
    {GR_METHOD_DIV_UI,              (gr_funcptr) _nfxc_div_ui_op}, \
    {GR_METHOD_NEG,                 (gr_funcptr) _nfxc_neg}

#define DEF_NFXC_INPUT(N) \
static gr_method_tab_input _nfxc_methods_input_##N[] = \
{ \
    NFXC_COMMON_METHODS, \
    {GR_METHOD_ADD,                 (gr_funcptr) _nfxc_add_op_##N}, \
    {GR_METHOD_SUB,                 (gr_funcptr) _nfxc_sub_op_##N}, \
    {GR_METHOD_MUL,                 (gr_funcptr) _nfxc_mul_op_##N}, \
    {0,                             (gr_funcptr) NULL} \
};

DEF_NFXC_INPUT(1)
DEF_NFXC_INPUT(2)
DEF_NFXC_INPUT(3)
DEF_NFXC_INPUT(4)

static gr_method_tab_input _nfxc_methods_input_g[] =
{
    NFXC_COMMON_METHODS,
    {GR_METHOD_ADD,                 (gr_funcptr) _nfxc_add_op},
    {GR_METHOD_SUB,                 (gr_funcptr) _nfxc_sub_op},
    {GR_METHOD_MUL,                 (gr_funcptr) _nfxc_mul_generic},
    {GR_METHOD_SQR,                 (gr_funcptr) _nfxc_sqr_generic},
    {0,                             (gr_funcptr) NULL}
};

static void
_nfx_init_method_tables(void)
{
    if (!_nfx_methods_initialized)
    {
        gr_method_tab_init(_nfx_methods_tab[0], _nfx_methods_input_1);
        gr_method_tab_init(_nfx_methods_tab[1], _nfx_methods_input_2);
        gr_method_tab_init(_nfx_methods_tab[2], _nfx_methods_input_3);
        gr_method_tab_init(_nfx_methods_tab[3], _nfx_methods_input_4);
        gr_method_tab_init(_nfx_methods_tab[4], _nfx_methods_input_g);
        gr_method_tab_init(_nfxc_methods_tab[0], _nfxc_methods_input_1);
        gr_method_tab_init(_nfxc_methods_tab[1], _nfxc_methods_input_2);
        gr_method_tab_init(_nfxc_methods_tab[2], _nfxc_methods_input_3);
        gr_method_tab_init(_nfxc_methods_tab[3], _nfxc_methods_input_4);
        gr_method_tab_init(_nfxc_methods_tab[4], _nfxc_methods_input_g);
        _nfx_methods_initialized = 1;
    }
}

static slong
_nfx_tab_index(slong nlimbs)
{
    return (nlimbs <= 4) ? nlimbs - 1 : 4;
}

int
gr_dft_ctx_init_nfixed(gr_ctx_t ctx, slong nlimbs)
{
    if (nlimbs < 1 || nlimbs > GR_DFT_NFIXED_MAX_NLIMBS)
        return GR_UNABLE;

    ctx->which_ring = GR_CTX_UNKNOWN_DOMAIN;
    ctx->sizeof_elem = sizeof(ulong) * (nlimbs + 1);
    ctx->size_limit = WORD_MAX;
    NFIXED_CTX_NLIMBS(ctx) = nlimbs;

    _nfx_init_method_tables();
    ctx->methods = _nfx_methods_tab[_nfx_tab_index(nlimbs)];

    return GR_SUCCESS;
}

int
gr_dft_ctx_init_nfixed_complex(gr_ctx_t ctx, slong nlimbs)
{
    if (nlimbs < 1 || nlimbs > GR_DFT_NFIXED_MAX_NLIMBS)
        return GR_UNABLE;

    ctx->which_ring = GR_CTX_UNKNOWN_DOMAIN;
    ctx->sizeof_elem = 2 * sizeof(ulong) * (nlimbs + 1);
    ctx->size_limit = WORD_MAX;
    NFIXED_CTX_NLIMBS(ctx) = nlimbs;

    _nfx_init_method_tables();
    ctx->methods = _nfxc_methods_tab[_nfx_tab_index(nlimbs)];

    return GR_SUCCESS;
}

int
_gr_dft_ctx_is_nfixed(gr_ctx_t ctx)
{
    slong i;
    if (!_nfx_methods_initialized)
        return 0;
    for (i = 0; i < NFX_TAB_COUNT; i++)
        if (ctx->methods == _nfx_methods_tab[i])
            return 1;
    return 0;
}

int
_gr_dft_ctx_is_nfixed_complex(gr_ctx_t ctx)
{
    slong i;
    if (!_nfx_methods_initialized)
        return 0;
    for (i = 0; i < NFX_TAB_COUNT; i++)
        if (ctx->methods == _nfxc_methods_tab[i])
            return 1;
    return 0;
}

/* Root table *****************************************************/

/* Error bound (in ulp, as a complex modulus and hence valid for both
   the real and imaginary parts) for the entries of a root table of
   length n built by the doubling products of _gr_dft_nfixed_roots,
   where one complex multiplication contributes a modulus rounding of
   at most kmul ulp. The recurrence E(1) = 2,
   E(j) = E(floor(j/2)) + E(ceil(j/2)) + kmul has the exact solution
   bound E(j) <= (2 + kmul) j - kmul, established by induction; the
   factor below absorbs the fudge multipliers of the tracked
   recurrence over up to 64 doubling levels. */
double
_gr_dft_nfixed_root_err_bound(ulong n, double kmul)
{
    if (n <= 1)
        return 0.0;
    return ((2.0 + kmul) * (double) (n - 1) - kmul) * (1.0 + 1e-7) + 2.0;
}

/* Fill roots[0], ..., roots[n-1] with fixed-point approximations of
   w^j, w = exp(-2 pi i / n). The primitive root w is computed from
   arb with elevated internal precision so that its rounding (toward
   zero) to fixed point has error below 2 ulp; entries up to n/8 (or
   n/4 or n/2, depending on the divisibility of n) are computed as
   w^j = w^floor(j/2) w^ceil(j/2) (squarings for even j), giving the
   error recurrence E(j) <= E(floor(j/2)) + E(ceil(j/2)) + K, hence
   E(j) = O(j); the remaining entries are derived from the eight-fold
   symmetry of the unit circle by sign flips and swaps of the parts,
   which are exact. The maximum error over the table (a bound on
   |computed - exact| for both the real and imaginary parts, measured
   in ulp) is written to err_ulps. */
int
_gr_dft_nfixed_roots(gr_ptr roots, ulong n, double * err_ulps, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong nlimbs = NFIXED_CTX_NLIMBS(ctx);
    slong sz = ctx->sizeof_elem;
    slong prec = nlimbs * FLINT_BITS + 64;
    arb_t s, c;
    fmpq_t q;
    ulong j, L;
    /* rounding of one complex multiplication in this context,
       converted from a per-component to a modulus bound */
    const double KMUL = _gr_dft_nfixed_cmul_err_ulps(ctx) * 1.4142135623730951;

    status |= _nfxc_one(roots, ctx);
    *err_ulps = 0.0;

    if (n == 1)
        return status;

    arb_init(s);
    arb_init(c);
    fmpq_init(q);

    /* w = cos(-2 pi / n) + sin(-2 pi / n) i */
    fmpq_set_si(q, -2, n);
    arb_sin_cos_pi_fmpq(s, c, q, prec);

    {
        gr_ptr w1 = GR_ENTRY(roots, 1, sz);
        status |= gr_dft_nfixed_set_arf((gr_ptr) w1,
                arb_midref(c), ctx);   /* real part; real layout prefix */
        status |= gr_dft_nfixed_set_arf(
                (gr_ptr) NFXC_IM(w1, nlimbs), arb_midref(s), ctx);
    }

    arb_clear(s);
    arb_clear(c);
    fmpq_clear(q);

    /* entries computed by products: up to an eighth of the circle
       when 4 | n, a quarter when 2 | n, half (rounded down) for odd
       n; the remaining entries follow by symmetry */
    if (n % 4 == 0)
        L = (n + 7) / 8;
    else if (n % 2 == 0)
        L = (n + 3) / 4;
    else
        L = n / 2;

    for (j = 2; j <= L; j++)
    {
        ulong a = j / 2, b = j - a;

        if (a == b)
            status |= gr_sqr(GR_ENTRY(roots, j, sz),
                    GR_ENTRY(roots, a, sz), ctx);
        else
            status |= gr_mul(GR_ENTRY(roots, j, sz),
                    GR_ENTRY(roots, a, sz), GR_ENTRY(roots, b, sz), ctx);
    }

    /* symmetry fills (exact): writing w^j = x + y i,
       w^(n/4 - j) = -y - x i, w^(n/2 - j) = -x + y i,
       w^(n - j) = x - y i */
    if (n % 4 == 0)
    {
        for (j = L + 1; j <= n / 4; j++)
        {
            nn_srcptr src = (nn_srcptr) GR_ENTRY(roots, n / 4 - j, sz);
            nn_ptr dst = (nn_ptr) GR_ENTRY(roots, j, sz);
            slong i;

            for (i = 0; i <= nlimbs; i++)
            {
                dst[i] = src[nlimbs + 1 + i];             /* re = -im */
                dst[nlimbs + 1 + i] = src[i];             /* im = -re */
            }
            dst[0] ^= 1;
            dst[nlimbs + 1] ^= 1;
        }
    }
    if (n % 2 == 0)
    {
        for (j = FLINT_MAX(L, n / 4) + 1; j <= n / 2; j++)
        {
            nn_srcptr src = (nn_srcptr) GR_ENTRY(roots, n / 2 - j, sz);
            nn_ptr dst = (nn_ptr) GR_ENTRY(roots, j, sz);
            slong i;

            for (i = 0; i <= nlimbs; i++)
            {
                dst[i] = src[i];                          /* re = -re */
                dst[nlimbs + 1 + i] = src[nlimbs + 1 + i];
            }
            dst[0] ^= 1;
        }
    }
    for (j = n / 2 + 1; j < n; j++)
    {
        nn_srcptr src = (nn_srcptr) GR_ENTRY(roots, n - j, sz);
        nn_ptr dst = (nn_ptr) GR_ENTRY(roots, j, sz);
        slong i;

        for (i = 0; i <= nlimbs; i++)
        {
            dst[i] = src[i];                              /* conjugate */
            dst[nlimbs + 1 + i] = src[nlimbs + 1 + i];
        }
        dst[nlimbs + 1] ^= 1;
    }

    *err_ulps = _gr_dft_nfixed_root_err_bound(L + 1, KMUL);

    return status;
}
