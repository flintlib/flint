#ifndef MP_REAL_HAND_MULHI_INC
#define MP_REAL_HAND_MULHI_INC

/* Inline high-multiply helpers shared by the fixed-point series of the
   dense-table kernels (*_opt_<n>.c) and the table-driven kernels of
   small.h.

   All operands are fractions with the binary point above the top limb.
   The products keep only the top limbs: the low partials fall below the
   precision that is retained, and the error budgets of the callers
   (thousands of ulp against outputs of hundreds of bits) absorb this
   many times over. */

#include "ulong_extras.h"

/* Top 2 limbs of (a1:a0) * (b1:b0), discarding the low words of both
   cross products and the a0*b0 partial entirely: three multiplies and
   two adds.  This mirrors the helper in ulong_extras/ll_is_prime.c,
   which is file-local there and would be worth promoting to a header. */
FLINT_FORCE_INLINE void
_mp_real_mulhi_2x2_sloppy(ulong * r1, ulong * r0,
    ulong a1, ulong a0, ulong b1, ulong b0)
{
    ulong s2, u2, v3, v2, r3, r2;

    s2 = n_mulhi(a0, b1);
    u2 = n_mulhi(a1, b0);
    add_ssaaaa(r3, r2, UWORD(0), s2, UWORD(0), u2);
    umul_ppmm(v3, v2, a1, b1);
    add_ssaaaa(r3, r2, r3, r2, v3, v2);
    *r0 = r2;
    *r1 = r3;
}

/* Top 2 limbs of (a1:a0) * b with b a ONE-limb fraction: the b0 column
   of the product above vanishes, leaving two multiplies. */
FLINT_FORCE_INLINE void
_mp_real_mulhi_2x1(ulong * r1, ulong * r0, ulong a1, ulong a0, ulong b)
{
    ulong rh, rl, ch, cl;

    umul_ppmm(rh, rl, a1, b);
    umul_ppmm(ch, cl, a0, b);
    add_ssaaaa(rh, rl, rh, rl, UWORD(0), ch);
    (void) cl;
    *r0 = rl;
    *r1 = rh;
}

/* Top 2 limbs of (a1:a0)^2: the two cross products coincide, so this
   takes two multiplies rather than three -- measurably better than
   routing the squaring through the generic 2x2 form above. */
FLINT_FORCE_INLINE void
_mp_real_sqrhi_2x2(ulong * r1, ulong * r0, ulong a1, ulong a0)
{
    ulong rh, rl, ch, cl;

    umul_ppmm(rh, rl, a1, a1);
    umul_ppmm(ch, cl, a1, a0);
    add_ssaaaa(rh, rl, rh, rl, UWORD(0), ch);
    add_ssaaaa(rh, rl, rh, rl, UWORD(0), ch);
    (void) cl;
    *r0 = rl;
    *r1 = rh;
}

/* Top 3 limbs of a (3 limbs) times a 2-limb fraction b, which embeds
   as (0, b0, b1): the whole b-limb-0 column of the product vanishes,
   so three full products and two high-only ones suffice, against the
   six of a general 3x3. */
FLINT_FORCE_INLINE void
_mp_real_mulhi_3x2(ulong * r2, ulong * r1, ulong * r0,
    ulong a2, ulong a1, ulong a0, ulong b1, ulong b0)
{
    ulong t2, t1, x2, x1, x0, h1, h0, ph, pl, qh, ql;

    h0 = n_mulhi(a0, b1);
    add_ssaaaa(h1, h0, UWORD(0), h0, UWORD(0), n_mulhi(a1, b0));
    umul_ppmm(ph, pl, a1, b1);
    umul_ppmm(qh, ql, a2, b0);
    add_sssaaaaaa(x2, x1, x0, UWORD(0), ph, pl, UWORD(0), qh, ql);
    add_sssaaaaaa(x2, x1, x0, x2, x1, x0, UWORD(0), h1, h0);
    umul_ppmm(t2, t1, a2, b1);
    add_ssaaaa(*r2, *r1, t2, t1, x2, x1);
    *r0 = x0;
}

#endif /* MP_REAL_HAND_MULHI_INC */
