#ifndef FIXED_HAND_MULHI_INC
#define FIXED_HAND_MULHI_INC

/* Inline high-multiply helpers shared by the hand-written fixed-point
   series (exp_rs_opt_hand.inc, trig_rs_opt_hand.inc).

   All operands are fractions with the binary point above the top limb.
   The products keep only the top limbs: the low partials fall below the
   precision that is retained, and the error budgets of the callers
   (thousands of ulp against outputs of hundreds of bits) absorb this
   many times over. */

/* n_mulhi is newer than FLINT 3.3; define it when building against an
   older release so that this file stands alone.  (This block can go
   once the module only targets releases that ship it.) */
#if __FLINT_RELEASE < __FLINT_RELEASE_NUM(3, 4, 0)
static inline ulong
n_mulhi(ulong a, ulong b)
{
#if FLINT_BITS == 64 && defined(__GNUC__)
    return ((__uint128_t) a * (__uint128_t) b) >> 64;
#else
    ulong hi, lo;
    umul_ppmm(hi, lo, a, b);
    return hi;
#endif
}
#endif

/* Top 2 limbs of (a1:a0) * (b1:b0), discarding the low words of both
   cross products and the a0*b0 partial entirely: three multiplies and
   two adds.  This mirrors the helper in ulong_extras/ll_is_prime.c,
   which is file-local there and would be worth promoting to a header. */
FLINT_FORCE_INLINE void
_fixed_mulhi_2x2_sloppy(ulong * r1, ulong * r0,
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
_fixed_mulhi_2x1(ulong * r1, ulong * r0, ulong a1, ulong a0, ulong b)
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
_fixed_sqrhi_2x2(ulong * r1, ulong * r0, ulong a1, ulong a0)
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
_fixed_mulhi_3x2(ulong * r2, ulong * r1, ulong * r0,
    ulong a2, ulong a1, ulong a0, ulong b1, ulong b0)
{
    ulong t2, t1, t0, ph, pl;

    umul_ppmm(t2, t1, a2, b1);          /* weight 4 */
    t0 = 0;
    umul_ppmm(ph, pl, a1, b1);          /* weight 3 */
    add_sssaaaaaa(t2, t1, t0, t2, t1, t0, UWORD(0), ph, pl);
    umul_ppmm(ph, pl, a2, b0);          /* weight 3 */
    add_sssaaaaaa(t2, t1, t0, t2, t1, t0, UWORD(0), ph, pl);
    ph = n_mulhi(a0, b1);               /* weight 2: high word only */
    add_sssaaaaaa(t2, t1, t0, t2, t1, t0, UWORD(0), UWORD(0), ph);
    ph = n_mulhi(a1, b0);               /* weight 2: high word only */
    add_sssaaaaaa(t2, t1, t0, t2, t1, t0, UWORD(0), UWORD(0), ph);

    *r0 = t0;
    *r1 = t1;
    *r2 = t2;
}

#endif /* FIXED_HAND_MULHI_INC */
