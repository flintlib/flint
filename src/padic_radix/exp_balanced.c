/*
    Copyright (C) 2012 Sebastian Pancratz
    Copyright (C) 2012, 2014, 2022, 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "fmpz.h"
#include "arb/impl.h"
#include "padic_radix.h"
#include "gr.h"

/*
    Port of padic_exp_balanced (src/padic/exp_balanced.c) to the padic_radix
    module, with the multiprecision arithmetic carried out in radix_integer_t.

    The input is split additively at the digit boundaries p, p^2, p^4, p^8, ...
    (blocks covering the digit ranges [0,2), [2,4), [4,8), ...) and exp is the
    product of the pieces,

        exp(x) = prod_i exp(r_i),   x = sum_i r_i  (mod p^N).

    The exact cut sequence is a tuning factor; only the block extraction would
    change to retune it.

    Each exp(r_i) is evaluated by binary splitting of its Taylor series.  A block
    at digit valuation w has only n = O((N - w)/w) surviving terms, and the
    series evaluation exploits that valuation directly: writing r = p^w ru with
    ru a unit, binary splitting runs on ru while the geometric factor p^w is
    tracked as a digit shift, so a block's working precision is only (N - w) plus
    a small guard.

    The binary splitting uses the power-table scheme of _arb_exp_sum_bs_powtab:
    only (T, Q) are threaded through the tree, while the powers ru^step needed by
    the merges (step = (b-a)/2) are precomputed once into a table.  This avoids
    recomputing repeated powers and never forms the unused top-level power.
*/

/* res = digits [a, b) of x, shifted down to position 0
   (= floor(|x| / p^a) mod p^{b-a}), as a nonnegative residue.  Requires
   0 <= a <= b.  Used to pull a block out of the argument without first
   materialising its low zero digits. */
static void
_radix_integer_slice_digits(radix_integer_t res, const radix_integer_t x,
    slong a, slong b, const radix_t radix)
{
    radix_integer_mod_digits(res, x, b, radix);
    if (a > 0)
        radix_integer_rshift_digits(res, res, a, radix);
}

/* For 2 <= b - a < this bound, and when x = p^vr ru fits in a single limb, the
   binary-splitting leaf is evaluated by an iterative mpn accumulation (below)
   rather than by recursion.  Tuning parameter. */
#define PADIC_RADIX_EXP_BSPLIT_BASECASE 64
#define PADIC_RADIX_EXP_BC_LIMBS (2 * PADIC_RADIX_EXP_BSPLIT_BASECASE + 4)

/* Digit width of the smallest block, i.e. the last (lowest-valuation) series
   evaluated in the reverse-order accumulation.  Must be a power of two. */
#define PADIC_RADIX_EXP_SMALLEST_BLOCK 4

/*
    Iterative basecase for x = p^vr ru that fits in a single limb (so xl = x and
    ru_l = ru are single limbs).  Computes the unit-factored (Tu, Q) for the
    interval [a, b) by accumulating in plain mpn limbs and converting to radix
    once at the end -- the same technique as the e binary splitting in
    radix/profile/p-constants.c.

    With U = sum_{j=0}^{b-a-1} (b-1)!/(a+j)! * x^j, one has Tu = ru * U, and U, Q
    satisfy the single-limb-multiplier recurrence (Xp tracks x^{J-a})

        U <- U*J + Xp,   Q <- Q*J,   Xp <- Xp*x   (J = a+1, ..., b-1),

    started from U = 1, Q = a, Xp = x.  Every multiplier (J, x, ru) is one limb,
    so mpn_mul_1 / mpn_add suffice.
*/
static void
_bsplit_basecase_mpn(radix_integer_t Tu, radix_integer_t Q, slong a, slong b,
    ulong xl, ulong ru_l, slong l, const radix_t radix)
{
    ulong U[PADIC_RADIX_EXP_BC_LIMBS];
    ulong QQ[PADIC_RADIX_EXP_BC_LIMBS];
    ulong Xp[PADIC_RADIX_EXP_BC_LIMBS];
    slong Un, Qn, Xpn, J, need;
    ulong hi, cy;
    nn_ptr d;

    U[0] = 1;       Un = 1;
    QQ[0] = (ulong) a;  Qn = 1;
    Xp[0] = xl;     Xpn = 1;

    for (J = a + 1; J < b; J++)
    {
        /* U = U*J + Xp */
        hi = mpn_mul_1(U, U, Un, (ulong) J);
        U[Un] = hi; Un += (hi != 0);
        if (Un >= Xpn)
            cy = mpn_add(U, U, Un, Xp, Xpn);
        else
        {
            cy = mpn_add(U, Xp, Xpn, U, Un);
            Un = Xpn;
        }
        U[Un] = cy; Un += (cy != 0);

        /* Q = Q*J */
        hi = mpn_mul_1(QQ, QQ, Qn, (ulong) J);
        QQ[Qn] = hi; Qn += (hi != 0);

        /* Xp = Xp*x  (skip the final, unused update) */
        if (J + 1 < b)
        {
            hi = mpn_mul_1(Xp, Xp, Xpn, xl);
            Xp[Xpn] = hi; Xpn += (hi != 0);
        }
    }

    /* Tu = ru * U */
    hi = mpn_mul_1(U, U, Un, ru_l);
    U[Un] = hi; Un += (hi != 0);

    need = radix_set_mpn_need_alloc(Un, radix);
    d = radix_integer_fit_limbs(Tu, need, radix);
    Tu->size = radix_set_mpn(d, U, Un, radix);
    radix_integer_mod_limbs(Tu, Tu, l, radix);

    need = radix_set_mpn_need_alloc(Qn, radix);
    d = radix_integer_fit_limbs(Q, need, radix);
    Q->size = radix_set_mpn(d, QQ, Qn, radix);
    radix_integer_mod_limbs(Q, Q, l, radix);
}

/*
    Binary splitting of sum_{i=a}^{b-1} x^i / i! with x = p^vr ru, returning the
    unit-factored (Tu, Q): Q = (b-1)!/(a-1)! and Tu = T / p^vr, where T is the
    upstream numerator (T has valuation >= vr).  All reduced modulo B^l.  The
    merge

        Tu = Tu_L Q_R + p^{vr*step} (ru^step) Tu_R,   step = (b-a)/2,

    reads ru^step from the precomputed table xpow; the second term is skipped
    once the shift reaches the working precision.  P is never threaded.
*/
static void
_bsplit(radix_integer_t Tu, radix_integer_t Q, slong a, slong b,
    const slong * xexp, const radix_integer_struct * xpow, slong vr, slong l,
    ulong xl, ulong ru_l, int xl_fits, const radix_t radix)
{
    slong e = radix->exp;

    if (xl_fits && (b - a) < PADIC_RADIX_EXP_BSPLIT_BASECASE)
    {
        /* 1 <= b - a < tuning, and x fits in a limb: iterative mpn leaf.  This
           also covers b - a == 1 (the loop is then empty, giving Tu = ru), so
           the power table is never consulted along this path. */
        _bsplit_basecase_mpn(Tu, Q, a, b, xl, ru_l, l, radix);
    }
    else if (b - a == 1)
    {
        radix_integer_set(Tu, xpow + 0, radix);         /* ru^1 */
        radix_integer_set_ui(Q, (ulong) a, radix);
    }
    else if (b - a == 2)
    {
        slong i2 = _arb_get_exp_pos(xexp, 2);
        radix_integer_t tmp;
        radix_integer_init(tmp, radix);

        /* Tu = ru*(a+1) + p^vr * ru^2 */
        radix_integer_set_ui(tmp, (ulong) (a + 1), radix);
        radix_integer_mullow_limbs(Tu, xpow + 0, tmp, l, radix);
        radix_integer_lshift_digits(tmp, xpow + i2, vr, radix);
        radix_integer_mod_limbs(tmp, tmp, l, radix);
        radix_integer_add(Tu, Tu, tmp, radix);
        radix_integer_mod_limbs(Tu, Tu, l, radix);

        /* Q = a(a+1) */
        radix_integer_set_ui(Q, (ulong) a, radix);
        radix_integer_set_ui(tmp, (ulong) (a + 1), radix);
        radix_integer_mullow_limbs(Q, Q, tmp, l, radix);

        radix_integer_clear(tmp, radix);
    }
    else
    {
        slong step = (b - a) / 2;
        slong m = a + step;
        slong el = e * l;
        slong shift;
        int include;
        radix_integer_t TuR, QR, t2;

        radix_integer_init(TuR, radix);
        radix_integer_init(QR, radix);
        radix_integer_init(t2, radix);

        _bsplit(Tu, Q, a, m, xexp, xpow, vr, l, xl, ru_l, xl_fits, radix);
        _bsplit(TuR, QR, m, b, xexp, xpow, vr, l, xl, ru_l, xl_fits, radix);

        /* Tu = Tu*QR + p^{vr*step} * ru^step * TuR */
        radix_integer_mullow_limbs(Tu, Tu, QR, l, radix);

        if (step != 0 && vr > el / step)
            include = 0;
        else
        {
            shift = vr * step;
            include = (shift < el);
        }

        if (include)
        {
            slong is = _arb_get_exp_pos(xexp, step);
            radix_integer_mullow_limbs(t2, xpow + is, TuR, l, radix);
            radix_integer_lshift_digits(t2, t2, shift, radix);
            radix_integer_mod_limbs(t2, t2, l, radix);
            radix_integer_add(Tu, Tu, t2, radix);
            radix_integer_mod_limbs(Tu, Tu, l, radix);
        }

        radix_integer_mullow_limbs(Q, Q, QR, l, radix);

        radix_integer_clear(TuR, radix);
        radix_integer_clear(QR, radix);
        radix_integer_clear(t2, radix);
    }
}

/* exp(p^voff * s) = pnum / pden modulo p^N, where s is the slice of block
   digits already shifted down to position 0 and voff is the digit offset of the
   block within the argument (so the block value is p^voff * s).  Both pnum and
   pden are returned as units (the common power of p stripped).  Requires s != 0
   and that exp converges at the block.  Caller divides the accumulated products
   once at the end. */
static void
_bsplit_block_numden(radix_integer_t pnum, radix_integer_t pden,
    const radix_integer_t s, slong voff, slong N, const radix_t radix)
{
    ulong p = DIGIT_RADIX(radix);
    slong e = radix->exp;
    slong vs, vr, n, k, l, w, length, i;
    slong * xexp;
    radix_integer_struct * xpow;
    radix_integer_t ru, Q, Tu;
    ulong xl = 0, ru_l = 0;
    int xl_fits, use_table;

    vs = radix_integer_valuation_digits(s, radix);    /* valuation within slice */
    vr = voff + vs;                                   /* true block valuation */
    if (vr >= N)
    {
        radix_integer_one(pnum, radix);
        radix_integer_one(pden, radix);
        return;
    }

    n = _padic_radix_exp_bound(vr, N, p);
    if (n == 1)
    {
        radix_integer_one(pnum, radix);
        radix_integer_one(pden, radix);
        return;
    }

    k = (n >= 2) ? (n - 2) / (slong) (p - 1) : 0;     /* guard >= v_p((n-1)!) */
    /* Deferred division: the numerator and denominator are needed to the full
       target precision (not just the N - vr relevant for a divided block), so
       work to N + k digits. */
    l = (N + k + e - 1) / e;

    radix_integer_init(ru, radix);
    radix_integer_init(Q, radix);
    radix_integer_init(Tu, radix);

    radix_integer_rshift_digits(ru, s, vs, radix);    /* unit part of the block */
    radix_integer_mod_limbs(ru, ru, l, radix);

    /* The block value p^vr * ru fits in a single radix limb exactly when
       vr < e and the unit is below p^{e-vr}; then both the block value xl and
       the unit ru_l are single limbs and the iterative mpn basecase can be used
       for the small splitting leaves. */
    ru_l = (ru->size == 0) ? 0 : ru->d[0];
    xl_fits = (FLINT_ABS(ru->size) <= 1) && (vr < e) && (ru_l < radix->bpow[e - vr]);
    if (xl_fits)
        xl = radix->bpow[vr] * ru_l;                  /* block value, < B */
    else
        ru_l = 0;

    /* When the block fits in a limb and the whole series [1, n) is shorter than
       the basecase threshold, _bsplit takes the iterative leaf immediately and
       never consults the power table, so we skip building it entirely. */
    use_table = !(xl_fits && (n - 1) < PADIC_RADIX_EXP_BSPLIT_BASECASE);

    xexp = NULL;
    xpow = NULL;
    length = 0;

    if (use_table)
    {
        /* Power table ru^xexp[i] for the binary splitting over [1, n) (width n-1). */
        xexp = flint_calloc(2 * FLINT_BITS, sizeof(slong));
        length = _arb_compute_bs_exponents(xexp, n - 1);
        xpow = flint_malloc(sizeof(radix_integer_struct) * length);
        for (i = 0; i < length; i++)
            radix_integer_init(xpow + i, radix);

        radix_integer_set(xpow + 0, ru, radix);           /* ru^1 */
        for (i = 1; i < length; i++)
        {
            if (xexp[i] == 2 * xexp[i - 1])
                radix_integer_mullow_limbs(xpow + i, xpow + i - 1, xpow + i - 1, l, radix);
            else if (xexp[i] == 2 * xexp[i - 2])
                radix_integer_mullow_limbs(xpow + i, xpow + i - 2, xpow + i - 2, l, radix);
            else if (xexp[i] == 2 * xexp[i - 1] + 1)
            {
                radix_integer_mullow_limbs(xpow + i, xpow + i - 1, xpow + i - 1, l, radix);
                radix_integer_mullow_limbs(xpow + i, xpow + i, ru, l, radix);
            }
            else if (xexp[i] == 2 * xexp[i - 2] + 1)
            {
                radix_integer_mullow_limbs(xpow + i, xpow + i - 2, xpow + i - 2, l, radix);
                radix_integer_mullow_limbs(xpow + i, xpow + i, ru, l, radix);
            }
            else
                flint_throw(FLINT_ERROR, "padic_radix exp: power table malformed\n");
        }
    }

    _bsplit(Tu, Q, 1, n, xexp, xpow, vr, l, xl, ru_l, xl_fits, radix);

    /* Strip the common power of p (= v_p(Q) = v_p(Tu)) so that the numerator
       Q + p^vr Tu and denominator Q are returned as units. */
    w = radix_integer_valuation_digits(Tu, radix);
    if (w > 0)
    {
        radix_integer_rshift_digits(Tu, Tu, w, radix);
        radix_integer_rshift_digits(Q, Q, radix_integer_valuation_digits(Q, radix), radix);
    }

    /* pden = Q (unit), pnum = Q + p^vr Tu (unit), both mod p^N */
    radix_integer_mod_digits(pden, Q, N, radix);
    radix_integer_lshift_digits(pnum, Tu, vr, radix);
    radix_integer_add(pnum, pnum, Q, radix);
    radix_integer_mod_digits(pnum, pnum, N, radix);

    if (use_table)
    {
        for (i = 0; i < length; i++)
            radix_integer_clear(xpow + i, radix);
        flint_free(xpow);
        flint_free(xexp);
    }
    radix_integer_clear(ru, radix);
    radix_integer_clear(Q, radix);
    radix_integer_clear(Tu, radix);
}

/*
    rop = unit of exp(p^v u) modulo p^N, with u treated as an exact integer.

    The result is a 1-unit (valuation 0); rop is its nonnegative residue modulo
    p^N.  Requires v >= 1 (odd p) or v >= 2 (p = 2) and N >= 1.  Always returns
    GR_SUCCESS (an internal invariant failure aborts via flint_throw).
*/
int
_padic_radix_exp_balanced(radix_integer_t rop, const radix_integer_t u,
    slong v, slong N, const radix_t radix)
{
    const slong S = PADIC_RADIX_EXP_SMALLEST_BLOCK;
    ulong p = DIGIT_RADIX(radix);
    slong e = radix->exp;
    slong lN, D, jmax, j, tl, td;
    ulong top;
    radix_integer_t t, s, Pnum, Pden, pnum, pden;

    if (N <= 0)
    {
        radix_integer_zero(rop, radix);
        return GR_SUCCESS;
    }

    if (v >= N || radix_integer_is_zero(u, radix))
    {
        radix_integer_one(rop, radix);
        return GR_SUCCESS;
    }

    lN = (N + e - 1) / e;

    radix_integer_init(t, radix);

    /* t = p^v u  (mod p^N) */
    radix_integer_lshift_digits(t, u, v, radix);
    radix_integer_mod_digits(t, t, N, radix);

    if (t->size == 0)            /* exp of a value with valuation >= N is 1 */
    {
        radix_integer_one(rop, radix);
        radix_integer_clear(t, radix);
        return GR_SUCCESS;
    }

    /* Digit length D of t.  Blocks are [0,S), [S,2S), [2S,4S), [4S,8S), ...
       (smallest width S, then doubling), so they need to cover [0, S*2^{jmax-1})
       with S*2^{jmax-1} >= D. */
    tl = FLINT_ABS(t->size);
    top = t->d[tl - 1];
    td = 0;
    while (top != 0)
    {
        top /= p;
        td++;
    }
    D = (tl - 1) * e + td;

    jmax = 1;
    while ((S << (jmax - 1)) < D)
        jmax++;

    radix_integer_init(s, radix);
    radix_integer_init(Pnum, radix);
    radix_integer_init(Pden, radix);
    radix_integer_init(pnum, radix);
    radix_integer_init(pden, radix);
    radix_integer_one(Pnum, radix);
    radix_integer_one(Pden, radix);

    /* exp(x) = prod_j exp(block_j) = (prod_j pnum_j) / (prod_j pden_j).  Instead
       of dividing each block by its factorial, accumulate the unit numerator
       and denominator products and divide once.  Blocks are extracted on the
       fly highest valuation first (j = jmax down to 1), each shifted down to
       position 0 with the low zeros already removed; processing from the
       smallest contributions up keeps the running product growing downward. */
    for (j = jmax; j >= 1; j--)
    {
        slong a = (j == 1) ? 0 : (S << (j - 2));
        slong b = S << (j - 1);

        _radix_integer_slice_digits(s, t, a, b, radix);
        if (s->size != 0)
        {
            _bsplit_block_numden(pnum, pden, s, a, N, radix);
            radix_integer_mullow_limbs(Pnum, Pnum, pnum, lN, radix);
            radix_integer_mullow_limbs(Pden, Pden, pden, lN, radix);
        }
    }

    radix_integer_mod_digits(Pnum, Pnum, N, radix);
    radix_integer_mod_digits(Pden, Pden, N, radix);

    /* rop = Pnum / Pden  (mod p^N); both are units, so no stripping is needed. */
    if (!radix_integer_divmod_limbs(rop, Pnum, Pden, lN, radix))
        flint_throw(FLINT_ERROR, "_padic_radix_exp_balanced: denominator "
            "product is not invertible\n");

    radix_integer_mod_digits(rop, rop, N, radix);

    radix_integer_clear(t, radix);
    radix_integer_clear(s, radix);
    radix_integer_clear(Pnum, radix);
    radix_integer_clear(Pden, radix);
    radix_integer_clear(pnum, radix);
    radix_integer_clear(pden, radix);

    return GR_SUCCESS;
}

