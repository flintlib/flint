/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"
#include "fixed.h"
#include "arf.h"

/*
    Newton iteration backends for inversion, division and (reciprocal)
    square roots with the relaxed rounding modes ARF_RND_FAST (error at
    most 1 ulp in any direction) and ARF_RND_ACCURATE (at most 0.51 ulp),
    where an ulp is 2^(e - prec) with e the exponent of the exact result f
    (2^(e-1) <= |f| < 2^e).

    The mantissas are fed to the fixed_*_newton routines with n fraction
    limbs. Their results have absolute error at most 8 B^-n on values of
    magnitude at least 1/2, i.e. relative error below 2^(4 - 64 n); with
    64 n >= prec + g this is at most 2^-(prec + g - 4) <= ulp / 2^(g-4).
    The approximation Q is then rounded to nearest at prec bits, which adds
    at most half an ulp at the scale of Q. If Q has the same exponent as f,
    the total is (1/2 + 2^-(g-4)) ulp; if Q has straddled a power of two
    upwards, Q lies within ulp / 2^(g-4) above 2^e and rounds to exactly
    2^e (the nearest coarser-scale value), giving error below 2 ulp/2^(g-4);
    if Q has straddled downwards, the finer rounding contributes 1/4 ulp
    (and Q, within ulp/2^(g-4) below 2^(e-1), rounds up to 2^(e-1)). Hence
    g = 8 guard bits give at most 0.5625 ulp (ARF_RND_FAST) and g = 11 give
    at most 0.508 ulp (ARF_RND_ACCURATE). In all cases the returned value
    has the exponent of f, or is exactly 2^e, so callers adding one ulp of
    the result as an error bound (as arb does) remain correct.

    The cutoffs (in bits) were tuned with arb/profile/p-approx.c against
    arf_div / arf_ui_div / arf_sqrt / arf_rsqrt on full-precision operands.
    Division by a number much shorter than the precision is left to the
    standard algorithm, which is then cheaper.
*/

#ifndef ARF_INV_APPROX_NEWTON_CUTOFF
#define ARF_INV_APPROX_NEWTON_CUTOFF 8000
#endif
#ifndef ARF_DIV_APPROX_NEWTON_CUTOFF
#define ARF_DIV_APPROX_NEWTON_CUTOFF 12000
#endif
#ifndef ARF_SQRT_APPROX_NEWTON_CUTOFF
#define ARF_SQRT_APPROX_NEWTON_CUTOFF 100000
#endif
#ifndef ARF_RSQRT_APPROX_NEWTON_CUTOFF
#define ARF_RSQRT_APPROX_NEWTON_CUTOFF 30000
#endif

#define NEWTON_GUARD_BITS(rnd) ((rnd) == ARF_RND_ACCURATE ? 11 : 8)
#define NEWTON_LIMBS(prec, rnd) (((prec) + NEWTON_GUARD_BITS(rnd) + FLINT_BITS - 1) / FLINT_BITS)

/* set res from the n + 2 limb fixed-point value Q (n fraction limbs),
   times 2^exp, rounded to nearest at prec bits; Q may lie inside the
   buffer of res (_arf_set_round_mpn rounds in place) */
static int
_arf_set_fixed(arf_t res, mp_ptr Q, slong n, const fmpz_t exp, int sgnbit, slong prec)
{
    slong qn = n + 2, fix;
    int inexact;

    while (Q[qn - 1] == 0)
        qn--;

    inexact = _arf_set_round_mpn(res, &fix, Q, qn, sgnbit, prec, ARF_RND_NEAR);
    _fmpz_add_fast(ARF_EXPREF(res), exp, fix + (qn - n) * FLINT_BITS);
    return inexact;
}

int
_arf_want_newton_inv(const arf_t x, slong prec)
{
    return prec >= ARF_INV_APPROX_NEWTON_CUTOFF && !arf_is_special(x)
        && arf_bits(x) > prec / 2;
}

int
_arf_want_newton_div(const arf_t x, const arf_t y, slong prec)
{
    return prec >= ARF_DIV_APPROX_NEWTON_CUTOFF && !arf_is_special(x)
        && !arf_is_special(y) && arf_bits(y) > prec / 2;
}

int
_arf_want_newton_sqrt(const arf_t x, slong prec)
{
    return prec >= ARF_SQRT_APPROX_NEWTON_CUTOFF && !arf_is_special(x) && arf_sgn(x) > 0;
}

int
_arf_want_newton_rsqrt(const arf_t x, slong prec)
{
    return prec >= ARF_RSQRT_APPROX_NEWTON_CUTOFF && !arf_is_special(x) && arf_sgn(x) > 0;
}

/* Space for the n + 2 limb fixed-point result: when res does not alias an
   input, the limbs are computed directly in the buffer of res (sized to
   n + 2 limbs here, normalised afterwards by _arf_set_round_mpn), avoiding
   a temporary allocation; otherwise a temporary is used. */
static mp_ptr
_arf_approx_scratch(arf_t res, slong len, int aliased)
{
    mp_ptr Q;

    if (aliased)
        return flint_malloc(len * sizeof(mp_limb_t));

    ARF_GET_MPN_WRITE(Q, len, res);
    return Q;
}

int
_arf_inv_newton(arf_t res, const arf_t x, slong prec, arf_rnd_t rnd)
{
    slong xn, n;
    nn_srcptr xptr;
    mp_ptr Q;
    fmpz_t e;
    int inexact;

    ARF_GET_MPN_READONLY(xptr, xn, x);
    n = NEWTON_LIMBS(prec, rnd);

    fmpz_init(e);
    fmpz_neg(e, ARF_EXPREF(x));
    Q = _arf_approx_scratch(res, n + 2, res == x);
    fixed_inv_newton(Q, xptr, xn, n);
    inexact = _arf_set_fixed(res, Q, n, e, ARF_SGNBIT(x), prec);
    if (res == x)
        flint_free(Q);
    fmpz_clear(e);
    return inexact;
}

int
_arf_div_newton(arf_t res, const arf_t x, const arf_t y, slong prec, arf_rnd_t rnd)
{
    slong xn, yn, n;
    nn_srcptr xptr, yptr;
    mp_ptr Q;
    fmpz_t e;
    int inexact;

    if (arf_is_one(x))
        return _arf_inv_newton(res, y, prec, rnd);

    ARF_GET_MPN_READONLY(xptr, xn, x);
    ARF_GET_MPN_READONLY(yptr, yn, y);
    n = NEWTON_LIMBS(prec, rnd);

    fmpz_init(e);
    fmpz_sub(e, ARF_EXPREF(x), ARF_EXPREF(y));
    Q = _arf_approx_scratch(res, n + 2, res == x || res == y);
    fixed_div_newton(Q, xptr, xn, yptr, yn, n);
    inexact = _arf_set_fixed(res, Q, n, e, ARF_SGNBIT(x) ^ ARF_SGNBIT(y), prec);
    if (res == x || res == y)
        flint_free(Q);
    fmpz_clear(e);
    return inexact;
}

/*
    Operand for the square roots: the mantissa of x rescaled to an even
    exponent, i.e. X (exponent even) or X/2 (exponent odd), truncated to its
    top n + 2 limbs (the discarded limbs perturb the root by less than
    B^-(n+1), negligible against the 8 B^-n budget of the iteration). With
    an even exponent this is just a pointer into x; with an odd exponent
    the top limbs are shifted into the scratch buffer A (n + 3 limbs).
    Returns the limb count and writes the halved exponent to e.
*/
static mp_size_t
_arf_sqrt_operand(mp_srcptr * Aptr, mp_ptr A, fmpz_t e, const arf_t x, slong n)
{
    slong xn, An;
    nn_srcptr xptr;

    ARF_GET_MPN_READONLY(xptr, xn, x);
    An = FLINT_MIN(xn, n + 2);
    xptr += xn - An;

    if (fmpz_is_odd(ARF_EXPREF(x)))
    {
        A[An] = mpn_lshift(A, xptr, An, FLINT_BITS - 1);
        *Aptr = A;
        fmpz_add_ui(e, ARF_EXPREF(x), 1);
        fmpz_tdiv_q_2exp(e, e, 1);
        return An + 1;
    }
    else
    {
        *Aptr = xptr;
        fmpz_tdiv_q_2exp(e, ARF_EXPREF(x), 1);
        return An;
    }
}

#define SQRT_NEWTON_BODY(newton, negate_exp) \
    slong An, n; \
    mp_ptr A, Q; \
    mp_srcptr Aptr; \
    fmpz_t e; \
    int inexact; \
    TMP_INIT; \
    n = NEWTON_LIMBS(prec, rnd); \
    TMP_START; \
    /* scratch for the shifted operand only when the exponent is odd */ \
    A = fmpz_is_odd(ARF_EXPREF(x)) ? TMP_ALLOC((n + 3) * sizeof(mp_limb_t)) : NULL; \
    fmpz_init(e); \
    An = _arf_sqrt_operand(&Aptr, A, e, x, n); \
    if (negate_exp) \
        fmpz_neg(e, e); \
    Q = _arf_approx_scratch(res, n + 2, res == x); \
    newton(Q, Aptr, An, n); \
    inexact = _arf_set_fixed(res, Q, n, e, 0, prec); \
    if (res == x) \
        flint_free(Q); \
    fmpz_clear(e); \
    TMP_END; \
    return inexact;

int
_arf_sqrt_newton(arf_t res, const arf_t x, slong prec, arf_rnd_t rnd)
{
    SQRT_NEWTON_BODY(fixed_sqrt_newton, 0)
}

int
_arf_rsqrt_newton(arf_t res, const arf_t x, slong prec, arf_rnd_t rnd)
{
    SQRT_NEWTON_BODY(fixed_rsqrt_newton, 1)
}
