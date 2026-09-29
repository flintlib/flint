/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef DECIMAL_H
#define DECIMAL_H

#ifdef DECIMAL_INLINES_C
#define DECIMAL_INLINE
#else
#define DECIMAL_INLINE static inline
#endif

#include "radix.h"
#include "fmpz.h"
#include "fmpq_types.h"
#include "arf_types.h"
#include "arb_types.h"
#include "acf_types.h"
#include "acb_types.h"
#include "gr_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
    Decimal floating-point numbers (decfloat), decimal radius bounds (decmag)
    and decimal balls (decball).

    Representation
    --------------

    A decfloat is a value

        x = (-1)^s * M * B^v,      B = 10^e,

    where M is a nonnegative integer stored as radix_integer limbs in the
    limb radix B and v is an fmpz exponent counted in limbs. Nonzero values
    are kept canonical by stripping zero limbs at both ends of the mantissa
    (M0 != 0 and M_{n-1} != 0). Precision is measured in decimal digits:
    a value has sd(x) = digits(M) - val_10(M_0) significant digits, and all
    rounding operations produce results with sd(x) <= prec.

    Choosing limb-aligned exponents (rather than digit-aligned ones as in arf)
    means that addition never needs a digit shift and that rounding after any
    operation costs O(1) (drop whole limbs, then clear up to e - 1 digits in one
    boundary limb). The price is that a prec-digit value may occupy one limb
    more than ceil(prec/e) since the significant digits need not be aligned
    to a limb boundary. Rounding semantics (round to prec significant decimal
    digits in a given rounding mode) do not depend on the internal alignment.

    Special values (zero, +inf, -inf, nan) are encoded by M = 0 with a code
    in the exponent field.
*/

/* Rounding modes */
#define DECIMAL_RND_DOWN        0   /* toward zero */
#define DECIMAL_RND_UP          1   /* away from zero */
#define DECIMAL_RND_FLOOR       2   /* toward -inf */
#define DECIMAL_RND_CEIL        3   /* toward +inf */
#define DECIMAL_RND_NEAR        4   /* nearest, ties to even */
#define DECIMAL_RND_NEAR_AWAY   5   /* nearest, ties away from zero */
#define DECIMAL_RND_NEAR_ZERO   6   /* nearest, ties toward zero */
#define DECIMAL_RND_NUM         7
#define DECIMAL_RND_MASK        15
/* Internal flag which can be or'ed into a rounding mode to skip the
   exponent limit checks (used for intermediate computations). */
#define DECIMAL_RND_NOLIMITS    64
/* Flip floor/ceil (for negation) preserving flags */
#define DECIMAL_RND_NEGATE(rnd) ((((rnd) & DECIMAL_RND_MASK) == DECIMAL_RND_FLOOR) ? (((rnd) & ~DECIMAL_RND_MASK) | DECIMAL_RND_CEIL) : \
                                 ((((rnd) & DECIMAL_RND_MASK) == DECIMAL_RND_CEIL) ? (((rnd) & ~DECIMAL_RND_MASK) | DECIMAL_RND_FLOOR) : (rnd)))

/* Special precision meaning "no rounding" (exact ring Z[1/10]) */
#define DECIMAL_PREC_EXACT      WORD_MAX

/* Largest finite precision (digits) we allow */
#define DECIMAL_PREC_MAX        (WORD_MAX / 64)

/* Conversions to integers/rationals and exact additions refuse to
   materialize more than this many digits (GR_UNABLE is returned). */
#define DECIMAL_CONV_DIGITS_LIMIT (WORD(1) << 28)

/* Context flags */
#define DECIMAL_ALLOW_INF           1   /* allow +/-inf (else GR_UNABLE / GR_DOMAIN) */
#define DECIMAL_ALLOW_NAN           2   /* allow nan (else GR_UNABLE / GR_DOMAIN) */
#define DECIMAL_ALLOW_UNDERFLOW     4   /* flush to zero below emin (else GR_UNABLE) */
#define DECIMAL_SLOPPY_RADIUS       8   /* balls: use cheap ulp bounds for rounding errors */
#define DECIMAL_WRITE_SCIENTIFIC    16  /* always print in scientific notation */

/* DECIMAL_WORD_DIGITS is the largest k with 10^k < 2^FLINT_BITS. Radius
   (decmag) precision limits: the mantissa must fit a single limb with room
   for a product of two mantissas times ten (10^(2 rp + 1) < 2^FLINT_BITS). */
#if FLINT_BITS == 64
#define DECIMAL_WORD_DIGITS 19
#define DECMAG_MAX_PREC   9
#else
#define DECIMAL_WORD_DIGITS 9
#define DECMAG_MAX_PREC   4
#endif
#define DECMAG_MIN_PREC   1

/* largest k with 5^k < 2^FLINT_BITS, and 5^k */
#if FLINT_BITS == 64
#define DECIMAL_POW5_MAX_EXP 27
#define DECIMAL_POW5_MAX UWORD(7450580596923828125)
#else
#define DECIMAL_POW5_MAX_EXP 13
#define DECIMAL_POW5_MAX UWORD(1220703125)
#endif
#define DECMAG_DEFAULT_PREC 4

/* Special value codes stored in the exponent when the mantissa is zero */
#define DECFLOAT_EXP_ZERO     0
#define DECFLOAT_EXP_POS_INF  1
#define DECFLOAT_EXP_NEG_INF  (-1)
#define DECFLOAT_EXP_NAN      2

typedef struct
{
    radix_integer_struct m;   /* signed mantissa (limbs in radix B = 10^e) */
    fmpz exp;                 /* exponent in limbs; special code if m == 0 */
}
decfloat_struct;

typedef decfloat_struct decfloat_t[1];
typedef decfloat_struct * decfloat_ptr;
typedef const decfloat_struct * decfloat_srcptr;

typedef struct
{
    ulong m;                  /* 0 = zero, UWORD_MAX = inf, else 1 <= m < 10^DECMAG_MAX_PREC
                                 (normally with rad_prec digits, see below) */
    fmpz exp;                 /* value = m * 10^exp (exp measured in digits) */
}
decmag_struct;

typedef decmag_struct decmag_t[1];
typedef decmag_struct * decmag_ptr;
typedef const decmag_struct * decmag_srcptr;

typedef struct
{
    decfloat_struct mid;
    decmag_struct rad;
}
decball_struct;

typedef decball_struct decball_t[1];
typedef decball_struct * decball_ptr;
typedef const decball_struct * decball_srcptr;

#define DECFLOAT_MANT(x)      (&((x)->m))
#define DECFLOAT_EXPREF(x)    (&((x)->exp))
#define DECFLOAT_SIZE(x)      ((x)->m.size)
#define DECFLOAT_IS_SPECIAL(x) ((x)->m.size == 0)
#define DECFLOAT_IS_ZERO(x)   ((x)->m.size == 0 && (x)->exp == DECFLOAT_EXP_ZERO)
#define DECFLOAT_IS_POS_INF(x) ((x)->m.size == 0 && (x)->exp == DECFLOAT_EXP_POS_INF)
#define DECFLOAT_IS_NEG_INF(x) ((x)->m.size == 0 && (x)->exp == DECFLOAT_EXP_NEG_INF)
#define DECFLOAT_IS_INF(x)    (DECFLOAT_IS_POS_INF(x) || DECFLOAT_IS_NEG_INF(x))
#define DECFLOAT_IS_NAN(x)    ((x)->m.size == 0 && (x)->exp == DECFLOAT_EXP_NAN)
#define DECFLOAT_IS_FINITE(x) ((x)->m.size != 0 || (x)->exp == DECFLOAT_EXP_ZERO)
#define DECFLOAT_SGNBIT(x)    ((x)->m.size < 0)

#define DECMAG_IS_ZERO(x)     ((x)->m == 0)
#define DECMAG_IS_INF(x)      ((x)->m == UWORD_MAX)
#define DECMAG_IS_SPECIAL(x)  ((x)->m == 0 || (x)->m == UWORD_MAX)
#define DECMAG_IS_FINITE(x)   ((x)->m != UWORD_MAX)
/* whether the mantissa is special or has exactly rad_prec digits */
#define DECMAG_IS_NORMALIZED(x, ctx) (DECMAG_IS_SPECIAL(x) || ((x)->m >= DECIMAL_CTX_RAD_POW1(ctx) && (x)->m < DECIMAL_CTX_RAD_POW(ctx)))

#define DECBALL_MIDREF(x)     (&((x)->mid))
#define DECBALL_RADREF(x)     (&((x)->rad))

typedef struct
{
    decfloat_struct re;
    decfloat_struct im;
}
deccfloat_struct;

typedef deccfloat_struct deccfloat_t[1];
typedef deccfloat_struct * deccfloat_ptr;
typedef const deccfloat_struct * deccfloat_srcptr;

typedef struct
{
    decball_struct re;
    decball_struct im;
}
deccball_struct;

typedef deccball_struct deccball_t[1];
typedef deccball_struct * deccball_ptr;
typedef const deccball_struct * deccball_srcptr;

#define DECCFLOAT_REALREF(z)  (&((z)->re))
#define DECCFLOAT_IMAGREF(z)  (&((z)->im))
#define DECCBALL_REALREF(z)   (&((z)->re))
#define DECCBALL_IMAGREF(z)   (&((z)->im))

/* Context */

typedef enum { DECIMAL_CTX_FLOAT = 0, DECIMAL_CTX_BALL = 1, DECIMAL_CTX_CFLOAT = 2, DECIMAL_CTX_CBALL = 3 } decimal_ctx_which;

typedef struct
{
    radix_struct radix;       /* b = 10, B = 10^e */
    slong prec;               /* default precision in digits, or DECIMAL_PREC_EXACT */
    int rnd;                  /* default rounding mode (real parts) */
    int rnd_im;               /* rounding mode for imaginary parts (complex floats) */
    int flags;
    slong rad_prec;           /* radius precision in digits, in [1, DECMAG_MAX_PREC] */
    ulong rad_pow;            /* 10^rad_prec */
    ulong rad_pow1;           /* 10^(rad_prec-1) */
    slong emin;               /* lower bound on the scientific exponent (WORD_MIN: none) */
    slong emax;               /* upper bound on the scientific exponent (WORD_MAX: none) */
    decimal_ctx_which which;
}
decimal_ctx_struct;

#define DECIMAL_CTX(ctx)        ((decimal_ctx_struct *) (GR_CTX_DATA_AS_PTR(ctx)))
#define DECIMAL_CTX_RADIX(ctx)  (&(DECIMAL_CTX(ctx)->radix))
#define DECIMAL_CTX_PREC(ctx)   (DECIMAL_CTX(ctx)->prec)
#define DECIMAL_CTX_RND(ctx)    (DECIMAL_CTX(ctx)->rnd)
#define DECIMAL_CTX_RND_IM(ctx) (DECIMAL_CTX(ctx)->rnd_im)
#define DECIMAL_CTX_WHICH(ctx)  (DECIMAL_CTX(ctx)->which)
#define DECIMAL_CTX_IS_BALL(ctx) (DECIMAL_CTX(ctx)->which & 1)
#define DECIMAL_CTX_IS_COMPLEX(ctx) (DECIMAL_CTX(ctx)->which & 2)
#define DECIMAL_CTX_FLAGS(ctx)  (DECIMAL_CTX(ctx)->flags)
#define DECIMAL_CTX_RAD_PREC(ctx) (DECIMAL_CTX(ctx)->rad_prec)
#define DECIMAL_CTX_RAD_POW(ctx) (DECIMAL_CTX(ctx)->rad_pow)
#define DECIMAL_CTX_RAD_POW1(ctx) (DECIMAL_CTX(ctx)->rad_pow1)
#define DECIMAL_CTX_EMIN(ctx)   (DECIMAL_CTX(ctx)->emin)
#define DECIMAL_CTX_EMAX(ctx)   (DECIMAL_CTX(ctx)->emax)
#define DECIMAL_CTX_E(ctx)      (DECIMAL_CTX(ctx)->radix.exp)
#define DECIMAL_CTX_B(ctx)      (DECIMAL_CTX(ctx)->radix.B.n)
#define DECIMAL_CTX_HAS_EXP_LIMITS(ctx) (DECIMAL_CTX(ctx)->emin != WORD_MIN || DECIMAL_CTX(ctx)->emax != WORD_MAX)
#define DECIMAL_CTX_IS_EXACT(ctx) (DECIMAL_CTX(ctx)->prec == DECIMAL_PREC_EXACT)
/* balls represent real numbers only: their midpoints are never infinite or
   nan (an infinite radius denotes the whole real line), so the flags have
   no effect in ball contexts. TODO: a flag to make a ball context represent
   extended real intervals. */
#define DECIMAL_CTX_ALLOW_INF(ctx) ((DECIMAL_CTX(ctx)->flags & DECIMAL_ALLOW_INF) && !DECIMAL_CTX_IS_BALL(ctx))
#define DECIMAL_CTX_ALLOW_NAN(ctx) ((DECIMAL_CTX(ctx)->flags & DECIMAL_ALLOW_NAN) && !DECIMAL_CTX_IS_BALL(ctx))
#define DECIMAL_CTX_INF_ON_OVERFLOW(ctx) DECIMAL_CTX_ALLOW_INF(ctx)
#define DECIMAL_CTX_ALLOW_UNDERFLOW(ctx) (DECIMAL_CTX(ctx)->flags & DECIMAL_ALLOW_UNDERFLOW)
#define DECIMAL_CTX_PRECISE_RADIUS(ctx) (!(DECIMAL_CTX(ctx)->flags & DECIMAL_SLOPPY_RADIUS))

void gr_ctx_init_decfloat(gr_ctx_t ctx, slong prec, int flags);
void gr_ctx_init_decball(gr_ctx_t ctx, slong prec, int flags);
void gr_ctx_init_deccfloat(gr_ctx_t ctx, slong prec, int flags);
void gr_ctx_init_deccball(gr_ctx_t ctx, slong prec, int flags);
void _gr_ctx_init_decimal(gr_ctx_t ctx, decimal_ctx_which which, unsigned int e, slong prec, int rnd, int flags);
void gr_ctx_init_decfloat_randtest(gr_ctx_t ctx, flint_rand_t state, slong max_prec);
void gr_ctx_init_decball_randtest(gr_ctx_t ctx, flint_rand_t state, slong max_prec);
void gr_ctx_init_deccfloat_randtest(gr_ctx_t ctx, flint_rand_t state, slong max_prec);
void gr_ctx_init_deccball_randtest(gr_ctx_t ctx, flint_rand_t state, slong max_prec);

void decimal_ctx_clear(gr_ctx_t ctx);
int decimal_ctx_write(gr_stream_t out, gr_ctx_t ctx);
truth_t _decimal_ctx_is_vector_space(gr_ctx_t ctx);
void decimal_ctx_set_prec(gr_ctx_t ctx, slong prec);
void decimal_ctx_set_rnd(gr_ctx_t ctx, int rnd);
void decimal_ctx_set_rnd_im(gr_ctx_t ctx, int rnd);
void decimal_ctx_set_rad_prec(gr_ctx_t ctx, slong rad_prec);
void decimal_ctx_set_exp_limits(gr_ctx_t ctx, slong emin, slong emax);
void decimal_ctx_set_flags(gr_ctx_t ctx, int flags);
slong decimal_ctx_get_prec(gr_ctx_t ctx);
int decimal_ctx_get_rnd(gr_ctx_t ctx);
int decimal_ctx_get_rnd_im(gr_ctx_t ctx);
slong decimal_ctx_get_rad_prec(gr_ctx_t ctx);
int decimal_ctx_get_flags(gr_ctx_t ctx);
slong decimal_ctx_get_limb_digits(gr_ctx_t ctx);
void decimal_ctx_get_exp_limits(slong * emin, slong * emax, gr_ctx_t ctx);
int _decimal_ctx_set_real_prec(gr_ctx_t ctx, slong prec_bits);
int _decimal_ctx_get_real_prec(slong * res, gr_ctx_t ctx);

/* Rounding information reported by the internal rounding routines. */
typedef struct
{
    int inexact;        /* whether rounding changed the value */
    int increased;      /* whether |result| > |exact value| */
    int underflow;      /* result flushed to zero (exponent below emin) */
    int overflow;       /* result replaced by an infinity (exponent above emax) */
}
decimal_rounding_info;

/* ------------------------------------------------------------------------- */
/* decfloat                                                                  */
/* ------------------------------------------------------------------------- */

DECIMAL_INLINE void
decfloat_init(decfloat_t res, gr_ctx_t ctx)
{
    res->m.d = NULL;
    res->m.alloc = 0;
    res->m.size = 0;
    res->exp = DECFLOAT_EXP_ZERO;
}

DECIMAL_INLINE void
decfloat_clear(decfloat_t res, gr_ctx_t ctx)
{
    if (res->m.d != NULL)
        flint_free(res->m.d);
    fmpz_clear(&res->exp);
}

DECIMAL_INLINE void
decfloat_swap(decfloat_t x, decfloat_t y, gr_ctx_t ctx)
{
    FLINT_SWAP(decfloat_struct, *x, *y);
}

DECIMAL_INLINE void
decfloat_set_shallow(decfloat_t res, const decfloat_t x, gr_ctx_t ctx)
{
    *res = *x;
}

DECIMAL_INLINE int
decfloat_zero(decfloat_t res, gr_ctx_t ctx)
{
    res->m.size = 0;
    fmpz_zero(&res->exp);
    return GR_SUCCESS;
}

int decfloat_pos_inf(decfloat_t res, gr_ctx_t ctx);
int decfloat_neg_inf(decfloat_t res, gr_ctx_t ctx);
int decfloat_nan(decfloat_t res, gr_ctx_t ctx);
int decfloat_one(decfloat_t res, gr_ctx_t ctx);
int decfloat_neg_one(decfloat_t res, gr_ctx_t ctx);

/* Unconditionally create special values (ignoring context flags). */
void _decfloat_pos_inf(decfloat_t res);
void _decfloat_neg_inf(decfloat_t res);
void _decfloat_nan(decfloat_t res);

/* Whether x is admitted as an operand: GR_SUCCESS if x is finite or a
   special value allowed by the context, otherwise GR_DOMAIN (infinity)
   or GR_UNABLE (nan). */
DECIMAL_INLINE int
_decfloat_operand_status(const decfloat_t x, gr_ctx_t ctx)
{
    if (DECFLOAT_IS_FINITE(x))
        return GR_SUCCESS;
    if (DECFLOAT_IS_NAN(x))
        return DECIMAL_CTX_ALLOW_NAN(ctx) ? GR_SUCCESS : GR_UNABLE;
    return DECIMAL_CTX_ALLOW_INF(ctx) ? GR_SUCCESS : GR_DOMAIN;
}

#define DECFLOAT_CHECK_OPERAND(x, ctx) \
    do { int _status = _decfloat_operand_status(x, ctx); if (_status != GR_SUCCESS) return _status; } while (0)

truth_t decfloat_is_zero(const decfloat_t x, gr_ctx_t ctx);
truth_t decfloat_is_one(const decfloat_t x, gr_ctx_t ctx);
truth_t decfloat_is_neg_one(const decfloat_t x, gr_ctx_t ctx);
truth_t decfloat_equal(const decfloat_t x, const decfloat_t y, gr_ctx_t ctx);

DECIMAL_INLINE int _decfloat_is_special(const decfloat_t x) { return DECFLOAT_IS_SPECIAL(x); }
DECIMAL_INLINE int _decfloat_is_finite(const decfloat_t x) { return DECFLOAT_IS_FINITE(x); }
DECIMAL_INLINE int _decfloat_is_nan(const decfloat_t x) { return DECFLOAT_IS_NAN(x); }
DECIMAL_INLINE int _decfloat_is_inf(const decfloat_t x) { return DECFLOAT_IS_INF(x); }
DECIMAL_INLINE int _decfloat_is_pos_inf(const decfloat_t x) { return DECFLOAT_IS_POS_INF(x); }
DECIMAL_INLINE int _decfloat_is_neg_inf(const decfloat_t x) { return DECFLOAT_IS_NEG_INF(x); }

/* Whether x is an integer (special values are not integers, except zero). */
DECIMAL_INLINE int
_decfloat_is_int(const decfloat_t x, gr_ctx_t ctx)
{
    if (DECFLOAT_IS_SPECIAL(x))
        return DECFLOAT_IS_ZERO(x);
    return fmpz_sgn(&x->exp) >= 0;
}

truth_t decfloat_is_integer(const decfloat_t x, gr_ctx_t ctx);

/* Comparisons following the gr conventions (GR_UNABLE for nan), and raw
   versions returning the result directly (0 for comparisons with nan). */
int decfloat_cmp(int * res, const decfloat_t x, const decfloat_t y, gr_ctx_t ctx);
int decfloat_cmpabs(int * res, const decfloat_t x, const decfloat_t y, gr_ctx_t ctx);
int decfloat_sgn(decfloat_t res, const decfloat_t x, gr_ctx_t ctx);
int _decfloat_sgn(const decfloat_t x, gr_ctx_t ctx);
int _decfloat_cmp(const decfloat_t x, const decfloat_t y, gr_ctx_t ctx);
int _decfloat_cmpabs(const decfloat_t x, const decfloat_t y, gr_ctx_t ctx);
int _decfloat_cmp_ui(const decfloat_t x, ulong c, gr_ctx_t ctx);
int _decfloat_cmp_si(const decfloat_t x, slong c, gr_ctx_t ctx);
int _decfloat_cmpabs_ui(const decfloat_t x, ulong c, gr_ctx_t ctx);
int _decfloat_cmpabs_finite(const decfloat_t x, const decfloat_t y);

/* Number of digits of M where |x| = M 10^v with M not divisible by 10,
   and the number of limbs of the mantissa (0 for special values). */
slong decfloat_digits(const decfloat_t x, gr_ctx_t ctx);
slong decfloat_limbs(const decfloat_t x, gr_ctx_t ctx);
/* Number of digits of the limb mantissa (including trailing zeros of the
   lowest limb), and the exponent v above clamped to +/- DECFLOAT_EXP_CLAMP. */
slong _decfloat_mant_digits(const decfloat_t x, gr_ctx_t ctx);
slong _decfloat_val10_clamped(const decfloat_t x, gr_ctx_t ctx);

/* Digits and limbs of |x| at absolute positions 10^k and B^k. */
ulong decfloat_get_digit_si(const decfloat_t x, slong k, gr_ctx_t ctx);
ulong decfloat_get_limb_si(const decfloat_t x, slong k, gr_ctx_t ctx);
int decfloat_set_digit_si(decfloat_t res, const decfloat_t x, slong k, ulong c, gr_ctx_t ctx);
int decfloat_set_limb_si(decfloat_t res, const decfloat_t x, slong k, ulong c, gr_ctx_t ctx);

/* Conversions to and from integers in the limb radix of the context. */
int decfloat_get_radix_integer(radix_integer_t res, const decfloat_t x, gr_ctx_t ctx);
int decfloat_get_radix_integer_Bexp_fmpz(radix_integer_t m, fmpz_t exp, const decfloat_t x, gr_ctx_t ctx);
int decfloat_set_radix_integer(decfloat_t res, const radix_integer_t m, gr_ctx_t ctx);
int decfloat_set_radix_integer_Bexp_fmpz(decfloat_t res, const radix_integer_t m, const fmpz_t exp, gr_ctx_t ctx);
/* Scientific exponent E such that 10^E <= |x| < 10^(E+1). Requires x nonzero finite. */
void decfloat_get_sci_exp(fmpz_t E, const decfloat_t x, gr_ctx_t ctx);
/* Whether the scientific exponent fits in a slong; if so, returns it. */
int decfloat_get_sci_exp_si(slong * E, const decfloat_t x, gr_ctx_t ctx);

/* Rounding */
int decfloat_set_round(decfloat_t res, const decfloat_t x, slong prec, int rnd, gr_ctx_t ctx);
int decfloat_set_round_info(decfloat_t res, const decfloat_t x, slong prec, int rnd, decimal_rounding_info * info, decmag_ptr err, gr_ctx_t ctx);

/* Assignment (rounded to the context precision) */
int decfloat_set(decfloat_t res, const decfloat_t x, gr_ctx_t ctx);
int decfloat_set_si(decfloat_t res, slong x, gr_ctx_t ctx);
int decfloat_set_ui(decfloat_t res, ulong x, gr_ctx_t ctx);
int decfloat_set_fmpz(decfloat_t res, const fmpz_t x, gr_ctx_t ctx);
int decfloat_set_fmpq(decfloat_t res, const fmpq_t x, gr_ctx_t ctx);
int decfloat_set_d(decfloat_t res, double x, gr_ctx_t ctx);
int decfloat_set_str(decfloat_t res, const char * s, gr_ctx_t ctx);
int decfloat_set_other(decfloat_t res, gr_srcptr x, gr_ctx_t x_ctx, gr_ctx_t ctx);
int _decfloat_set_decfloat_other(decfloat_t res, const decfloat_t y, gr_ctx_t x_ctx, gr_ctx_t ctx);
int decfloat_set_fmpz_10exp_fmpz(decfloat_t res, const fmpz_t m, const fmpz_t e, gr_ctx_t ctx);
int decfloat_set_fmpz_10exp_si(decfloat_t res, const fmpz_t m, slong e, gr_ctx_t ctx);
int decfloat_set_si_10exp_si(decfloat_t res, slong m, slong e, gr_ctx_t ctx);
int decfloat_set_fmpz_2exp_fmpz(decfloat_t res, const fmpz_t m, const fmpz_t e, gr_ctx_t ctx);
int decfloat_set_arf(decfloat_t res, const arf_t x, gr_ctx_t ctx);

/* Explicit-precision assignment */
int decfloat_set_round_si(decfloat_t res, slong x, slong prec, int rnd, gr_ctx_t ctx);
int decfloat_set_round_ui(decfloat_t res, ulong x, slong prec, int rnd, gr_ctx_t ctx);
int decfloat_set_round_fmpz(decfloat_t res, const fmpz_t x, slong prec, int rnd, gr_ctx_t ctx);
int decfloat_set_round_fmpq(decfloat_t res, const fmpq_t x, slong prec, int rnd, gr_ctx_t ctx);
int decfloat_set_round_fmpq_reference(decfloat_t res, const fmpq_t x, slong prec, int rnd, gr_ctx_t ctx);
int decfloat_set_round_d(decfloat_t res, double x, slong prec, int rnd, gr_ctx_t ctx);
/* exact conversion (no rounding, no exponent limits) of *y, of type
   DECIMAL_SCALAR_UI (ulong), _SI (slong), _FMPZ (fmpz) or _D (double) */
#define DECIMAL_SCALAR_UI   0
#define DECIMAL_SCALAR_SI   1
#define DECIMAL_SCALAR_FMPZ 2
#define DECIMAL_SCALAR_D    3
int _decfloat_set_scalar_exact(decfloat_t res, const void * y, int type, gr_ctx_t ctx);
int decfloat_set_round_str(decfloat_t res, const char * s, slong prec, int rnd, gr_ctx_t ctx);
int _decfloat_set_str_literal(decfloat_t res, const char * s, slong prec, int rnd, gr_ctx_t ctx);
int decfloat_set_round_fmpz_10exp_fmpz(decfloat_t res, const fmpz_t m, const fmpz_t e, slong prec, int rnd, gr_ctx_t ctx);
int decfloat_set_round_fmpz_2exp_fmpz(decfloat_t res, const fmpz_t m, const fmpz_t e, slong prec, int rnd, gr_ctx_t ctx);
int _decfloat_set_round_fmpz_2exp_err(decfloat_t res, const fmpz_t m, const fmpz_t t, slong prec, int rnd, decimal_rounding_info * info, decmag_ptr err, gr_ctx_t ctx);
int _decfloat_set_round_fmpz_2exp_scaled(decfloat_t res, const fmpz_t m, slong t, slong prec, int rnd, decimal_rounding_info * info, decmag_ptr err, gr_ctx_t ctx);
int _decfloat_set_exact_decimal_fmpz_2exp(decfloat_t res, const fmpz_t m, const fmpz_t t, slong prec, gr_ctx_t ctx);
int decfloat_set_round_arf(decfloat_t res, const arf_t x, slong prec, int rnd, gr_ctx_t ctx);

/* Conversion out */
int decfloat_get_fmpz(fmpz_t res, const decfloat_t x, gr_ctx_t ctx);
int decfloat_get_fmpq(fmpq_t res, const decfloat_t x, gr_ctx_t ctx);
int decfloat_get_si(slong * res, const decfloat_t x, gr_ctx_t ctx);
int decfloat_get_ui(ulong * res, const decfloat_t x, gr_ctx_t ctx);
int decfloat_get_d(double * res, const decfloat_t x, gr_ctx_t ctx);
int decfloat_get_fmpz_10exp_fmpz(fmpz_t m, fmpz_t e, const decfloat_t x, gr_ctx_t ctx);
int decfloat_get_arf(arf_t res, const decfloat_t x, slong prec_bits, int rnd, gr_ctx_t ctx);
int _decfloat_maybe_dyadic(const decfloat_t x, gr_ctx_t ctx);
int decfloat_get_arb(arb_t res, const decfloat_t x, slong prec_bits, gr_ctx_t ctx);
int decfloat_get_fmpz_fixed_si(fmpz_t res, const decfloat_t x, slong e, int rnd, gr_ctx_t ctx);

/* Strings */
char * decfloat_get_str(const decfloat_t x, gr_ctx_t ctx);
char * decfloat_get_str_sci(const decfloat_t x, gr_ctx_t ctx);
char * _decfloat_get_str_mode(const decfloat_t x, int mode, gr_ctx_t ctx);
slong _decimal_write_number_bound(slong L, const fmpz_t t);
slong _decimal_write_number(char * out, int negative, const char * digits, slong L, const fmpz_t t, int mode);
slong _decfloat_get_digits(char * buf, fmpz_t t, const decfloat_t x, gr_ctx_t ctx);
slong _decmag_get_digits(char * buf, fmpz_t t, const decmag_t x, gr_ctx_t ctx);
int decfloat_write(gr_stream_t out, const decfloat_t x, gr_ctx_t ctx);
int decfloat_write_sci(gr_stream_t out, const decfloat_t x, gr_ctx_t ctx);

int decfloat_randtest(decfloat_t res, flint_rand_t state, gr_ctx_t ctx);
int decfloat_randtest_special(decfloat_t res, flint_rand_t state, gr_ctx_t ctx);

/* Arithmetic (context precision and rounding mode) */
int decfloat_neg(decfloat_t res, const decfloat_t x, gr_ctx_t ctx);
int decfloat_abs(decfloat_t res, const decfloat_t x, gr_ctx_t ctx);
int decfloat_add(decfloat_t res, const decfloat_t x, const decfloat_t y, gr_ctx_t ctx);
int decfloat_sub(decfloat_t res, const decfloat_t x, const decfloat_t y, gr_ctx_t ctx);
int decfloat_mul(decfloat_t res, const decfloat_t x, const decfloat_t y, gr_ctx_t ctx);
int decfloat_sqr(decfloat_t res, const decfloat_t x, gr_ctx_t ctx);
int decfloat_div(decfloat_t res, const decfloat_t x, const decfloat_t y, gr_ctx_t ctx);
int decfloat_inv(decfloat_t res, const decfloat_t x, gr_ctx_t ctx);
int decfloat_sqrt(decfloat_t res, const decfloat_t x, gr_ctx_t ctx);
int decfloat_rsqrt(decfloat_t res, const decfloat_t x, gr_ctx_t ctx);
int decfloat_add_ui(decfloat_t res, const decfloat_t x, ulong y, gr_ctx_t ctx);
int decfloat_add_si(decfloat_t res, const decfloat_t x, slong y, gr_ctx_t ctx);
int decfloat_add_fmpz(decfloat_t res, const decfloat_t x, const fmpz_t y, gr_ctx_t ctx);
int decfloat_sub_ui(decfloat_t res, const decfloat_t x, ulong y, gr_ctx_t ctx);
int decfloat_sub_si(decfloat_t res, const decfloat_t x, slong y, gr_ctx_t ctx);
int decfloat_sub_fmpz(decfloat_t res, const decfloat_t x, const fmpz_t y, gr_ctx_t ctx);
int decfloat_mul_ui(decfloat_t res, const decfloat_t x, ulong y, gr_ctx_t ctx);
int decfloat_mul_si(decfloat_t res, const decfloat_t x, slong y, gr_ctx_t ctx);
int decfloat_mul_fmpz(decfloat_t res, const decfloat_t x, const fmpz_t y, gr_ctx_t ctx);
int decfloat_div_ui(decfloat_t res, const decfloat_t x, ulong y, gr_ctx_t ctx);
int decfloat_div_si(decfloat_t res, const decfloat_t x, slong y, gr_ctx_t ctx);
int decfloat_div_fmpz(decfloat_t res, const decfloat_t x, const fmpz_t y, gr_ctx_t ctx);
int decfloat_mul_two(decfloat_t res, const decfloat_t x, gr_ctx_t ctx);
int decfloat_mul_10exp_si(decfloat_t res, const decfloat_t x, slong e, gr_ctx_t ctx);
int decfloat_mul_10exp_fmpz(decfloat_t res, const decfloat_t x, const fmpz_t e, gr_ctx_t ctx);
int decfloat_mul_2exp_si(decfloat_t res, const decfloat_t x, slong e, gr_ctx_t ctx);
int decfloat_mul_2exp_fmpz(decfloat_t res, const decfloat_t x, const fmpz_t e, gr_ctx_t ctx);
int decfloat_floor(decfloat_t res, const decfloat_t x, gr_ctx_t ctx);
int decfloat_ceil(decfloat_t res, const decfloat_t x, gr_ctx_t ctx);
int decfloat_trunc(decfloat_t res, const decfloat_t x, gr_ctx_t ctx);
int decfloat_nint(decfloat_t res, const decfloat_t x, gr_ctx_t ctx);
int decfloat_pow_ui(decfloat_t res, const decfloat_t x, ulong y, gr_ctx_t ctx);
int decfloat_pow_si(decfloat_t res, const decfloat_t x, slong y, gr_ctx_t ctx);
int decfloat_pow_fmpz(decfloat_t res, const decfloat_t x, const fmpz_t y, gr_ctx_t ctx);
int _decfloat_pow_ziv(decfloat_t res, const decfloat_t x, const decfloat_t y, gr_ctx_t ctx);
int _decfloat_pow_int_exact(decfloat_t res, const decfloat_t x, const fmpz_t n, slong prec, int rnd, gr_ctx_t ctx);

/* Arithmetic with explicit precision and rounding mode. All functions
   return GR_SUCCESS or an error code; info and err (both optional, may be
   NULL) receive the rounding information and an upper bound for the
   absolute rounding error. */
int decfloat_neg_round(decfloat_t res, const decfloat_t x, slong prec, int rnd, gr_ctx_t ctx);
int decfloat_abs_round(decfloat_t res, const decfloat_t x, slong prec, int rnd, gr_ctx_t ctx);
int decfloat_add_round(decfloat_t res, const decfloat_t x, const decfloat_t y, slong prec, int rnd, gr_ctx_t ctx);
int decfloat_sub_round(decfloat_t res, const decfloat_t x, const decfloat_t y, slong prec, int rnd, gr_ctx_t ctx);
int decfloat_mul_round(decfloat_t res, const decfloat_t x, const decfloat_t y, slong prec, int rnd, gr_ctx_t ctx);
int decfloat_sqr_round(decfloat_t res, const decfloat_t x, slong prec, int rnd, gr_ctx_t ctx);
int decfloat_div_round(decfloat_t res, const decfloat_t x, const decfloat_t y, slong prec, int rnd, gr_ctx_t ctx);
int decfloat_inv_round(decfloat_t res, const decfloat_t x, slong prec, int rnd, gr_ctx_t ctx);
int decfloat_sqrt_round(decfloat_t res, const decfloat_t x, slong prec, int rnd, gr_ctx_t ctx);
int decfloat_rsqrt_round(decfloat_t res, const decfloat_t x, slong prec, int rnd, gr_ctx_t ctx);
int decfloat_mul_10exp_fmpz_round(decfloat_t res, const decfloat_t x, const fmpz_t e, slong prec, int rnd, gr_ctx_t ctx);
int decfloat_mul_10exp_si_round(decfloat_t res, const decfloat_t x, slong e, slong prec, int rnd, gr_ctx_t ctx);
int decfloat_mul_2exp_fmpz_round(decfloat_t res, const decfloat_t x, const fmpz_t e, slong prec, int rnd, gr_ctx_t ctx);
int decfloat_floor_round(decfloat_t res, const decfloat_t x, slong prec, int rnd, gr_ctx_t ctx);
int decfloat_ceil_round(decfloat_t res, const decfloat_t x, slong prec, int rnd, gr_ctx_t ctx);
int decfloat_trunc_round(decfloat_t res, const decfloat_t x, slong prec, int rnd, gr_ctx_t ctx);
int decfloat_nint_round(decfloat_t res, const decfloat_t x, slong prec, int rnd, gr_ctx_t ctx);

int _decfloat_add(decfloat_t res, const decfloat_t x, const decfloat_t y, int negate_y, slong prec, int rnd, decimal_rounding_info * info, decmag_ptr err, gr_ctx_t ctx);
int _decfloat_mul(decfloat_t res, const decfloat_t x, const decfloat_t y, slong prec, int rnd, decimal_rounding_info * info, decmag_ptr err, gr_ctx_t ctx);
int _decfloat_div(decfloat_t res, const decfloat_t x, const decfloat_t y, slong prec, int rnd, decimal_rounding_info * info, decmag_ptr err, gr_ctx_t ctx);
int _decfloat_sqrt(decfloat_t res, const decfloat_t x, slong prec, int rnd, decimal_rounding_info * info, decmag_ptr err, gr_ctx_t ctx);
int _decfloat_mul_10exp(decfloat_t res, const decfloat_t x, const fmpz_t e, slong prec, int rnd, decimal_rounding_info * info, decmag_ptr err, gr_ctx_t ctx);
int _decfloat_round_to_int(decfloat_t res, const decfloat_t x, int int_rnd, slong prec, int rnd, decimal_rounding_info * info, decmag_ptr err, gr_ctx_t ctx);

/* Low-level: set res to the value (-1)^negative * (d, n) * B^exp [+ eps]
   rounded to prec digits, where 0 <= eps < 1 with eps > 0 indicating an
   unrepresented positive tail (sticky). Limbs of d may be zero at either end.
   The buffer d may alias res->m.d. Exponent overflow/underflow are handled
   according to the context. */
int _decfloat_set_round_limbs(decfloat_t res, nn_srcptr d, slong n, int negative, const fmpz_t exp, int eps, slong prec, int rnd, decimal_rounding_info * info, decmag_ptr err, gr_ctx_t ctx);

/* Rounds a mantissa in place. Returns the number of low limbs to remove
   (the caller must shift the result down by this amount and adjust the
   exponent) and writes the new limb count to *newn. */
slong _decimal_round_mantissa(nn_ptr d, slong n, int negative, int eps, slong prec, int rnd, slong * newn, decimal_rounding_info * info, decmag_ptr err, const fmpz_t exp, gr_ctx_t ctx);

/* Checks exponent limits and special-value permissions after a computation. */
int _decfloat_finalize(decfloat_t res, gr_ctx_t ctx);
int _decfloat_finalize_info(decfloat_t res, decimal_rounding_info * info, gr_ctx_t ctx);

/* Vector operations */
int decfloat_vec_dot(decfloat_t res, const decfloat_t initial, int subtract, decfloat_srcptr vec1, decfloat_srcptr vec2, slong len, gr_ctx_t ctx);
int decfloat_vec_dot_rev(decfloat_t res, const decfloat_t initial, int subtract, decfloat_srcptr vec1, decfloat_srcptr vec2, slong len, gr_ctx_t ctx);

/* ------------------------------------------------------------------------- */
/* decmag                                                                    */
/* ------------------------------------------------------------------------- */

DECIMAL_INLINE void
_decmag_init(decmag_t res, gr_ctx_t ctx)
{
    res->m = 0;
    res->exp = 0;
}

DECIMAL_INLINE void
_decmag_clear(decmag_t res, gr_ctx_t ctx)
{
    fmpz_clear(&res->exp);
}

DECIMAL_INLINE void
_decmag_swap(decmag_t x, decmag_t y, gr_ctx_t ctx)
{
    FLINT_SWAP(decmag_struct, *x, *y);
}

DECIMAL_INLINE void
_decmag_zero(decmag_t res, gr_ctx_t ctx)
{
    res->m = 0;
    fmpz_zero(&res->exp);
}

DECIMAL_INLINE void
_decmag_inf(decmag_t res, gr_ctx_t ctx)
{
    res->m = UWORD_MAX;
    fmpz_zero(&res->exp);
}

DECIMAL_INLINE int _decmag_is_zero(const decmag_t x, gr_ctx_t ctx) { return DECMAG_IS_ZERO(x); }
DECIMAL_INLINE int _decmag_is_inf(const decmag_t x, gr_ctx_t ctx) { return DECMAG_IS_INF(x); }
DECIMAL_INLINE int _decmag_is_finite(const decmag_t x, gr_ctx_t ctx) { return DECMAG_IS_FINITE(x); }

void _decmag_set(decmag_t res, const decmag_t x, gr_ctx_t ctx);
void _decmag_set_lower(decmag_t res, const decmag_t x, gr_ctx_t ctx);
void _decmag_set_round(decmag_t res, const decmag_t x, slong prec, gr_ctx_t ctx);
void _decmag_one(decmag_t res, gr_ctx_t ctx);
/* Number of mantissa digits (0 for special values), and the scientific
   exponent E such that 10^E <= x < 10^(E+1) (x finite nonzero). */
slong _decmag_digits(const decmag_t x);
void _decmag_get_sci_exp(fmpz_t E, const decmag_t x);
int _decmag_equal(const decmag_t x, const decmag_t y, gr_ctx_t ctx);
int _decmag_cmp(const decmag_t x, const decmag_t y, gr_ctx_t ctx);
int _decmag_cmp_10exp_si(const decmag_t x, slong e, gr_ctx_t ctx);
int _decmag_is_10exp(const decmag_t x, gr_ctx_t ctx);

/* Normalize a mantissa m * 10^exp to rp digits, rounding up (or down). */
void _decmag_set_ui_10exp_fmpz(decmag_t res, ulong m, const fmpz_t exp, gr_ctx_t ctx);
void _decmag_set_ui_10exp_si(decmag_t res, ulong m, slong exp, gr_ctx_t ctx);
void _decmag_set_ui_10exp_fmpz_lower(decmag_t res, ulong m, const fmpz_t exp, gr_ctx_t ctx);
void _decmag_set_uiui_10exp_fmpz(decmag_t res, ulong hi, ulong lo, ulong hi_radix, int sticky, const fmpz_t exp, gr_ctx_t ctx);
void _decmag_set_ui(decmag_t res, ulong x, gr_ctx_t ctx);
void _decmag_set_ui_lower(decmag_t res, ulong x, gr_ctx_t ctx);
void _decmag_set_fmpz(decmag_t res, const fmpz_t x, gr_ctx_t ctx);
void _decmag_set_fmpz_lower(decmag_t res, const fmpz_t x, gr_ctx_t ctx);
int _decmag_set_d(decmag_t res, double x, gr_ctx_t ctx);
void _decmag_set_10exp_si(decmag_t res, slong e, gr_ctx_t ctx);
void _decmag_set_10exp_fmpz(decmag_t res, const fmpz_t e, gr_ctx_t ctx);

/* Upper (lower) bounds for |x|. */
void _decmag_set_decfloat(decmag_t res, const decfloat_t x, gr_ctx_t ctx);
void _decmag_set_decfloat_lower(decmag_t res, const decfloat_t x, gr_ctx_t ctx);
/* Exact conversion of a finite decmag (or inf if allowed). */
int _decmag_get_decfloat(decfloat_t res, const decmag_t x, gr_ctx_t ctx);
void _decmag_get_fmpq(fmpq_t res, const decmag_t x, gr_ctx_t ctx);
int _decmag_get_d(double * res, const decmag_t x, gr_ctx_t ctx);
void _decmag_get_mag(mag_t res, const decmag_t x, gr_ctx_t ctx);
void _decmag_set_mag(decmag_t res, const mag_t x, gr_ctx_t ctx);

/* ulp(x) at precision prec: 10^(E - prec + 1) where E is the scientific
   exponent of x. Requires x nonzero finite. */
void _decmag_set_ulp(decmag_t res, const decfloat_t x, slong prec, gr_ctx_t ctx);

/* Arithmetic with upper-bound rounding unless stated otherwise. */
void _decmag_add(decmag_t res, const decmag_t x, const decmag_t y, gr_ctx_t ctx);
void _decmag_add_lower(decmag_t res, const decmag_t x, const decmag_t y, gr_ctx_t ctx);
void _decmag_sub_lower(decmag_t res, const decmag_t x, const decmag_t y, gr_ctx_t ctx);
void _decmag_mul(decmag_t res, const decmag_t x, const decmag_t y, gr_ctx_t ctx);
void _decmag_mul_lower(decmag_t res, const decmag_t x, const decmag_t y, gr_ctx_t ctx);
void _decmag_addmul(decmag_t res, const decmag_t x, const decmag_t y, gr_ctx_t ctx);
void _decmag_div(decmag_t res, const decmag_t x, const decmag_t y, gr_ctx_t ctx);
void _decmag_div_lower(decmag_t res, const decmag_t x, const decmag_t y, gr_ctx_t ctx);
void _decmag_inv(decmag_t res, const decmag_t x, gr_ctx_t ctx);
void _decmag_sqrt(decmag_t res, const decmag_t x, gr_ctx_t ctx);
void _decmag_sqrt_lower(decmag_t res, const decmag_t x, gr_ctx_t ctx);
void _decmag_rsqrt(decmag_t res, const decmag_t x, gr_ctx_t ctx);
void _decmag_mul_ui(decmag_t res, const decmag_t x, ulong y, gr_ctx_t ctx);
void _decmag_div_ui(decmag_t res, const decmag_t x, ulong y, gr_ctx_t ctx);
void _decmag_mul_10exp_si(decmag_t res, const decmag_t x, slong e, gr_ctx_t ctx);
void _decmag_mul_10exp_fmpz(decmag_t res, const decmag_t x, const fmpz_t e, gr_ctx_t ctx);
void _decmag_max(decmag_t res, const decmag_t x, const decmag_t y, gr_ctx_t ctx);
void _decmag_min(decmag_t res, const decmag_t x, const decmag_t y, gr_ctx_t ctx);
void _decmag_pow_ui(decmag_t res, const decmag_t x, ulong e, gr_ctx_t ctx);
void _decmag_hypot(decmag_t res, const decmag_t x, const decmag_t y, gr_ctx_t ctx);

char * _decmag_get_str(const decmag_t x, gr_ctx_t ctx);
int _decmag_write(gr_stream_t out, const decmag_t x, gr_ctx_t ctx);
void _decmag_randtest(decmag_t res, flint_rand_t state, gr_ctx_t ctx);
void _decmag_randtest_special(decmag_t res, flint_rand_t state, gr_ctx_t ctx);

/* ------------------------------------------------------------------------- */
/* decball                                                                   */
/* ------------------------------------------------------------------------- */

DECIMAL_INLINE void
decball_init(decball_t res, gr_ctx_t ctx)
{
    decfloat_init(&res->mid, ctx);
    _decmag_init(&res->rad, ctx);
}

/* out of line: expanded at many call sites */
void decball_clear(decball_t res, gr_ctx_t ctx);

DECIMAL_INLINE void
decball_swap(decball_t x, decball_t y, gr_ctx_t ctx)
{
    FLINT_SWAP(decball_struct, *x, *y);
}

DECIMAL_INLINE void
decball_set_shallow(decball_t res, const decball_t x, gr_ctx_t ctx)
{
    *res = *x;
}

/* Accessors (function versions of DECBALL_MIDREF / DECBALL_RADREF) */
DECIMAL_INLINE decfloat_ptr decball_midref(decball_t x) { return &x->mid; }
DECIMAL_INLINE decmag_ptr decball_radref(decball_t x) { return &x->rad; }

DECIMAL_INLINE int _decball_is_exact(const decball_t x, gr_ctx_t ctx) { return DECMAG_IS_ZERO(&x->rad); }
DECIMAL_INLINE int _decball_is_finite(const decball_t x, gr_ctx_t ctx) { return DECFLOAT_IS_FINITE(&x->mid) && DECMAG_IS_FINITE(&x->rad); }

int decball_zero(decball_t res, gr_ctx_t ctx);
int decball_one(decball_t res, gr_ctx_t ctx);
int decball_neg_one(decball_t res, gr_ctx_t ctx);
int decball_zero_pm_inf(decball_t res, gr_ctx_t ctx);

int decball_set(decball_t res, const decball_t x, gr_ctx_t ctx);
int decball_set_round(decball_t res, const decball_t x, slong prec, gr_ctx_t ctx);
int decball_set_decfloat(decball_t res, const decfloat_t x, gr_ctx_t ctx);
int decball_set_si(decball_t res, slong x, gr_ctx_t ctx);
int decball_set_ui(decball_t res, ulong x, gr_ctx_t ctx);
int decball_set_fmpz(decball_t res, const fmpz_t x, gr_ctx_t ctx);
int decball_set_fmpq(decball_t res, const fmpq_t x, gr_ctx_t ctx);
int decball_set_d(decball_t res, double x, gr_ctx_t ctx);
int decball_set_str(decball_t res, const char * s, gr_ctx_t ctx);
int decball_set_other(decball_t res, gr_srcptr x, gr_ctx_t x_ctx, gr_ctx_t ctx);
int _decball_set_decfloat_other(decball_t res, const decfloat_t x, gr_ctx_t x_ctx, gr_ctx_t ctx);
int _decball_set_decball_other(decball_t res, const decball_t y, gr_ctx_t x_ctx, gr_ctx_t ctx);
int decball_set_fmpz_10exp_fmpz(decball_t res, const fmpz_t m, const fmpz_t e, gr_ctx_t ctx);
int decball_set_round2(decball_t res, const decball_t x, slong prec, slong rad_prec, gr_ctx_t ctx);
int decball_set_interval_mid_rad(decball_t res, const decball_t m, const decball_t r, gr_ctx_t ctx);
int decball_set_interval(decball_t res, const decball_t lo, const decball_t hi, gr_ctx_t ctx);
int decball_add_rad(decball_t res, const decball_t x, const decball_t r, gr_ctx_t ctx);
int decball_set_interval_mid_inf(decball_t res, const decball_t m, gr_ctx_t ctx);
int decball_mid(decball_t res, const decball_t x, gr_ctx_t ctx);
int decball_rad(decball_t res, const decball_t x, gr_ctx_t ctx);
int decball_shell(decball_t res, const decball_t x, gr_ctx_t ctx);
int decball_lower(decball_t res, const decball_t x, gr_ctx_t ctx);
int decball_upper(decball_t res, const decball_t x, gr_ctx_t ctx);
int decball_abs_lower(decball_t res, const decball_t x, gr_ctx_t ctx);
int decball_abs_upper(decball_t res, const decball_t x, gr_ctx_t ctx);
int decball_set_arb(decball_t res, const arb_t x, gr_ctx_t ctx);
int decball_get_arb(arb_t res, const decball_t x, slong prec_bits, gr_ctx_t ctx);
int decball_get_fmpz(fmpz_t res, const decball_t x, gr_ctx_t ctx);
int decball_get_fmpq(fmpq_t res, const decball_t x, gr_ctx_t ctx);
int decball_get_d(double * res, const decball_t x, gr_ctx_t ctx);
int decball_get_mid(decfloat_t res, const decball_t x, gr_ctx_t ctx);
void decball_get_rad(decmag_t res, const decball_t x, gr_ctx_t ctx);
void decball_get_abs_ubound(decmag_t res, const decball_t x, gr_ctx_t ctx);
void decball_get_abs_lbound(decmag_t res, const decball_t x, gr_ctx_t ctx);

int decball_randtest(decball_t res, flint_rand_t state, gr_ctx_t ctx);
char * decball_get_str(const decball_t x, gr_ctx_t ctx);
int decball_write(gr_stream_t out, const decball_t x, gr_ctx_t ctx);

truth_t decball_is_zero(const decball_t x, gr_ctx_t ctx);
truth_t decball_is_one(const decball_t x, gr_ctx_t ctx);
truth_t decball_is_neg_one(const decball_t x, gr_ctx_t ctx);
truth_t decball_equal(const decball_t x, const decball_t y, gr_ctx_t ctx);
truth_t decball_is_integer(const decball_t x, gr_ctx_t ctx);
int _decball_contains_zero(const decball_t x, gr_ctx_t ctx);
int _decball_contains_decfloat(const decball_t x, const decfloat_t y, gr_ctx_t ctx);
int _decball_contains_fmpq(const decball_t x, const fmpq_t y, gr_ctx_t ctx);
int _decball_contains_fmpz(const decball_t x, const fmpz_t y, gr_ctx_t ctx);
int _decball_contains_si(const decball_t x, slong y, gr_ctx_t ctx);
int _decball_contains(const decball_t x, const decball_t y, gr_ctx_t ctx);
int _decball_overlaps(const decball_t x, const decball_t y, gr_ctx_t ctx);
int _decball_is_positive(const decball_t x, gr_ctx_t ctx);
int _decball_is_negative(const decball_t x, gr_ctx_t ctx);
int _decball_is_nonnegative(const decball_t x, gr_ctx_t ctx);
int _decball_is_nonpositive(const decball_t x, gr_ctx_t ctx);
int _decball_contains_negative(const decball_t x, gr_ctx_t ctx);
int _decball_contains_positive(const decball_t x, gr_ctx_t ctx);
int _decball_contains_nonnegative(const decball_t x, gr_ctx_t ctx);
int _decball_contains_nonpositive(const decball_t x, gr_ctx_t ctx);
int decball_cmp(int * res, const decball_t x, const decball_t y, gr_ctx_t ctx);
int decball_cmpabs(int * res, const decball_t x, const decball_t y, gr_ctx_t ctx);
int decball_sgn(decball_t res, const decball_t x, gr_ctx_t ctx);

int decball_neg(decball_t res, const decball_t x, gr_ctx_t ctx);
int decball_abs(decball_t res, const decball_t x, gr_ctx_t ctx);
int decball_add(decball_t res, const decball_t x, const decball_t y, gr_ctx_t ctx);
int decball_sub(decball_t res, const decball_t x, const decball_t y, gr_ctx_t ctx);
int decball_mul(decball_t res, const decball_t x, const decball_t y, gr_ctx_t ctx);
int decball_sqr(decball_t res, const decball_t x, gr_ctx_t ctx);
int decball_div(decball_t res, const decball_t x, const decball_t y, gr_ctx_t ctx);
int decball_inv(decball_t res, const decball_t x, gr_ctx_t ctx);
int decball_sqrt(decball_t res, const decball_t x, gr_ctx_t ctx);
int decball_rsqrt(decball_t res, const decball_t x, gr_ctx_t ctx);
int decball_add_ui(decball_t res, const decball_t x, ulong y, gr_ctx_t ctx);
int decball_add_si(decball_t res, const decball_t x, slong y, gr_ctx_t ctx);
int decball_add_fmpz(decball_t res, const decball_t x, const fmpz_t y, gr_ctx_t ctx);
int decball_sub_ui(decball_t res, const decball_t x, ulong y, gr_ctx_t ctx);
int decball_sub_si(decball_t res, const decball_t x, slong y, gr_ctx_t ctx);
int decball_sub_fmpz(decball_t res, const decball_t x, const fmpz_t y, gr_ctx_t ctx);
int decball_mul_ui(decball_t res, const decball_t x, ulong y, gr_ctx_t ctx);
int decball_mul_si(decball_t res, const decball_t x, slong y, gr_ctx_t ctx);
int decball_mul_fmpz(decball_t res, const decball_t x, const fmpz_t y, gr_ctx_t ctx);
int decball_div_ui(decball_t res, const decball_t x, ulong y, gr_ctx_t ctx);
int decball_div_si(decball_t res, const decball_t x, slong y, gr_ctx_t ctx);
int decball_div_fmpz(decball_t res, const decball_t x, const fmpz_t y, gr_ctx_t ctx);
int decball_mul_two(decball_t res, const decball_t x, gr_ctx_t ctx);
int decball_mul_10exp_si(decball_t res, const decball_t x, slong e, gr_ctx_t ctx);
int decball_mul_10exp_fmpz(decball_t res, const decball_t x, const fmpz_t e, gr_ctx_t ctx);
int decball_mul_2exp_si(decball_t res, const decball_t x, slong e, gr_ctx_t ctx);
int decball_mul_2exp_fmpz(decball_t res, const decball_t x, const fmpz_t e, gr_ctx_t ctx);
int decball_floor(decball_t res, const decball_t x, gr_ctx_t ctx);
int decball_ceil(decball_t res, const decball_t x, gr_ctx_t ctx);
int decball_trunc(decball_t res, const decball_t x, gr_ctx_t ctx);
int decball_nint(decball_t res, const decball_t x, gr_ctx_t ctx);
int decball_pow_ui(decball_t res, const decball_t x, ulong y, gr_ctx_t ctx);
int decball_pow_si(decball_t res, const decball_t x, slong y, gr_ctx_t ctx);
int decball_pow_fmpz(decball_t res, const decball_t x, const fmpz_t y, gr_ctx_t ctx);
int _decball_pow_fmpz_binexp(decball_t res, const decball_t x, const fmpz_t y, gr_ctx_t ctx);

int decball_add_round(decball_t res, const decball_t x, const decball_t y, slong prec, gr_ctx_t ctx);
int decball_sub_round(decball_t res, const decball_t x, const decball_t y, slong prec, gr_ctx_t ctx);
int decball_mul_round(decball_t res, const decball_t x, const decball_t y, slong prec, gr_ctx_t ctx);
int decball_div_round(decball_t res, const decball_t x, const decball_t y, slong prec, gr_ctx_t ctx);
int decball_inv_round(decball_t res, const decball_t x, slong prec, gr_ctx_t ctx);
int decball_sqrt_round(decball_t res, const decball_t x, slong prec, gr_ctx_t ctx);

int decball_add_error_decmag(decball_t res, const decmag_t err, gr_ctx_t ctx);
void _decball_add_rounding_error(decball_t res, const decimal_rounding_info * info, const decmag_t precise_err, slong prec, gr_ctx_t ctx);
int decball_add_error_10exp_si(decball_t res, slong e, gr_ctx_t ctx);
int decball_add_error_decfloat(decball_t res, const decfloat_t err, gr_ctx_t ctx);
int decball_trim(decball_t res, const decball_t x, gr_ctx_t ctx);
slong decball_rel_accuracy_digits(const decball_t x, gr_ctx_t ctx);

/* Transcendental functions via arb roundtrip */
int decball_via_arb(decball_t res, const decball_t x, int (*func)(arb_t, const arb_t, slong), gr_ctx_t ctx);
int decball_pow(decball_t res, const decball_t x, const decball_t y, gr_ctx_t ctx);

int decfloat_pow(decfloat_t res, const decfloat_t x, const decfloat_t y, gr_ctx_t ctx);

/* Common elementary functions (functions.c); the other elementary and
   special functions are available through the generic-ring interface */
int decfloat_exp(decfloat_t res, const decfloat_t x, gr_ctx_t ctx);
int decball_exp(decball_t res, const decball_t x, gr_ctx_t ctx);
int decfloat_expm1(decfloat_t res, const decfloat_t x, gr_ctx_t ctx);
int decball_expm1(decball_t res, const decball_t x, gr_ctx_t ctx);
int decfloat_log1p(decfloat_t res, const decfloat_t x, gr_ctx_t ctx);
int decball_log1p(decball_t res, const decball_t x, gr_ctx_t ctx);
int decfloat_sin(decfloat_t res, const decfloat_t x, gr_ctx_t ctx);
int decball_sin(decball_t res, const decball_t x, gr_ctx_t ctx);
int decfloat_cos(decfloat_t res, const decfloat_t x, gr_ctx_t ctx);
int decball_cos(decball_t res, const decball_t x, gr_ctx_t ctx);
int decfloat_tan(decfloat_t res, const decfloat_t x, gr_ctx_t ctx);
int decball_tan(decball_t res, const decball_t x, gr_ctx_t ctx);
int decfloat_asin(decfloat_t res, const decfloat_t x, gr_ctx_t ctx);
int decball_asin(decball_t res, const decball_t x, gr_ctx_t ctx);
int decfloat_atan(decfloat_t res, const decfloat_t x, gr_ctx_t ctx);
int decball_atan(decball_t res, const decball_t x, gr_ctx_t ctx);
int decfloat_sinh(decfloat_t res, const decfloat_t x, gr_ctx_t ctx);
int decball_sinh(decball_t res, const decball_t x, gr_ctx_t ctx);
int decfloat_cosh(decfloat_t res, const decfloat_t x, gr_ctx_t ctx);
int decball_cosh(decball_t res, const decball_t x, gr_ctx_t ctx);
int decfloat_tanh(decfloat_t res, const decfloat_t x, gr_ctx_t ctx);
int decball_tanh(decball_t res, const decball_t x, gr_ctx_t ctx);
int decfloat_asinh(decfloat_t res, const decfloat_t x, gr_ctx_t ctx);
int decball_asinh(decball_t res, const decball_t x, gr_ctx_t ctx);
int decfloat_atanh(decfloat_t res, const decfloat_t x, gr_ctx_t ctx);
int decball_atanh(decball_t res, const decball_t x, gr_ctx_t ctx);
int decfloat_acos(decfloat_t res, const decfloat_t x, gr_ctx_t ctx);
int decball_acos(decball_t res, const decball_t x, gr_ctx_t ctx);
int decfloat_acosh(decfloat_t res, const decfloat_t x, gr_ctx_t ctx);
int decball_acosh(decball_t res, const decball_t x, gr_ctx_t ctx);
int decfloat_log(decfloat_t res, const decfloat_t x, gr_ctx_t ctx);
int decball_log(decball_t res, const decball_t x, gr_ctx_t ctx);
int decfloat_sin_cos(decfloat_t res1, decfloat_t res2, const decfloat_t x, gr_ctx_t ctx);
int decball_sin_cos(decball_t res1, decball_t res2, const decball_t x, gr_ctx_t ctx);
int decfloat_sinh_cosh(decfloat_t res1, decfloat_t res2, const decfloat_t x, gr_ctx_t ctx);
int decball_sinh_cosh(decball_t res1, decball_t res2, const decball_t x, gr_ctx_t ctx);
int decfloat_atan2(decfloat_t res, const decfloat_t y, const decfloat_t x, gr_ctx_t ctx);
int decball_atan2(decball_t res, const decball_t y, const decball_t x, gr_ctx_t ctx);
int decfloat_pi(decfloat_t res, gr_ctx_t ctx);
int decball_pi(decball_t res, gr_ctx_t ctx);
int _decfloat_arb_gr(decfloat_ptr * res, slong nres, decfloat_srcptr * args, slong nargs, int has_flag, int flag, int method, gr_ctx_t ctx);
int _decball_arb_gr(decball_ptr * res, slong nres, decball_srcptr * args, slong nargs, int has_flag, int flag, int method, gr_ctx_t ctx);
int _decfloat_round_with_tail(decfloat_t res, const decfloat_t S, int tail_sign, slong tail_exp, slong prec, int rnd, gr_ctx_t ctx);
slong _decfloat_sci_exp_clamped(const decfloat_t x, gr_ctx_t ctx);
int _decfloat_round_ball(decfloat_t res, const decball_t Y, slong prec, int rnd, gr_ctx_t bctx, gr_ctx_t ctx);
slong _decimal_digits_to_bits(slong prec);
int _decimal_arb_args_infeasible(int method, decfloat_srcptr * re, decfloat_srcptr * im, slong nargs, slong prec, gr_ctx_t ctx);

/* special functions by method (functions.c, cfunctions.c) */
#define DECIMAL_ARGSPEC(method, nargs, has_flag) ((method) * 16 + (nargs) * 2 + (has_flag))
int _decfloat_unary_method(decfloat_t res, const decfloat_t x, int method, gr_ctx_t ctx);
int _decfloat_unary2_method(decfloat_t res1, decfloat_t res2, const decfloat_t x, int method1, int method2, gr_ctx_t ctx);
int _decball_unary_fn(decball_t res, const decball_t x, int method, gr_ctx_t ctx);
int _decball_unary2_method(decball_t res1, decball_t res2, const decball_t x, int method, gr_ctx_t ctx);
int _decfloat_arb_args(decfloat_ptr res, decfloat_srcptr x, decfloat_srcptr y, decfloat_srcptr z, decfloat_srcptr w, int flag, int spec, gr_ctx_t ctx);
int _decball_arb_args(decball_ptr res, decball_srcptr x, decball_srcptr y, decball_srcptr z, decball_srcptr w, int flag, int spec, gr_ctx_t ctx);
int _decfloat_constant(decfloat_t res, int method, gr_ctx_t ctx);
int _decball_constant(decball_t res, int method, gr_ctx_t ctx);
int _decfloat_gamma_fmpz(decfloat_t res, const fmpz_t n, gr_ctx_t ctx);
int _decball_gamma_fmpz(decball_t res, const fmpz_t n, gr_ctx_t ctx);
int _decfloat_gamma_fmpq(decfloat_t res, const fmpq_t x, gr_ctx_t ctx);
int _decball_gamma_fmpq(decball_t res, const fmpq_t x, gr_ctx_t ctx);
int _decfloat_fac_ui(decfloat_t res, ulong n, gr_ctx_t ctx);
int _decball_fac_ui(decball_t res, ulong n, gr_ctx_t ctx);
int _decfloat_fac_fmpz(decfloat_t res, const fmpz_t n, gr_ctx_t ctx);
int _decball_fac_fmpz(decball_t res, const fmpz_t n, gr_ctx_t ctx);
int _decfloat_rising_ui(decfloat_t res, const decfloat_t x, ulong n, gr_ctx_t ctx);
int _decball_rising_ui(decball_t res, const decball_t x, ulong n, gr_ctx_t ctx);
int _decfloat_lambertw_fmpz(decfloat_t res, const decfloat_t x, const fmpz_t k, gr_ctx_t ctx);
int _decball_lambertw_fmpz(decball_t res, const decball_t x, const fmpz_t k, gr_ctx_t ctx);
#define DECFLOAT_EXP_CLAMP (WORD_MAX / 16)

/* For the internal drivers shared by many thin wrappers (special
   functions, conversions): a single out-of-line copy, neither inlined
   into the wrappers nor cloned into specialized versions. These are not
   performance-critical relative to the work they dispatch. */
#if defined(__GNUC__) && !defined(__clang__)
# define DECIMAL_DRIVER static __attribute__((noinline, noclone))
#elif defined(__GNUC__)
# define DECIMAL_DRIVER static __attribute__((noinline))
#else
# define DECIMAL_DRIVER static
#endif

/* Largest working precision (digits) in Ziv loops at precision prec with
   operands of up to operand_digits digits: arguments longer than the
   working precision are only seen exactly once the working precision
   exceeds their length. */
DECIMAL_INLINE slong
_decimal_ziv_wp_max(slong prec, slong operand_digits)
{
    return 64 * prec + 1000 + 2 * FLINT_MIN(operand_digits, DECIMAL_CONV_DIGITS_LIMIT);
}

/* Whether a gr method is an elementary function, whose evaluation at high
   precision is cheap enough (quasi-linear) to go up to the length of the
   operands in Ziv loops; other functions stop at 64 prec + 1000 digits. */
int _decimal_method_is_elementary(int method);

int decfloat_via_arb(decfloat_t res, const decfloat_t x, int (*func)(arb_t, const arb_t, slong), gr_ctx_t ctx);
void _decimal_arb_10_pow_fmpz(arb_t res, const fmpz_t t, slong wp);
int _decfloat_set_round_arb_err(decfloat_t res, decmag_ptr err, const arb_t x, slong prec, int rnd, gr_ctx_t ctx);


/* ------------------------------------------------------------------------- */
/* deccfloat (complex floating-point numbers)                                */
/* ------------------------------------------------------------------------- */

DECIMAL_INLINE void
deccfloat_init(deccfloat_t res, gr_ctx_t ctx)
{
    decfloat_init(&res->re, ctx);
    decfloat_init(&res->im, ctx);
}

/* out of line: expanded at many call sites */
void deccfloat_clear(deccfloat_t res, gr_ctx_t ctx);

DECIMAL_INLINE void
deccfloat_swap(deccfloat_t x, deccfloat_t y, gr_ctx_t ctx)
{
    FLINT_SWAP(deccfloat_struct, *x, *y);
}

DECIMAL_INLINE void
deccfloat_set_shallow(deccfloat_t res, const deccfloat_t x, gr_ctx_t ctx)
{
    *res = *x;
}

DECIMAL_INLINE int
deccfloat_zero(deccfloat_t res, gr_ctx_t ctx)
{
    decfloat_zero(&res->re, ctx);
    decfloat_zero(&res->im, ctx);
    return GR_SUCCESS;
}

DECIMAL_INLINE decfloat_ptr deccfloat_realref(deccfloat_t x) { return &x->re; }
DECIMAL_INLINE decfloat_ptr deccfloat_imagref(deccfloat_t x) { return &x->im; }

DECIMAL_INLINE int _deccfloat_is_real(const deccfloat_t x) { return DECFLOAT_IS_ZERO(&x->im); }
DECIMAL_INLINE int _deccfloat_is_imaginary(const deccfloat_t x) { return DECFLOAT_IS_ZERO(&x->re) && !DECFLOAT_IS_ZERO(&x->im); }
DECIMAL_INLINE int _deccfloat_is_finite(const deccfloat_t x) { return DECFLOAT_IS_FINITE(&x->re) && DECFLOAT_IS_FINITE(&x->im); }
DECIMAL_INLINE int _deccfloat_is_nan(const deccfloat_t x) { return DECFLOAT_IS_NAN(&x->re) || DECFLOAT_IS_NAN(&x->im); }

int deccfloat_one(deccfloat_t res, gr_ctx_t ctx);
int deccfloat_neg_one(deccfloat_t res, gr_ctx_t ctx);
int deccfloat_i(deccfloat_t res, gr_ctx_t ctx);
int deccfloat_nan(deccfloat_t res, gr_ctx_t ctx);
int deccfloat_pos_inf(deccfloat_t res, gr_ctx_t ctx);
int deccfloat_neg_inf(deccfloat_t res, gr_ctx_t ctx);

truth_t deccfloat_is_zero(const deccfloat_t x, gr_ctx_t ctx);
truth_t deccfloat_is_one(const deccfloat_t x, gr_ctx_t ctx);
truth_t deccfloat_is_neg_one(const deccfloat_t x, gr_ctx_t ctx);
truth_t deccfloat_equal(const deccfloat_t x, const deccfloat_t y, gr_ctx_t ctx);
truth_t deccfloat_is_integer(const deccfloat_t x, gr_ctx_t ctx);

/* Assignment (rounded to the context precision; the imaginary part uses
   the rounding mode for imaginary parts) */
int deccfloat_set(deccfloat_t res, const deccfloat_t x, gr_ctx_t ctx);
int deccfloat_set_round(deccfloat_t res, const deccfloat_t x, slong prec, int rnd, int rnd_im, gr_ctx_t ctx);
int deccfloat_set_decfloat(deccfloat_t res, const decfloat_t x, gr_ctx_t ctx);
int deccfloat_set_decfloat_decfloat(deccfloat_t res, const decfloat_t re, const decfloat_t im, gr_ctx_t ctx);
int deccfloat_set_si(deccfloat_t res, slong x, gr_ctx_t ctx);
int deccfloat_set_ui(deccfloat_t res, ulong x, gr_ctx_t ctx);
int deccfloat_set_fmpz(deccfloat_t res, const fmpz_t x, gr_ctx_t ctx);
int deccfloat_set_fmpq(deccfloat_t res, const fmpq_t x, gr_ctx_t ctx);
int deccfloat_set_d(deccfloat_t res, double x, gr_ctx_t ctx);
int deccfloat_set_str(deccfloat_t res, const char * s, gr_ctx_t ctx);
int deccfloat_set_other(deccfloat_t res, gr_srcptr x, gr_ctx_t x_ctx, gr_ctx_t ctx);
int deccfloat_set_fmpz_10exp_fmpz(deccfloat_t res, const fmpz_t m, const fmpz_t e, gr_ctx_t ctx);
int deccfloat_set_fmpz_2exp_fmpz(deccfloat_t res, const fmpz_t m, const fmpz_t e, gr_ctx_t ctx);
int deccfloat_set_acf(deccfloat_t res, const acf_t x, gr_ctx_t ctx);
int deccfloat_set_acb(deccfloat_t res, const acb_t x, gr_ctx_t ctx);
int deccfloat_get_acf(acf_t res, const deccfloat_t x, slong prec_bits, int rnd, gr_ctx_t ctx);
int deccfloat_get_acb(acb_t res, const deccfloat_t x, slong prec_bits, gr_ctx_t ctx);
int deccfloat_get_fmpz(fmpz_t res, const deccfloat_t x, gr_ctx_t ctx);
int deccfloat_get_fmpq(fmpq_t res, const deccfloat_t x, gr_ctx_t ctx);
int deccfloat_get_si(slong * res, const deccfloat_t x, gr_ctx_t ctx);
int deccfloat_get_ui(ulong * res, const deccfloat_t x, gr_ctx_t ctx);
int deccfloat_get_d(double * res, const deccfloat_t x, gr_ctx_t ctx);
int deccfloat_get_re(decfloat_t res, const deccfloat_t x, gr_ctx_t ctx);
int deccfloat_get_im(decfloat_t res, const deccfloat_t x, gr_ctx_t ctx);

char * deccfloat_get_str(const deccfloat_t x, gr_ctx_t ctx);
int deccfloat_write(gr_stream_t out, const deccfloat_t x, gr_ctx_t ctx);
int deccfloat_randtest(deccfloat_t res, flint_rand_t state, gr_ctx_t ctx);
int deccfloat_randtest_special(deccfloat_t res, flint_rand_t state, gr_ctx_t ctx);

/* Arithmetic: correctly rounded componentwise */
int deccfloat_neg(deccfloat_t res, const deccfloat_t x, gr_ctx_t ctx);
int deccfloat_conj(deccfloat_t res, const deccfloat_t x, gr_ctx_t ctx);
int deccfloat_re(deccfloat_t res, const deccfloat_t x, gr_ctx_t ctx);
int deccfloat_im(deccfloat_t res, const deccfloat_t x, gr_ctx_t ctx);
int deccfloat_abs(deccfloat_t res, const deccfloat_t x, gr_ctx_t ctx);
int deccfloat_arg(deccfloat_t res, const deccfloat_t x, gr_ctx_t ctx);
int deccfloat_sgn(deccfloat_t res, const deccfloat_t x, gr_ctx_t ctx);
int deccfloat_csgn(deccfloat_t res, const deccfloat_t x, gr_ctx_t ctx);
int deccfloat_add(deccfloat_t res, const deccfloat_t x, const deccfloat_t y, gr_ctx_t ctx);
int deccfloat_sub(deccfloat_t res, const deccfloat_t x, const deccfloat_t y, gr_ctx_t ctx);
int deccfloat_mul(deccfloat_t res, const deccfloat_t x, const deccfloat_t y, gr_ctx_t ctx);
int deccfloat_sqr(deccfloat_t res, const deccfloat_t x, gr_ctx_t ctx);
int deccfloat_div(deccfloat_t res, const deccfloat_t x, const deccfloat_t y, gr_ctx_t ctx);
int deccfloat_inv(deccfloat_t res, const deccfloat_t x, gr_ctx_t ctx);
int deccfloat_sqrt(deccfloat_t res, const deccfloat_t x, gr_ctx_t ctx);
int deccfloat_rsqrt(deccfloat_t res, const deccfloat_t x, gr_ctx_t ctx);
int deccfloat_add_ui(deccfloat_t res, const deccfloat_t x, ulong y, gr_ctx_t ctx);
int deccfloat_add_si(deccfloat_t res, const deccfloat_t x, slong y, gr_ctx_t ctx);
int deccfloat_add_fmpz(deccfloat_t res, const deccfloat_t x, const fmpz_t y, gr_ctx_t ctx);
int deccfloat_sub_ui(deccfloat_t res, const deccfloat_t x, ulong y, gr_ctx_t ctx);
int deccfloat_sub_si(deccfloat_t res, const deccfloat_t x, slong y, gr_ctx_t ctx);
int deccfloat_sub_fmpz(deccfloat_t res, const deccfloat_t x, const fmpz_t y, gr_ctx_t ctx);
int deccfloat_mul_ui(deccfloat_t res, const deccfloat_t x, ulong y, gr_ctx_t ctx);
int deccfloat_mul_si(deccfloat_t res, const deccfloat_t x, slong y, gr_ctx_t ctx);
int deccfloat_mul_fmpz(deccfloat_t res, const deccfloat_t x, const fmpz_t y, gr_ctx_t ctx);
int deccfloat_div_ui(deccfloat_t res, const deccfloat_t x, ulong y, gr_ctx_t ctx);
int deccfloat_div_si(deccfloat_t res, const deccfloat_t x, slong y, gr_ctx_t ctx);
int deccfloat_div_fmpz(deccfloat_t res, const deccfloat_t x, const fmpz_t y, gr_ctx_t ctx);
int deccfloat_mul_decfloat(deccfloat_t res, const deccfloat_t x, const decfloat_t y, gr_ctx_t ctx);
int deccfloat_div_decfloat(deccfloat_t res, const deccfloat_t x, const decfloat_t y, gr_ctx_t ctx);
int deccfloat_mul_i(deccfloat_t res, const deccfloat_t x, gr_ctx_t ctx);
int deccfloat_div_i(deccfloat_t res, const deccfloat_t x, gr_ctx_t ctx);
int deccfloat_mul_two(deccfloat_t res, const deccfloat_t x, gr_ctx_t ctx);
int deccfloat_mul_10exp_si(deccfloat_t res, const deccfloat_t x, slong e, gr_ctx_t ctx);
int deccfloat_mul_10exp_fmpz(deccfloat_t res, const deccfloat_t x, const fmpz_t e, gr_ctx_t ctx);
int deccfloat_mul_2exp_si(deccfloat_t res, const deccfloat_t x, slong e, gr_ctx_t ctx);
int deccfloat_mul_2exp_fmpz(deccfloat_t res, const deccfloat_t x, const fmpz_t e, gr_ctx_t ctx);
int deccfloat_floor(deccfloat_t res, const deccfloat_t x, gr_ctx_t ctx);
int deccfloat_ceil(deccfloat_t res, const deccfloat_t x, gr_ctx_t ctx);
int deccfloat_trunc(deccfloat_t res, const deccfloat_t x, gr_ctx_t ctx);
int deccfloat_nint(deccfloat_t res, const deccfloat_t x, gr_ctx_t ctx);
int deccfloat_pow_ui(deccfloat_t res, const deccfloat_t x, ulong y, gr_ctx_t ctx);
int deccfloat_pow_si(deccfloat_t res, const deccfloat_t x, slong y, gr_ctx_t ctx);
int deccfloat_pow_fmpz(deccfloat_t res, const deccfloat_t x, const fmpz_t y, gr_ctx_t ctx);
int deccfloat_pow(deccfloat_t res, const deccfloat_t x, const deccfloat_t y, gr_ctx_t ctx);
int deccfloat_cmp(int * res, const deccfloat_t x, const deccfloat_t y, gr_ctx_t ctx);
int deccfloat_cmpabs(int * res, const deccfloat_t x, const deccfloat_t y, gr_ctx_t ctx);
int deccfloat_vec_dot(deccfloat_t res, const deccfloat_t initial, int subtract, deccfloat_srcptr vec1, deccfloat_srcptr vec2, slong len, gr_ctx_t ctx);
int deccfloat_vec_dot_rev(deccfloat_t res, const deccfloat_t initial, int subtract, deccfloat_srcptr vec1, deccfloat_srcptr vec2, slong len, gr_ctx_t ctx);

/* Exact arithmetic (no rounding, no exponent limits); returns GR_UNABLE
   if a special value would be produced. */
int _deccfloat_mul_exact(deccfloat_t res, const deccfloat_t x, const deccfloat_t y, gr_ctx_t ctx);
int _deccfloat_add_exact(deccfloat_t res, const deccfloat_t x, const deccfloat_t y, int subtract, gr_ctx_t ctx);
int _deccfloat_inv_exact(deccfloat_t res, const deccfloat_t x, gr_ctx_t ctx);
int _deccfloat_pow_ui_exact(deccfloat_t res, const deccfloat_t x, ulong n, gr_ctx_t ctx);
int _deccfloat_sqrt_exact(deccfloat_t res, const deccfloat_t x, gr_ctx_t ctx);
int _deccfloat_finalize(deccfloat_t res, gr_ctx_t ctx);
/* x1 y1 +/- x2 y2 computed exactly and rounded once */
int _decfloat_dot2(decfloat_t res, const decfloat_t x1, const decfloat_t y1, const decfloat_t x2, const decfloat_t y2, int subtract, slong prec, int rnd, decimal_rounding_info * info, decmag_ptr err, gr_ctx_t ctx);

/* ------------------------------------------------------------------------- */
/* deccball (complex balls)                                                  */
/* ------------------------------------------------------------------------- */

DECIMAL_INLINE void
deccball_init(deccball_t res, gr_ctx_t ctx)
{
    decball_init(&res->re, ctx);
    decball_init(&res->im, ctx);
}

/* out of line: expanded at many call sites */
void deccball_clear(deccball_t res, gr_ctx_t ctx);

DECIMAL_INLINE void
deccball_swap(deccball_t x, deccball_t y, gr_ctx_t ctx)
{
    FLINT_SWAP(deccball_struct, *x, *y);
}

DECIMAL_INLINE void
deccball_set_shallow(deccball_t res, const deccball_t x, gr_ctx_t ctx)
{
    *res = *x;
}

DECIMAL_INLINE decball_ptr deccball_realref(deccball_t x) { return &x->re; }
DECIMAL_INLINE decball_ptr deccball_imagref(deccball_t x) { return &x->im; }

DECIMAL_INLINE int _deccball_is_exact(const deccball_t x, gr_ctx_t ctx) { return DECMAG_IS_ZERO(&x->re.rad) && DECMAG_IS_ZERO(&x->im.rad); }
DECIMAL_INLINE int _deccball_is_finite(const deccball_t x, gr_ctx_t ctx) { return _decball_is_finite(&x->re, ctx) && _decball_is_finite(&x->im, ctx); }
DECIMAL_INLINE int _deccball_is_real(const deccball_t x, gr_ctx_t ctx) { return DECFLOAT_IS_ZERO(&x->im.mid) && DECMAG_IS_ZERO(&x->im.rad); }

int deccball_zero(deccball_t res, gr_ctx_t ctx);
int deccball_one(deccball_t res, gr_ctx_t ctx);
int deccball_neg_one(deccball_t res, gr_ctx_t ctx);
int deccball_i(deccball_t res, gr_ctx_t ctx);

int deccball_set(deccball_t res, const deccball_t x, gr_ctx_t ctx);
int deccball_set_round(deccball_t res, const deccball_t x, slong prec, gr_ctx_t ctx);
int deccball_set_round2(deccball_t res, const deccball_t x, slong prec, slong rad_prec, gr_ctx_t ctx);
int deccball_mid(deccball_t res, const deccball_t x, gr_ctx_t ctx);
int deccball_shell(deccball_t res, const deccball_t x, gr_ctx_t ctx);
int deccball_add_rad(deccball_t res, const deccball_t x, const deccball_t r, gr_ctx_t ctx);
int deccball_set_interval_mid_inf(deccball_t res, const deccball_t m, gr_ctx_t ctx);
int deccball_set_decball(deccball_t res, const decball_t x, gr_ctx_t ctx);
int deccball_set_decball_decball(deccball_t res, const decball_t re, const decball_t im, gr_ctx_t ctx);
int deccball_set_decfloat(deccball_t res, const decfloat_t x, gr_ctx_t ctx);
int deccball_set_deccfloat(deccball_t res, const deccfloat_t x, gr_ctx_t ctx);
int deccball_set_si(deccball_t res, slong x, gr_ctx_t ctx);
int deccball_set_ui(deccball_t res, ulong x, gr_ctx_t ctx);
int deccball_set_fmpz(deccball_t res, const fmpz_t x, gr_ctx_t ctx);
int deccball_set_fmpq(deccball_t res, const fmpq_t x, gr_ctx_t ctx);
int deccball_set_d(deccball_t res, double x, gr_ctx_t ctx);
int deccball_set_str(deccball_t res, const char * s, gr_ctx_t ctx);
int deccball_set_other(deccball_t res, gr_srcptr x, gr_ctx_t x_ctx, gr_ctx_t ctx);
int deccball_set_fmpz_10exp_fmpz(deccball_t res, const fmpz_t m, const fmpz_t e, gr_ctx_t ctx);
int deccball_set_fmpz_2exp_fmpz(deccball_t res, const fmpz_t m, const fmpz_t e, gr_ctx_t ctx);
int deccball_set_interval_mid_rad(deccball_t res, const deccball_t m, const deccball_t r, gr_ctx_t ctx);
int deccball_set_acb(deccball_t res, const acb_t x, gr_ctx_t ctx);
int deccball_get_acb(acb_t res, const deccball_t x, slong prec_bits, gr_ctx_t ctx);
int deccball_get_fmpz(fmpz_t res, const deccball_t x, gr_ctx_t ctx);
int deccball_get_fmpq(fmpq_t res, const deccball_t x, gr_ctx_t ctx);
int deccball_get_si(slong * res, const deccball_t x, gr_ctx_t ctx);
int deccball_get_ui(ulong * res, const deccball_t x, gr_ctx_t ctx);
int deccball_get_d(double * res, const deccball_t x, gr_ctx_t ctx);
int deccball_get_mid(deccfloat_t res, const deccball_t x, gr_ctx_t ctx);

int deccball_randtest(deccball_t res, flint_rand_t state, gr_ctx_t ctx);
char * deccball_get_str(const deccball_t x, gr_ctx_t ctx);
int deccball_write(gr_stream_t out, const deccball_t x, gr_ctx_t ctx);

truth_t deccball_is_zero(const deccball_t x, gr_ctx_t ctx);
truth_t deccball_is_one(const deccball_t x, gr_ctx_t ctx);
truth_t deccball_is_neg_one(const deccball_t x, gr_ctx_t ctx);
truth_t deccball_equal(const deccball_t x, const deccball_t y, gr_ctx_t ctx);
truth_t deccball_is_integer(const deccball_t x, gr_ctx_t ctx);
int _deccball_contains_zero(const deccball_t x, gr_ctx_t ctx);
int _deccball_contains(const deccball_t x, const deccball_t y, gr_ctx_t ctx);
int _deccball_overlaps(const deccball_t x, const deccball_t y, gr_ctx_t ctx);
int _deccball_contains_deccfloat(const deccball_t x, const deccfloat_t y, gr_ctx_t ctx);
int deccball_cmp(int * res, const deccball_t x, const deccball_t y, gr_ctx_t ctx);
int deccball_cmpabs(int * res, const deccball_t x, const deccball_t y, gr_ctx_t ctx);

int deccball_neg(deccball_t res, const deccball_t x, gr_ctx_t ctx);
int deccball_conj(deccball_t res, const deccball_t x, gr_ctx_t ctx);
int deccball_re(deccball_t res, const deccball_t x, gr_ctx_t ctx);
int deccball_im(deccball_t res, const deccball_t x, gr_ctx_t ctx);
int deccball_abs(deccball_t res, const deccball_t x, gr_ctx_t ctx);
int deccball_arg(deccball_t res, const deccball_t x, gr_ctx_t ctx);
int deccball_sgn(deccball_t res, const deccball_t x, gr_ctx_t ctx);
int deccball_csgn(deccball_t res, const deccball_t x, gr_ctx_t ctx);
int deccball_add(deccball_t res, const deccball_t x, const deccball_t y, gr_ctx_t ctx);
int deccball_sub(deccball_t res, const deccball_t x, const deccball_t y, gr_ctx_t ctx);
int deccball_mul(deccball_t res, const deccball_t x, const deccball_t y, gr_ctx_t ctx);
int deccball_sqr(deccball_t res, const deccball_t x, gr_ctx_t ctx);
int deccball_div(deccball_t res, const deccball_t x, const deccball_t y, gr_ctx_t ctx);
int deccball_inv(deccball_t res, const deccball_t x, gr_ctx_t ctx);
int deccball_sqrt(deccball_t res, const deccball_t x, gr_ctx_t ctx);
int deccball_rsqrt(deccball_t res, const deccball_t x, gr_ctx_t ctx);
int deccball_add_ui(deccball_t res, const deccball_t x, ulong y, gr_ctx_t ctx);
int deccball_add_si(deccball_t res, const deccball_t x, slong y, gr_ctx_t ctx);
int deccball_add_fmpz(deccball_t res, const deccball_t x, const fmpz_t y, gr_ctx_t ctx);
int deccball_sub_ui(deccball_t res, const deccball_t x, ulong y, gr_ctx_t ctx);
int deccball_sub_si(deccball_t res, const deccball_t x, slong y, gr_ctx_t ctx);
int deccball_sub_fmpz(deccball_t res, const deccball_t x, const fmpz_t y, gr_ctx_t ctx);
int deccball_mul_ui(deccball_t res, const deccball_t x, ulong y, gr_ctx_t ctx);
int deccball_mul_si(deccball_t res, const deccball_t x, slong y, gr_ctx_t ctx);
int deccball_mul_fmpz(deccball_t res, const deccball_t x, const fmpz_t y, gr_ctx_t ctx);
int deccball_div_ui(deccball_t res, const deccball_t x, ulong y, gr_ctx_t ctx);
int deccball_div_si(deccball_t res, const deccball_t x, slong y, gr_ctx_t ctx);
int deccball_div_fmpz(deccball_t res, const deccball_t x, const fmpz_t y, gr_ctx_t ctx);
int deccball_mul_decball(deccball_t res, const deccball_t x, const decball_t y, gr_ctx_t ctx);
int deccball_div_decball(deccball_t res, const deccball_t x, const decball_t y, gr_ctx_t ctx);
int deccball_mul_i(deccball_t res, const deccball_t x, gr_ctx_t ctx);
int deccball_div_i(deccball_t res, const deccball_t x, gr_ctx_t ctx);
int deccball_mul_two(deccball_t res, const deccball_t x, gr_ctx_t ctx);
int deccball_mul_10exp_si(deccball_t res, const deccball_t x, slong e, gr_ctx_t ctx);
int deccball_mul_10exp_fmpz(deccball_t res, const deccball_t x, const fmpz_t e, gr_ctx_t ctx);
int deccball_mul_2exp_si(deccball_t res, const deccball_t x, slong e, gr_ctx_t ctx);
int deccball_mul_2exp_fmpz(deccball_t res, const deccball_t x, const fmpz_t e, gr_ctx_t ctx);
int deccball_floor(deccball_t res, const deccball_t x, gr_ctx_t ctx);
int deccball_ceil(deccball_t res, const deccball_t x, gr_ctx_t ctx);
int deccball_trunc(deccball_t res, const deccball_t x, gr_ctx_t ctx);
int deccball_nint(deccball_t res, const deccball_t x, gr_ctx_t ctx);
int deccball_pow_ui(deccball_t res, const deccball_t x, ulong y, gr_ctx_t ctx);
int deccball_pow_si(deccball_t res, const deccball_t x, slong y, gr_ctx_t ctx);
int deccball_pow_fmpz(deccball_t res, const deccball_t x, const fmpz_t y, gr_ctx_t ctx);
int deccball_pow(deccball_t res, const deccball_t x, const deccball_t y, gr_ctx_t ctx);
int deccball_add_error_decmag(deccball_t res, const decmag_t err, gr_ctx_t ctx);
int deccball_trim(deccball_t res, const deccball_t x, gr_ctx_t ctx);
slong deccball_rel_accuracy_digits(const deccball_t x, gr_ctx_t ctx);
int _decball_dot2(decball_t res, const decball_t x1, const decball_t y1, const decball_t x2, const decball_t y2, int subtract, slong prec, gr_ctx_t ctx);

/* Common elementary functions of complex numbers (cfunctions.c) */
int deccfloat_exp(deccfloat_t res, const deccfloat_t x, gr_ctx_t ctx);
int deccball_exp(deccball_t res, const deccball_t x, gr_ctx_t ctx);
int deccfloat_expm1(deccfloat_t res, const deccfloat_t x, gr_ctx_t ctx);
int deccball_expm1(deccball_t res, const deccball_t x, gr_ctx_t ctx);
int deccfloat_log(deccfloat_t res, const deccfloat_t x, gr_ctx_t ctx);
int deccball_log(deccball_t res, const deccball_t x, gr_ctx_t ctx);
int deccfloat_log1p(deccfloat_t res, const deccfloat_t x, gr_ctx_t ctx);
int deccball_log1p(deccball_t res, const deccball_t x, gr_ctx_t ctx);
int deccfloat_sin(deccfloat_t res, const deccfloat_t x, gr_ctx_t ctx);
int deccball_sin(deccball_t res, const deccball_t x, gr_ctx_t ctx);
int deccfloat_cos(deccfloat_t res, const deccfloat_t x, gr_ctx_t ctx);
int deccball_cos(deccball_t res, const deccball_t x, gr_ctx_t ctx);
int deccfloat_tan(deccfloat_t res, const deccfloat_t x, gr_ctx_t ctx);
int deccball_tan(deccball_t res, const deccball_t x, gr_ctx_t ctx);
int deccfloat_sinh(deccfloat_t res, const deccfloat_t x, gr_ctx_t ctx);
int deccball_sinh(deccball_t res, const deccball_t x, gr_ctx_t ctx);
int deccfloat_cosh(deccfloat_t res, const deccfloat_t x, gr_ctx_t ctx);
int deccball_cosh(deccball_t res, const deccball_t x, gr_ctx_t ctx);
int deccfloat_tanh(deccfloat_t res, const deccfloat_t x, gr_ctx_t ctx);
int deccball_tanh(deccball_t res, const deccball_t x, gr_ctx_t ctx);
int deccfloat_asin(deccfloat_t res, const deccfloat_t x, gr_ctx_t ctx);
int deccball_asin(deccball_t res, const deccball_t x, gr_ctx_t ctx);
int deccfloat_acos(deccfloat_t res, const deccfloat_t x, gr_ctx_t ctx);
int deccball_acos(deccball_t res, const deccball_t x, gr_ctx_t ctx);
int deccfloat_atan(deccfloat_t res, const deccfloat_t x, gr_ctx_t ctx);
int deccball_atan(deccball_t res, const deccball_t x, gr_ctx_t ctx);
int deccfloat_asinh(deccfloat_t res, const deccfloat_t x, gr_ctx_t ctx);
int deccball_asinh(deccball_t res, const deccball_t x, gr_ctx_t ctx);
int deccfloat_acosh(deccfloat_t res, const deccfloat_t x, gr_ctx_t ctx);
int deccball_acosh(deccball_t res, const deccball_t x, gr_ctx_t ctx);
int deccfloat_atanh(deccfloat_t res, const deccfloat_t x, gr_ctx_t ctx);
int deccball_atanh(deccball_t res, const deccball_t x, gr_ctx_t ctx);
int deccfloat_sin_cos(deccfloat_t res1, deccfloat_t res2, const deccfloat_t x, gr_ctx_t ctx);
int deccball_sin_cos(deccball_t res1, deccball_t res2, const deccball_t x, gr_ctx_t ctx);
int deccfloat_sinh_cosh(deccfloat_t res1, deccfloat_t res2, const deccfloat_t x, gr_ctx_t ctx);
int deccball_sinh_cosh(deccball_t res1, deccball_t res2, const deccball_t x, gr_ctx_t ctx);
int deccfloat_pi(deccfloat_t res, gr_ctx_t ctx);
int deccball_pi(deccball_t res, gr_ctx_t ctx);

/* acb bridge (cfunctions.c) */
int _deccfloat_acb_gr(deccfloat_ptr * res, slong nres, deccfloat_srcptr * args, slong nargs, int has_flag, int flag, int method, gr_ctx_t ctx);
int _deccball_acb_gr(deccball_ptr * res, slong nres, deccball_srcptr * args, slong nargs, int has_flag, int flag, int method, gr_ctx_t ctx);
int _deccfloat_unary_method(deccfloat_t res, const deccfloat_t x, int method, gr_ctx_t ctx);
int _deccball_unary_method(deccball_t res, const deccball_t x, int method, gr_ctx_t ctx);

#ifdef __cplusplus
}
#endif

#endif

/* Conversion from algebraic numbers, declared when qqbar.h has been
   included (before or after this header) */
#if defined(QQBAR_H) && !defined(DECIMAL_QQBAR_H)
#define DECIMAL_QQBAR_H

#ifdef __cplusplus
extern "C" {
#endif

int decfloat_set_round_qqbar(decfloat_t res, const qqbar_t x, slong prec, int rnd, gr_ctx_t ctx);
int decfloat_set_qqbar(decfloat_t res, const qqbar_t x, gr_ctx_t ctx);
int decball_set_qqbar(decball_t res, const qqbar_t x, gr_ctx_t ctx);
int deccfloat_set_qqbar(deccfloat_t res, const qqbar_t x, gr_ctx_t ctx);
int deccball_set_qqbar(deccball_t res, const qqbar_t x, gr_ctx_t ctx);

#ifdef __cplusplus
}
#endif

#endif
