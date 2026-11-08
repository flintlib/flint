/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef ARF_H
#define ARF_H

#ifdef ARF_INLINES_C
#define ARF_INLINE
#else
#define ARF_INLINE static inline
#endif

#include <gmp.h>
#include "longlong.h"
#include "mpn_extras.h"
#include "mag.h"
#include "arf_types.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef enum
{
    ARF_RND_DOWN  = 0,
    ARF_RND_UP    = 1,
    ARF_RND_FLOOR = 2,
    ARF_RND_CEIL  = 3,
    ARF_RND_NEAR  = 4,
    /* relaxed modes, supported by arf_div, arf_sqrt, arf_rsqrt and their
       variants: the result may be off by at most 1 ulp in any direction
       (FAST) or at most 0.51 ulp (ACCURATE); the return value is not a
       reliable inexact flag with these modes */
    ARF_RND_FAST = 5,
    ARF_RND_ACCURATE = 6
}
arf_rnd_t;

ARF_INLINE int
arf_rounds_down(arf_rnd_t rnd, int sgnbit)
{
    if (rnd == ARF_RND_DOWN) return 1;
    if (rnd == ARF_RND_UP) return 0;
    if (rnd == ARF_RND_FLOOR) return !sgnbit;
    return sgnbit;
}

ARF_INLINE int
arf_rounds_up(arf_rnd_t rnd, int sgnbit)
{
    if (rnd == ARF_RND_DOWN) return 0;
    if (rnd == ARF_RND_UP) return 1;
    if (rnd == ARF_RND_FLOOR) return sgnbit;
    return !sgnbit;
}

#ifdef __MPFR_H
ARF_INLINE mpfr_rnd_t
arf_rnd_to_mpfr(arf_rnd_t rnd)
{
    if (rnd == ARF_RND_DOWN) return MPFR_RNDZ;
    else if (rnd == ARF_RND_UP) return MPFR_RNDA;
    else if (rnd == ARF_RND_FLOOR) return MPFR_RNDD;
    else if (rnd == ARF_RND_CEIL) return MPFR_RNDU;
    else return MPFR_RNDN;
}
#endif

/* Allow 'infinite' precision for operations where a result can be computed exactly. */
#define ARF_PREC_EXACT WORD_MAX

#define ARF_PREC_ADD(prec,extra) ((prec) == ARF_PREC_EXACT ? ARF_PREC_EXACT : (prec) + (extra))

#define ARF_RESULT_EXACT 0
#define ARF_RESULT_INEXACT 1

/* Range where we can skip fmpz overflow checks for exponent manipulation. */
#define ARF_MAX_LAGOM_EXP MAG_MAX_LAGOM_EXP
#define ARF_MIN_LAGOM_EXP MAG_MIN_LAGOM_EXP

/* Exponent values used to encode special values. */
#define ARF_EXP_ZERO 0
#define ARF_EXP_NAN COEFF_MIN
#define ARF_EXP_POS_INF (COEFF_MIN+1)
#define ARF_EXP_NEG_INF (COEFF_MIN+2)

/* Direct access to the exponent. */
#define ARF_EXP(x) ((x)->exp)
#define ARF_EXPREF(x) (&(x)->exp)

/* Finite and with lagom big exponents. */
#define ARF_IS_LAGOM(x) (ARF_EXP(x) >= ARF_MIN_LAGOM_EXP && \
                         ARF_EXP(x) <= ARF_MAX_LAGOM_EXP)

/* More than two limbs (needs pointer). */
#define ARF_HAS_PTR(x) ((x)->size > ((2 << 1) + 1))

/* Raw size field (encodes both limb size and sign). */
#define ARF_XSIZE(x) ((x)->size)

/* Construct size field value from size in limbs and sign bit. */
#define ARF_MAKE_XSIZE(size, sgnbit) ((((slong) size) << 1) | (sgnbit))

/* The limb size, and the sign bit. */
#define ARF_SIZE(x) (ARF_XSIZE(x) >> 1)
#define ARF_SGNBIT(x) (ARF_XSIZE(x) & 1)

/* Assumes non-special value */
#define ARF_NEG(x) (ARF_XSIZE(x) ^= 1)

/* Direct access to the limb data. */
#define ARF_NOPTR_D(x)   ((x)->d.noptr.d)
#define ARF_PTR_D(x)     ((x)->d.ptr.d)
#define ARF_PTR_ALLOC(x) ((x)->d.ptr.alloc)

/* Encoding for special values. */
#define ARF_IS_SPECIAL(x) (ARF_XSIZE(x) == 0)

/* Value is +/- a power of two */
#define ARF_IS_POW2(x) (ARF_SIZE(x) == 1) && (ARF_NOPTR_D(x)[0] == LIMB_TOP)

/* To set a special value, first call this and then set the exponent. */
#define ARF_MAKE_SPECIAL(x)         \
    do {                            \
        fmpz_clear(ARF_EXPREF(x));  \
        ARF_DEMOTE(x);              \
        ARF_XSIZE(x) = 0;           \
    } while (0)


void _arf_promote(arf_t x, slong n);

void _arf_demote(arf_t x);


/* Warning: does not set size! -- also doesn't demote exponent. */
#define ARF_DEMOTE(x)                 \
    do {                              \
        if (ARF_HAS_PTR(x))           \
            _arf_demote(x);           \
    } while (0)

/* Get mpn pointer and size (xptr, xn) for read-only use. */
#define ARF_GET_MPN_READONLY(xptr, xn, x)   \
    do {                                    \
        xn = ARF_SIZE(x);                   \
        if (xn <= ARF_NOPTR_LIMBS)          \
            xptr = ARF_NOPTR_D(x);          \
        else                                \
            xptr = ARF_PTR_D(x);            \
    } while (0)

/* Assumes non-special! */
#define ARF_GET_TOP_LIMB(lmb, x)                     \
    do {                                             \
        nn_srcptr __xptr;                            \
        slong __xn;                              \
        ARF_GET_MPN_READONLY(__xptr, __xn, (x));     \
        (lmb) = __xptr[__xn - 1];                    \
    } while (0)

/* Get mpn pointer xptr for writing *exactly* xn limbs to x. */
#define ARF_GET_MPN_WRITE(xptr, xn, x)                      \
    do {                                                    \
        slong __xn = (xn);                              \
        if ((__xn) <= ARF_NOPTR_LIMBS)                      \
        {                                                   \
            ARF_DEMOTE(x);                                  \
            xptr = ARF_NOPTR_D(x);                          \
        }                                                   \
        else                                                \
        {                                                   \
            if (!ARF_HAS_PTR(x))                            \
            {                                               \
                _arf_promote(x, __xn);                      \
            }                                               \
            else if (ARF_PTR_ALLOC(x) < (__xn))             \
            {                                               \
                ARF_PTR_D(x) = (nn_ptr)                     \
                        flint_realloc(ARF_PTR_D(x),         \
                        (xn) * sizeof(ulong));          \
                ARF_PTR_ALLOC(x) = (__xn);                  \
            }                                               \
            xptr = ARF_PTR_D(x);                            \
        }                                                   \
        ARF_XSIZE(x) = ARF_MAKE_XSIZE(__xn, 0);             \
    } while (0)

ARF_INLINE void
arf_init(arf_t x)
{
    fmpz_init(ARF_EXPREF(x));
    ARF_XSIZE(x) = 0;
}

void arf_clear(arf_t x);

arf_ptr _arf_vec_init(slong n);
void _arf_vec_clear(arf_ptr v, slong n);

ARF_INLINE void
arf_zero(arf_t x)
{
    ARF_MAKE_SPECIAL(x);
    ARF_EXP(x) = ARF_EXP_ZERO;
}

ARF_INLINE void
arf_pos_inf(arf_t x)
{
    ARF_MAKE_SPECIAL(x);
    ARF_EXP(x) = ARF_EXP_POS_INF;
}

ARF_INLINE void
arf_neg_inf(arf_t x)
{
    ARF_MAKE_SPECIAL(x);
    ARF_EXP(x) = ARF_EXP_NEG_INF;
}

ARF_INLINE void
arf_nan(arf_t x)
{
    ARF_MAKE_SPECIAL(x);
    ARF_EXP(x) = ARF_EXP_NAN;
}

ARF_INLINE int
arf_is_special(const arf_t x)
{
    return ARF_IS_SPECIAL(x);
}

ARF_INLINE int
arf_is_zero(const arf_t x)
{
    return ARF_IS_SPECIAL(x) && (ARF_EXP(x) == ARF_EXP_ZERO);
}

ARF_INLINE int
arf_is_pos_inf(const arf_t x)
{
    return ARF_IS_SPECIAL(x) && (ARF_EXP(x) == ARF_EXP_POS_INF);
}

ARF_INLINE int
arf_is_neg_inf(const arf_t x)
{
    return ARF_IS_SPECIAL(x) && (ARF_EXP(x) == ARF_EXP_NEG_INF);
}

ARF_INLINE int
arf_is_nan(const arf_t x)
{
    return ARF_IS_SPECIAL(x) && (ARF_EXP(x) == ARF_EXP_NAN);
}

ARF_INLINE int
arf_is_normal(const arf_t x)
{
    return !ARF_IS_SPECIAL(x);
}

ARF_INLINE int
arf_is_finite(const arf_t x)
{
    return !ARF_IS_SPECIAL(x) || (ARF_EXP(x) == ARF_EXP_ZERO);
}

ARF_INLINE int
arf_is_inf(const arf_t x)
{
    return ARF_IS_SPECIAL(x) && (ARF_EXP(x) == ARF_EXP_POS_INF ||
                                 ARF_EXP(x) == ARF_EXP_NEG_INF);
}

ARF_INLINE void
arf_one(arf_t x)
{
    fmpz_clear(ARF_EXPREF(x));
    ARF_DEMOTE(x);
    ARF_EXP(x) = 1;
    ARF_XSIZE(x) = ARF_MAKE_XSIZE(1, 0);
    ARF_NOPTR_D(x)[0] = LIMB_TOP;
}

ARF_INLINE int
arf_is_one(const arf_t x)
{
    return (ARF_EXP(x) == 1) && (ARF_XSIZE(x) == ARF_MAKE_XSIZE(1, 0))
                             && ARF_NOPTR_D(x)[0] == LIMB_TOP;
}

ARF_INLINE int
arf_sgn(const arf_t x)
{
    if (arf_is_special(x))
    {
        if (arf_is_zero(x) || arf_is_nan(x))
            return 0;
        return arf_is_pos_inf(x) ? 1 : -1;
    }
    else
    {
        return ARF_SGNBIT(x) ? -1 : 1;
    }
}

int arf_cmp(const arf_t x, const arf_t y);

int arf_cmpabs(const arf_t x, const arf_t y);

int arf_cmpabs_ui(const arf_t x, ulong y);

int arf_cmpabs_d(const arf_t x, double y);

int arf_cmp_si(const arf_t x, slong y);

int arf_cmp_ui(const arf_t x, ulong y);

int arf_cmp_d(const arf_t x, double y);

ARF_INLINE void
arf_swap(arf_t y, arf_t x)
{
    FLINT_SWAP(arf_struct, *x, *y);
}

void arf_set(arf_t y, const arf_t x);

ARF_INLINE void
arf_neg(arf_t y, const arf_t x)
{
    arf_set(y, x);

    if (arf_is_special(y))
    {
        if (arf_is_pos_inf(y))
            arf_neg_inf(y);
        else if (arf_is_neg_inf(y))
            arf_pos_inf(y);
    }
    else
    {
        ARF_NEG(y);
    }
}

ARF_INLINE void
arf_init_set_ui(arf_t x, ulong v)
{
    if (v == 0)
    {
        ARF_EXP(x) = ARF_EXP_ZERO;
        ARF_XSIZE(x) = 0;
    }
    else
    {
        unsigned int c;
        c = flint_clz(v);
        ARF_EXP(x) = FLINT_BITS - c;
        ARF_NOPTR_D(x)[0] = v << c;
        ARF_XSIZE(x) = ARF_MAKE_XSIZE(1, 0);
    }
}

ARF_INLINE void
arf_init_set_si(arf_t x, slong v)
{
    arf_init_set_ui(x, FLINT_UABS(v));
    if (v < 0)
        ARF_NEG(x);
}

ARF_INLINE void
arf_set_ui(arf_t x, ulong v)
{
    ARF_DEMOTE(x);
    _fmpz_demote(ARF_EXPREF(x));

    if (v == 0)
    {
        ARF_EXP(x) = ARF_EXP_ZERO;
        ARF_XSIZE(x) = 0;
    }
    else
    {
        unsigned int c;
        c = flint_clz(v);
        ARF_EXP(x) = FLINT_BITS - c;
        ARF_NOPTR_D(x)[0] = v << c;
        ARF_XSIZE(x) = ARF_MAKE_XSIZE(1, 0);
    }
}

ARF_INLINE void
arf_set_si(arf_t x, slong v)
{
    arf_set_ui(x, FLINT_UABS(v));
    if (v < 0)
        ARF_NEG(x);
}

ARF_INLINE void
arf_init_set_shallow(arf_t z, const arf_t x)
{
    *z = *x;
}

ARF_INLINE void
arf_init_neg_shallow(arf_t z, const arf_t x)
{
    *z = *x;
    arf_neg(z, z);
}

ARF_INLINE void
arf_init_set_mag_shallow(arf_t y, const mag_t x)
{
    ulong t = MAG_MAN(x);
    ARF_XSIZE(y) = ARF_MAKE_XSIZE(t != 0, 0);
    ARF_EXP(y) = MAG_EXP(x);
    ARF_NOPTR_D(y)[0] = t << (FLINT_BITS - MAG_BITS);
}

ARF_INLINE void
arf_init_neg_mag_shallow(arf_t z, const mag_t x)
{
    arf_init_set_mag_shallow(z, x);
    arf_neg(z, z);
}

ARF_INLINE int
arf_cmpabs_mag(const arf_t x, const mag_t y)
{
    arf_t t;
    arf_init_set_mag_shallow(t, y);  /* no need to free */
    return arf_cmpabs(x, t);
}

ARF_INLINE int
arf_mag_cmpabs(const mag_t x, const arf_t y)
{
    arf_t t;
    arf_init_set_mag_shallow(t, x);  /* no need to free */
    return arf_cmpabs(t, y);
}

/* Assumes xn > 0, x[xn-1] != 0. */
/* TBD: 1, 2 limb versions */
void arf_set_mpn(arf_t y, nn_srcptr x, slong xn, int sgnbit);

ARF_INLINE void
arf_set_mpz(arf_t y, const mpz_t x)
{
    slong size = x->_mp_size;

    if (size == 0)
        arf_zero(y);
    else
        arf_set_mpn(y, x->_mp_d, FLINT_ABS(size), size < 0);
}

ARF_INLINE void
arf_set_fmpz(arf_t y, const fmpz_t x)
{
    if (!COEFF_IS_MPZ(*x))
        arf_set_si(y, *x);
    else
        arf_set_mpz(y, COEFF_TO_PTR(*x));
}

int _arf_set_round_ui(arf_t x, ulong v, int sgnbit, slong prec, arf_rnd_t rnd);

int _arf_set_round_uiui(arf_t z, slong * fix, ulong hi, ulong lo, int sgnbit, slong prec, arf_rnd_t rnd);

int
_arf_set_round_mpn(arf_t y, slong * exp_shift, nn_srcptr x, slong xn,
    int sgnbit, slong prec, arf_rnd_t rnd);

ARF_INLINE int
arf_set_round_ui(arf_t x, ulong v, slong prec, arf_rnd_t rnd)
{
    return _arf_set_round_ui(x, v, 0, prec, rnd);
}

ARF_INLINE int
arf_set_round_si(arf_t x, slong v, slong prec, arf_rnd_t rnd)
{
    return _arf_set_round_ui(x, FLINT_UABS(v), v < 0, prec, rnd);
}

ARF_INLINE int
arf_set_round_mpz(arf_t y, const mpz_t x, slong prec, arf_rnd_t rnd)
{
    int inexact;
    slong size = x->_mp_size;
    slong fix;

    if (size == 0)
    {
        arf_zero(y);
        return 0;
    }

    inexact = _arf_set_round_mpn(y, &fix, x->_mp_d, FLINT_ABS(size),
        (size < 0), prec, rnd);
    _fmpz_demote(ARF_EXPREF(y));
    ARF_EXP(y) = FLINT_ABS(size) * FLINT_BITS + fix;
    return inexact;
}

ARF_INLINE int
arf_set_round_fmpz(arf_t y, const fmpz_t x, slong prec, arf_rnd_t rnd)
{
    if (!COEFF_IS_MPZ(*x))
        return arf_set_round_si(y, *x, prec, rnd);
    else
        return arf_set_round_mpz(y, COEFF_TO_PTR(*x), prec, rnd);
}

int arf_set_round(arf_t y, const arf_t x, slong prec, arf_rnd_t rnd);

int arf_neg_round(arf_t y, const arf_t x, slong prec, arf_rnd_t rnd);

#ifdef __MPFR_H
int arf_get_mpfr(mpfr_t x, const arf_t y, mpfr_rnd_t rnd);
void arf_set_mpfr(arf_t x, const mpfr_t y);
int _arf_call_mpfr_func(arf_ptr r1, arf_ptr r2, int (*func)(void), arf_srcptr x, arf_srcptr y, slong prec, arf_rnd_t rnd);
#endif

int arf_equal(const arf_t x, const arf_t y);
int arf_equal_si(const arf_t x, slong y);
int arf_equal_ui(const arf_t x, ulong y);
int arf_equal_d(const arf_t x, double y);

ARF_INLINE void
arf_min(arf_t z, const arf_t a, const arf_t b)
{
    if (arf_cmp(a, b) <= 0)
        arf_set(z, a);
    else
        arf_set(z, b);
}

ARF_INLINE void
arf_max(arf_t z, const arf_t a, const arf_t b)
{
    if (arf_cmp(a, b) > 0)
        arf_set(z, a);
    else
        arf_set(z, b);
}

ARF_INLINE void
arf_abs(arf_t y, const arf_t x)
{
    if (arf_sgn(x) < 0)
        arf_neg(y, x);
    else
        arf_set(y, x);
}

ARF_INLINE slong
arf_bits(const arf_t x)
{
    if (arf_is_special(x))
        return 0;
    else
    {
        nn_srcptr xp;
        slong xn;
        slong c;

        ARF_GET_MPN_READONLY(xp, xn, x);
        c = flint_ctz(xp[0]);
        return xn * FLINT_BITS - c;
    }
}

ARF_INLINE void
arf_bot(fmpz_t e, const arf_t x)
{
    if (arf_is_special(x))
        fmpz_zero(e);
    else
        fmpz_sub_si(e, ARF_EXPREF(x), arf_bits(x));
}

int arf_is_int(const arf_t x);

int arf_is_int_2exp_si(const arf_t x, slong e);

int arf_cmp_2exp_si(const arf_t x, slong e);

int arf_cmpabs_2exp_si(const arf_t x, slong e);

ARF_INLINE void
arf_set_si_2exp_si(arf_t x, slong man, slong exp)
{
    arf_set_si(x, man);
    if (man != 0)
        fmpz_add_si_inline(ARF_EXPREF(x), ARF_EXPREF(x), exp);
}

ARF_INLINE void
arf_set_ui_2exp_si(arf_t x, ulong man, slong exp)
{
    arf_set_ui(x, man);
    if (man != 0)
        fmpz_add_si_inline(ARF_EXPREF(x), ARF_EXPREF(x), exp);
}

ARF_INLINE void
arf_mul_2exp_si(arf_t y, const arf_t x, slong e)
{
    arf_set(y, x);
    if (!arf_is_special(y))
        fmpz_add_si_inline(ARF_EXPREF(y), ARF_EXPREF(y), e);
}

ARF_INLINE void
arf_mul_2exp_fmpz(arf_t y, const arf_t x, const fmpz_t e)
{
    arf_set(y, x);
    if (!arf_is_special(y))
        fmpz_add_inline(ARF_EXPREF(y), ARF_EXPREF(y), e);
}

ARF_INLINE int
arf_set_round_fmpz_2exp(arf_t y, const fmpz_t x, const fmpz_t exp, slong prec, arf_rnd_t rnd)
{
    if (fmpz_is_zero(x))
    {
        arf_zero(y);
        return 0;
    }
    else
    {
        int r = arf_set_round_fmpz(y, x, prec, rnd);
        fmpz_add_inline(ARF_EXPREF(y), ARF_EXPREF(y), exp);
        return r;
    }
}

ARF_INLINE void
arf_abs_bound_lt_2exp_fmpz(fmpz_t b, const arf_t x)
{
    if (arf_is_special(x))
        fmpz_zero(b);
    else
        fmpz_set(b, ARF_EXPREF(x));
}

ARF_INLINE void
arf_abs_bound_le_2exp_fmpz(fmpz_t b, const arf_t x)
{
    if (arf_is_special(x))
        fmpz_zero(b);
    else if (ARF_IS_POW2(x))
        fmpz_sub_ui(b, ARF_EXPREF(x), 1);
    else
        fmpz_set(b, ARF_EXPREF(x));
}

slong arf_abs_bound_lt_2exp_si(const arf_t x);

void arf_frexp(arf_t man, fmpz_t exp, const arf_t x);

void arf_get_fmpz_2exp(fmpz_t man, fmpz_t exp, const arf_t x);

int _arf_get_integer_mpn(nn_ptr y, nn_srcptr x, slong xn, slong exp);

int _arf_set_mpn_fixed(arf_t z, nn_srcptr xp, slong xn,
        slong fixn, int negative, slong prec, arf_rnd_t rnd);

int arf_get_fmpz(fmpz_t z, const arf_t x, arf_rnd_t rnd);

slong arf_get_si(const arf_t x, arf_rnd_t rnd);

int arf_get_fmpz_fixed_fmpz(fmpz_t y, const arf_t x, const fmpz_t e);

int arf_get_fmpz_fixed_si(fmpz_t y, const arf_t x, slong e);

ARF_INLINE void
arf_set_fmpz_2exp(arf_t x, const fmpz_t man, const fmpz_t exp)
{
    arf_set_fmpz(x, man);
    if (!arf_is_zero(x))
        fmpz_add_inline(ARF_EXPREF(x), ARF_EXPREF(x), exp);
}

void arf_floor(arf_t z, const arf_t x);
void arf_ceil(arf_t z, const arf_t x);
void arf_nint(arf_t z, const arf_t x);

void arf_debug(const arf_t x);

char * arf_get_str(const arf_t x, slong d);
char * arf_get_str_rnd(const arf_t x, slong d, arf_rnd_t rnd);

#ifdef FLINT_HAVE_FILE
void arf_fprint(FILE * file, const arf_t x);
void arf_fprintd(FILE * file, const arf_t y, slong d);
#endif

void arf_print(const arf_t x);
void arf_printd(const arf_t y, slong d);

void arf_randtest(arf_t x, flint_rand_t state, slong bits, slong mag_bits);

void arf_randtest_not_zero(arf_t x, flint_rand_t state, slong bits, slong mag_bits);

void arf_randtest_special(arf_t x, flint_rand_t state, slong bits, slong mag_bits);

void arf_urandom(arf_t x, flint_rand_t state, slong bits, arf_rnd_t rnd);

#define ARF_MUL_STACK_ALLOC 40
#define ARF_MUL_TLS_ALLOC 1000

extern FLINT_TLS_PREFIX nn_ptr __arf_mul_tmp;
extern FLINT_TLS_PREFIX slong __arf_mul_alloc;

extern void _arf_mul_tmp_cleanup(void);

#define ARF_MUL_TMP_DECL \
    ulong tmp_stack[ARF_MUL_STACK_ALLOC]; \

#define ARF_MUL_TMP_ALLOC(tmp, alloc) \
    if ((alloc) <= ARF_MUL_STACK_ALLOC) \
    { \
        tmp = tmp_stack; \
    } \
    else if ((alloc) <= ARF_MUL_TLS_ALLOC) \
    { \
        if (__arf_mul_alloc < (alloc)) \
        { \
            if (__arf_mul_alloc == 0) \
            { \
                flint_register_cleanup_function(_arf_mul_tmp_cleanup); \
            } \
            __arf_mul_tmp = flint_realloc(__arf_mul_tmp, sizeof(ulong) * (alloc)); \
            __arf_mul_alloc = (alloc); \
        } \
        tmp = __arf_mul_tmp; \
    } \
    else \
    { \
        tmp = flint_malloc(sizeof(ulong) * (alloc)); \
    }

#define ARF_MUL_TMP_FREE(tmp, alloc) \
    if ((alloc) > ARF_MUL_TLS_ALLOC) \
        flint_free(tmp);

void arf_mul_special(arf_t z, const arf_t x, const arf_t y);

int arf_mul_rnd_any(arf_ptr z, arf_srcptr x, arf_srcptr y, slong prec, arf_rnd_t rnd);

int arf_mul_rnd_down(arf_ptr z, arf_srcptr x, arf_srcptr y, slong prec);

#define arf_mul(z, x, y, prec, rnd)              \
    ((rnd == ARF_RND_DOWN)                      \
        ? arf_mul_rnd_down(z, x, y, prec)        \
        : arf_mul_rnd_any(z, x, y, prec, rnd))

/* Products of at most this many limb multiplications are formed in full
   even when part of them is discarded: cheaper than the truncation and
   normalisation of a high product (measured: up to 7 x 7 limbs for an
   output of n - 1 limbs). */
#define ARF_MUL_FULL_SMALL 64

/* whether the product of xn x yn limbs of which the top n limbs are
   wanted is computed in full rather than as a high product: small
   products, and the shapes for which the general high product would
   form the full product anyway */
ARF_INLINE int
_arf_mul_prefer_full(slong xn, slong yn, slong n)
{
    if (xn + yn <= 14 || xn * yn <= ARF_MUL_FULL_SMALL)
        return 1;
    return _flint_mpn_mulhigh_prefer_full(xn, yn, n);
}

/* the normalised high product of rn limbs into z, with its guard limb
   and the shift (0 or 1) by which it was normalised: the balanced
   kernels directly (the general functions reduce to them after a dozen
   checks, which matter at small sizes) */
ARF_INLINE mp_limb_pair_t
_arf_mpn_mulhigh_normalised(nn_ptr z, nn_srcptr x, slong xn, nn_srcptr y, slong yn, slong rn)
{
    if (x == y && xn == yn)
        return (xn == rn) ? flint_mpn_sqrhigh_n_normalised(z, x, rn)
                          : flint_mpn_sqrhigh_normalised(z, x, xn, rn);
    else if (xn == rn && yn == rn)
        return flint_mpn_mulhigh_n_normalised(z, x, y, rn);
    else
        return flint_mpn_mulhigh_normalised(z, x, xn, y, yn, rn);
}

ARF_INLINE int
arf_neg_mul(arf_t z, const arf_t x, const arf_t y, slong prec, arf_rnd_t rnd)
{
    if (arf_is_special(y))
    {
        arf_mul(z, x, y, prec, rnd);
        arf_neg(z, z);
        return 0;
    }
    else
    {
        arf_t t;
        *t = *y;
        ARF_NEG(t);
        return arf_mul(z, x, t, prec, rnd);
    }
}

ARF_INLINE int
arf_mul_ui(arf_ptr z, arf_srcptr x, ulong y, slong prec, arf_rnd_t rnd)
{
    arf_t t;
    arf_init_set_ui(t, y); /* no need to free */
    return arf_mul(z, x, t, prec, rnd);
}

ARF_INLINE int
arf_mul_si(arf_ptr z, arf_srcptr x, slong y, slong prec, arf_rnd_t rnd)
{
    arf_t t;
    arf_init_set_si(t, y); /* no need to free */
    return arf_mul(z, x, t, prec, rnd);
}

int arf_mul_mpz(arf_ptr z, arf_srcptr x, const mpz_t y, slong prec, arf_rnd_t rnd);

ARF_INLINE int
arf_mul_fmpz(arf_ptr z, arf_srcptr x, const fmpz_t y, slong prec, arf_rnd_t rnd)
{
    if (!COEFF_IS_MPZ(*y))
        return arf_mul_si(z, x, *y, prec, rnd);
    else
        return arf_mul_mpz(z, x, COEFF_TO_PTR(*y), prec, rnd);
}

#define ARF_ADD_STACK_ALLOC 40
#define ARF_ADD_TLS_ALLOC 1000

extern FLINT_TLS_PREFIX nn_ptr __arf_add_tmp;
extern FLINT_TLS_PREFIX slong __arf_add_alloc;

extern void _arf_add_tmp_cleanup(void);

#define ARF_ADD_TMP_DECL \
    ulong tmp_stack[ARF_ADD_STACK_ALLOC]; \

#define ARF_ADD_TMP_ALLOC(tmp, alloc) \
    if ((alloc) <= ARF_ADD_STACK_ALLOC) \
    { \
        tmp = tmp_stack; \
    } \
    else if ((alloc) <= ARF_ADD_TLS_ALLOC) \
    { \
        if (__arf_add_alloc < (alloc)) \
        { \
            if (__arf_add_alloc == 0) \
            { \
                flint_register_cleanup_function(_arf_add_tmp_cleanup); \
            } \
            __arf_add_tmp = flint_realloc(__arf_add_tmp, sizeof(ulong) * (alloc)); \
            __arf_add_alloc = (alloc); \
        } \
        tmp = __arf_add_tmp; \
    } \
    else \
    { \
        tmp = flint_malloc(sizeof(ulong) * (alloc)); \
    }

#define ARF_ADD_TMP_FREE(tmp, alloc) \
    if ((alloc) > ARF_ADD_TLS_ALLOC) \
        flint_free(tmp);

int _arf_add_mpn(arf_t z, nn_srcptr xp, slong xn, int xsgnbit,
    const fmpz_t xexp, nn_srcptr yp, slong yn, int ysgnbit,
    flint_bitcnt_t shift, slong prec, arf_rnd_t rnd);

int arf_add(arf_ptr z, arf_srcptr x, arf_srcptr y, slong prec, arf_rnd_t rnd);
int arf_add_si(arf_ptr z, arf_srcptr x, slong y, slong prec, arf_rnd_t rnd);
int arf_add_ui(arf_ptr z, arf_srcptr x, ulong y, slong prec, arf_rnd_t rnd);
int arf_add_fmpz(arf_ptr z, arf_srcptr x, const fmpz_t y, slong prec, arf_rnd_t rnd);

int arf_add_fmpz_2exp(arf_ptr z, arf_srcptr x, const fmpz_t y, const fmpz_t exp, slong prec, arf_rnd_t rnd);

int arf_sub(arf_ptr z, arf_srcptr x, arf_srcptr y, slong prec, arf_rnd_t rnd);
int arf_sub_si(arf_ptr z, arf_srcptr x, slong y, slong prec, arf_rnd_t rnd);
int arf_sub_ui(arf_ptr z, arf_srcptr x, ulong y, slong prec, arf_rnd_t rnd);
int arf_sub_fmpz(arf_ptr z, arf_srcptr x, const fmpz_t y, slong prec, arf_rnd_t rnd);

int arf_addmul(arf_ptr z, arf_srcptr x, arf_srcptr y, slong prec, arf_rnd_t rnd);

ARF_INLINE int
arf_addmul_ui(arf_ptr z, arf_srcptr x, ulong y, slong prec, arf_rnd_t rnd)
{
    arf_t t;
    arf_init_set_ui(t, y); /* no need to free */
    return arf_addmul(z, x, t, prec, rnd);
}

ARF_INLINE int
arf_addmul_si(arf_ptr z, arf_srcptr x, slong y, slong prec, arf_rnd_t rnd)
{
    arf_t t;
    arf_init_set_si(t, y); /* no need to free */
    return arf_addmul(z, x, t, prec, rnd);
}

int arf_addmul_mpz(arf_ptr z, arf_srcptr x, const mpz_t y, slong prec, arf_rnd_t rnd);

ARF_INLINE int
arf_addmul_fmpz(arf_ptr z, arf_srcptr x, const fmpz_t y, slong prec, arf_rnd_t rnd)
{
    if (!COEFF_IS_MPZ(*y))
        return arf_addmul_si(z, x, *y, prec, rnd);
    else
        return arf_addmul_mpz(z, x, COEFF_TO_PTR(*y), prec, rnd);
}

int arf_submul(arf_ptr z, arf_srcptr x, arf_srcptr y, slong prec, arf_rnd_t rnd);

ARF_INLINE int
arf_submul_ui(arf_ptr z, arf_srcptr x, ulong y, slong prec, arf_rnd_t rnd)
{
    arf_t t;
    arf_init_set_ui(t, y); /* no need to free */
    return arf_submul(z, x, t, prec, rnd);
}

ARF_INLINE int
arf_submul_si(arf_ptr z, arf_srcptr x, slong y, slong prec, arf_rnd_t rnd)
{
    arf_t t;
    arf_init_set_si(t, y); /* no need to free */
    return arf_submul(z, x, t, prec, rnd);
}

int arf_submul_mpz(arf_ptr z, arf_srcptr x, const mpz_t y, slong prec, arf_rnd_t rnd);

ARF_INLINE int
arf_submul_fmpz(arf_ptr z, arf_srcptr x, const fmpz_t y, slong prec, arf_rnd_t rnd)
{
    if (!COEFF_IS_MPZ(*y))
        return arf_submul_si(z, x, *y, prec, rnd);
    else
        return arf_submul_mpz(z, x, COEFF_TO_PTR(*y), prec, rnd);
}

int arf_fma(arf_ptr res, arf_srcptr x, arf_srcptr y, arf_srcptr z, slong prec, arf_rnd_t rnd);

int arf_sosq(arf_t z, const arf_t x, const arf_t y, slong prec, arf_rnd_t rnd);

int arf_div(arf_ptr z, arf_srcptr x, arf_srcptr y, slong prec, arf_rnd_t rnd);

/*
    Floating-point kernels on limb arrays: msb-aligned mantissas with a
    radix-2 exponent (the representation of arf and nfloat). A mantissa
    is a normalised limb array {x, xn} (top bit set) standing for the
    value {x, xn} 2^(e - xn FLINT_BITS), B = 2^FLINT_BITS, for an
    exponent e which stays with the caller, as does the sign. A kernel
    writes its result as {z, rn} with rn = ceil(prec / FLINT_BITS) limbs
    (normalised, the bits below prec zero) and sets *zexp to the
    correction of the natural exponent of the operation (ex - ey for a
    quotient, and so on). The rounded kernels take a sign-free mode
    (ARF_RND_DOWN, ARF_RND_UP or ARF_RND_NEAR; _arf_rnd_resolve gives it)
    and return the inexact flag; the _approx kernels return rn limbs
    within an ulp of the magnitude and nothing else. Scratch beyond a
    small stack buffer lives in a thread-local buffer.
*/

/* the sign-free strict mode rounding a value of sign sgnbit in the mode
   rnd (the relaxed modes as arf_rnd_relaxed_to_strict: ARF_RND_FAST
   toward zero, ARF_RND_ACCURATE to nearest) */
ARF_INLINE arf_rnd_t
_arf_rnd_resolve(arf_rnd_t rnd, int sgnbit)
{
    if (rnd == ARF_RND_FLOOR || rnd == ARF_RND_CEIL)
        return ((rnd == ARF_RND_CEIL) ^ sgnbit) ? ARF_RND_UP : ARF_RND_DOWN;
    if (rnd >= ARF_RND_FAST)
        return (rnd == ARF_RND_ACCURATE) ? ARF_RND_NEAR : ARF_RND_DOWN;
    return rnd;
}

/* limbs of scratch held on the stack by the kernels; larger scratch comes
   from a thread-local buffer (one buffer: not for nested use) */
#define ARF_MPN_STACK_LIMBS 64

extern FLINT_TLS_PREFIX nn_ptr __arf_mpn_tmp;
extern FLINT_TLS_PREFIX slong __arf_mpn_tmp_alloc;
void _arf_mpn_tmp_cleanup(void);

ARF_INLINE nn_ptr
_arf_mpn_tmp(slong alloc)
{
    if (__arf_mpn_tmp_alloc < alloc)
    {
        if (__arf_mpn_tmp_alloc == 0)
            flint_register_cleanup_function(_arf_mpn_tmp_cleanup);
        __arf_mpn_tmp = (nn_ptr) flint_realloc(__arf_mpn_tmp, sizeof(ulong) * alloc);
        __arf_mpn_tmp_alloc = alloc;
    }
    return __arf_mpn_tmp;
}

#define ARF_MPN_TMP_DECL ulong __arf_mpn_stack[ARF_MPN_STACK_LIMBS];
#define ARF_MPN_TMP_ALLOC(tmp, alloc) \
    (tmp) = ((alloc) <= ARF_MPN_STACK_LIMBS) ? __arf_mpn_stack : _arf_mpn_tmp(alloc);

/* res (rn limbs, top bit set, the bits below sb zero) rounded up by one
   unit of 2^sb; returns the carry (the result is then B^rn / 2, a power
   of two: the exponent increases) */
ARF_INLINE int
_arf_mpn_increment(nn_ptr res, slong rn, slong sb)
{
    ulong cy = (sb == 0) ? 1 : (UWORD(1) << sb);
    slong i;

    for (i = 0; i < rn; i++)
    {
        res[i] += cy;
        cy = (res[i] < cy);
        if (!cy)
            return 0;
    }

    res[rn - 1] = UWORD(1) << (FLINT_BITS - 1);
    _flint_mpn_zero_small(res, rn - 1);
    return 1;
}

/* res decremented by one unit of 2^sb; when res was the power of two
   B^rn / 2, the result is the prec = rn FLINT_BITS - sb ones below it,
   normalised: returns -1 (the exponent decreases), else 0 */
ARF_INLINE int
_arf_mpn_decrement(nn_ptr res, slong rn, slong sb)
{
    ulong one = (sb == 0) ? 1 : (UWORD(1) << sb);
    slong i;

    if (res[rn - 1] == (UWORD(1) << (FLINT_BITS - 1)) && (rn == 1 || flint_mpn_zero_p(res, rn - 1)))
    {
        ulong mask = (sb == 0) ? ~UWORD(0) : (~UWORD(0) << sb);
        for (i = 0; i < rn; i++)
            res[i] = ~UWORD(0);
        res[0] &= mask;
        return -1;
    }

    for (i = 0; i < rn; i++)
    {
        ulong t = res[i];
        res[i] = t - one;
        if (t >= one)
            return 0;
        one = 1;
    }
    return 0;
}

/* {x, xn} (any length, x[xn - 1] != 0) rounded to prec bits into {res,
   rn}, *exp the correction of the exponent of {x, xn} B^-xn: returns
   the inexact flag; res may be x or a lower position in its buffer */
int _arf_mpn_set_round(nn_ptr res, slong * exp, nn_srcptr x, slong xn, slong prec, arf_rnd_t rnd);

/* An approximation {a, n} (top bit set) of a value T, with a guard limb g
   below it within 2^e units (two-sided, or one-sided when side = 1 (T
   above) or -1 (T below)), rounded to prec <= n FLINT_BITS bits into
   {res, rn} when the bits e .. FLINT_BITS - 2 of g certify that T rounds
   like the approximation and is inexact (neither all zero nor all
   ones): returns 1 with *exp the carry (0 or 1), or -1 when the guard
   limb cannot certify the rounding (the exact value must decide). */
int _arf_mpn_round_guarded(nn_ptr res, slong * exp, nn_srcptr a, slong n, ulong g, int e, slong prec, arf_rnd_t rnd, int side);

/* whether rounding {a, an} (an approximation of a value P within 2^e
   units of a[0], on either side) to prec bits gives the same result as
   rounding P, in every strict rounding mode, with P inexact: the bits of
   a from bit e of its lowest limb up to the bit below the rounding point
   are neither all zero nor all ones */
int _arf_rounding_certain(nn_srcptr a, slong an, slong e, slong prec);

/* x y truncated to rn = ceil(prec / FLINT_BITS) limbs, within an ulp of
   rn limbs; z must not alias an operand */
void _arf_mpn_mul_approx(nn_ptr z, slong * zexp, nn_srcptr x, slong xn, nn_srcptr y, slong yn, slong prec);

/* quotients and roots (arf/mpn_div.c): x / y and 1 / y; for the roots
   x = A 2^e with A = {x, xn} B^-xn in [1/2, 1) and t = e mod 2, so that
   sqrt(x) = {z, rn} B^-rn 2^((e + t) / 2 + *zexp) and
   1 / sqrt(x) = {z, rn} B^-rn 2^((1 - t) - (e - t) / 2 + *zexp) */
int _arf_mpn_div(nn_ptr z, slong * zexp, nn_srcptr x, slong xn, nn_srcptr y, slong yn, slong prec, arf_rnd_t rnd);
int _arf_mpn_inv(nn_ptr z, slong * zexp, nn_srcptr y, slong yn, slong prec, arf_rnd_t rnd);
int _arf_mpn_sqrt(nn_ptr z, slong * zexp, nn_srcptr x, slong xn, int t, slong prec, arf_rnd_t rnd);
int _arf_mpn_rsqrt(nn_ptr z, slong * zexp, nn_srcptr x, slong xn, int t, slong prec, arf_rnd_t rnd);
void _arf_mpn_div_approx(nn_ptr z, slong * zexp, nn_srcptr x, slong xn, nn_srcptr y, slong yn, slong prec);
void _arf_mpn_inv_approx(nn_ptr z, slong * zexp, nn_srcptr y, slong yn, slong prec);
void _arf_mpn_sqrt_approx(nn_ptr z, slong * zexp, nn_srcptr x, slong xn, int t, slong prec);
void _arf_mpn_rsqrt_approx(nn_ptr z, slong * zexp, nn_srcptr x, slong xn, int t, slong prec);

/* complex reciprocals, quotients and roots (arf/mpn_complex.c): the
   magnitudes of the parts of 1 / (a + b i), (a + b i) / (c + d i),
   sqrt(a + b i) and 1 / sqrt(a + b i), the exponents given as slongs
   and the signs of the operands as sa, ...; see the file for the
   signs of the results */
void _arf_mpn_cinv_approx(nn_ptr re, slong * reexp, nn_ptr im, slong * imexp, nn_srcptr a, slong an, slong ea, nn_srcptr b, slong bn, slong eb, slong prec);
void _arf_mpn_cdiv_approx(nn_ptr re, slong * reexp, int * reneg, nn_ptr im, slong * imexp, int * imneg, nn_srcptr a, slong an, slong ea, int sa, nn_srcptr b, slong bn, slong eb, int sb, nn_srcptr c, slong cn, slong ec, int sc, nn_srcptr d, slong dn, slong ed, int sd, slong prec);
void _arf_mpn_csqrt_approx(nn_ptr re, slong * reexp, nn_ptr im, slong * imexp, nn_srcptr a, slong an, slong ea, int sa, nn_srcptr b, slong bn, slong eb, slong prec);
void _arf_mpn_crsqrt_approx(nn_ptr re, slong * reexp, nn_ptr im, slong * imexp, nn_srcptr a, slong an, slong ea, int sa, nn_srcptr b, slong bn, slong eb, slong prec);

/* The buffer of alloc limbs for a kernel writing its result into z: the
   limbs of z (allocated), unless z is one of the operands x, y (through
   its struct or a shallow copy of its buffer) or the result may be
   short enough for the inline limbs, in which case the kernel works in
   scratch and the result is copied (ARF_KERNEL_FINISH). */
#define ARF_KERNEL_DECL \
    ulong __kstack[ARF_MUL_STACK_ALLOC]; \
    nn_ptr __ktmp = NULL;

#define ARF_KERNEL_ALLOC(zd, alloc, z, x, y) \
    do { \
        int __kalias = ((z) == (x) || (z) == (y) || (ARF_HAS_PTR(z) && \
            (ARF_PTR_D(z) == (nn_ptr) ARF_NOPTR_D(x) || ARF_PTR_D(z) == (nn_ptr) ARF_NOPTR_D(y) || \
             (ARF_HAS_PTR(x) && ARF_PTR_D(z) == ARF_PTR_D(x)) || (ARF_HAS_PTR(y) && ARF_PTR_D(z) == ARF_PTR_D(y))))); \
        if (__kalias || (alloc) <= ARF_NOPTR_LIMBS + 1) \
        { \
            __ktmp = ((alloc) <= ARF_MUL_STACK_ALLOC) ? __kstack : flint_malloc(sizeof(ulong) * (alloc)); \
            (zd) = __ktmp; \
        } \
        else \
            ARF_GET_MPN_WRITE(zd, alloc, z); \
    } while (0)

/* the n normalised limbs of the result into z (n = 0: zero) */
#define ARF_KERNEL_FINISH(z, zd, n, sgnbit) \
    do { \
        if ((n) == 0) \
            arf_zero(z); \
        else if (__ktmp != NULL) \
            _arf_set_mpn_canonical_copy(z, zd, n, sgnbit); \
        else \
            _arf_set_mpn_canonical(z, n, sgnbit); \
        if (__ktmp != NULL && __ktmp != __kstack) \
            flint_free(__ktmp); \
    } while (0)

/* z from the n normalised limbs a kernel wrote into its buffer (of n
   limbs): zero low limbs stripped, two or fewer limbs moved into the
   inline limbs */
ARF_INLINE void
_arf_set_mpn_canonical(arf_ptr z, slong n, int sgnbit)
{
    nn_ptr d = ARF_HAS_PTR(z) ? ARF_PTR_D(z) : ARF_NOPTR_D(z);
    slong k = 0;

    while (d[k] == 0)
        k++;
    n -= k;

    if (n <= ARF_NOPTR_LIMBS && ARF_HAS_PTR(z))
    {
        /* the kept limbs move into the struct, freeing the buffer */
        ulong t0 = d[k], t1 = (n == 2) ? d[k + 1] : 0;
        _arf_demote(z);
        ARF_NOPTR_D(z)[0] = t0;
        ARF_NOPTR_D(z)[1] = t1;
    }
    else if (k != 0)
        flint_mpn_copyi(d, d + k, n);

    ARF_XSIZE(z) = ARF_MAKE_XSIZE(n, sgnbit);
}

/* z from the n normalised limbs a kernel wrote in scratch (zero low
   limbs stripped) */
ARF_INLINE void
_arf_set_mpn_canonical_copy(arf_ptr z, nn_srcptr d, slong n, int sgnbit)
{
    nn_ptr zd;

    while (d[0] == 0)
    {
        d++;
        n--;
    }

    ARF_GET_MPN_WRITE(zd, n, z);
    flint_mpn_copyi(zd, d, n);
    ARF_XSIZE(z) = ARF_MAKE_XSIZE(n, sgnbit);
}

/* strict mode with which the standard algorithms realise a relaxed mode */
ARF_INLINE arf_rnd_t arf_rnd_relaxed_to_strict(arf_rnd_t rnd)
{
    if (rnd == ARF_RND_FAST) return ARF_RND_DOWN;
    if (rnd == ARF_RND_ACCURATE) return ARF_RND_NEAR;
    return rnd;
}

ARF_INLINE int
arf_div_ui(arf_ptr z, arf_srcptr x, ulong y, slong prec, arf_rnd_t rnd)
{
    arf_t t;
    arf_init_set_ui(t, y); /* no need to free */
    return arf_div(z, x, t, prec, rnd);
}

ARF_INLINE int
arf_ui_div(arf_ptr z, ulong x, arf_srcptr y, slong prec, arf_rnd_t rnd)
{
    arf_t t;
    arf_init_set_ui(t, x); /* no need to free */
    return arf_div(z, t, y, prec, rnd);
}

ARF_INLINE int
arf_div_si(arf_ptr z, arf_srcptr x, slong y, slong prec, arf_rnd_t rnd)
{
    arf_t t;
    arf_init_set_si(t, y); /* no need to free */
    return arf_div(z, x, t, prec, rnd);
}

ARF_INLINE int
arf_si_div(arf_ptr z, slong x, arf_srcptr y, slong prec, arf_rnd_t rnd)
{
    arf_t t;
    arf_init_set_si(t, x); /* no need to free */
    return arf_div(z, t, y, prec, rnd);
}

ARF_INLINE int
arf_div_fmpz(arf_ptr z, arf_srcptr x, const fmpz_t y, slong prec, arf_rnd_t rnd)
{
    arf_t t;
    int r;
    arf_init(t);
    arf_set_fmpz(t, y);
    r = arf_div(z, x, t, prec, rnd);
    arf_clear(t);
    return r;
}

ARF_INLINE int
arf_fmpz_div(arf_ptr z, const fmpz_t x, arf_srcptr y, slong prec, arf_rnd_t rnd)
{
    arf_t t;
    int r;
    arf_init(t);
    arf_set_fmpz(t, x);
    r = arf_div(z, t, y, prec, rnd);
    arf_clear(t);
    return r;
}

ARF_INLINE int
arf_fmpz_div_fmpz(arf_ptr z, const fmpz_t x, const fmpz_t y, slong prec, arf_rnd_t rnd)
{
    arf_t t, u;
    int r;
    arf_init(t);
    arf_init(u);
    arf_set_fmpz(t, x);
    arf_set_fmpz(u, y);
    r = arf_div(z, t, u, prec, rnd);
    arf_clear(t);
    arf_clear(u);
    return r;
}

int arf_sqrt(arf_ptr z, arf_srcptr x, slong prec, arf_rnd_t rnd);

int arf_sqrt_ui(arf_t z, ulong x, slong prec, arf_rnd_t rnd);

int arf_sqrt_fmpz(arf_t z, const fmpz_t x, slong prec, arf_rnd_t rnd);

int arf_rsqrt(arf_ptr z, arf_srcptr x, slong prec, arf_rnd_t rnd);

int arf_root(arf_ptr z, arf_srcptr x, ulong k, slong prec, arf_rnd_t rnd);

/* Magnitude bounds */

void arf_get_mag(mag_t y, const arf_t x);

void arf_get_mag_lower(mag_t y, const arf_t x);

ARF_INLINE void
arf_set_mag(arf_t y, const mag_t x)
{
    if (mag_is_zero(x))
    {
        arf_zero(y);
    }
    else if (mag_is_inf(x))
    {
        arf_pos_inf(y);
    }
    else
    {
        _fmpz_set_fast(ARF_EXPREF(y), MAG_EXPREF(x));
        ARF_DEMOTE(y);
        ARF_XSIZE(y) = ARF_MAKE_XSIZE(1, 0);
        ARF_NOPTR_D(y)[0] = MAG_MAN(x) << (FLINT_BITS - MAG_BITS);
    }
}

ARF_INLINE void
mag_init_set_arf(mag_t y, const arf_t x)
{
    mag_init(y);
    arf_get_mag(y, x);
}

ARF_INLINE void
mag_fast_init_set_arf(mag_t y, const arf_t x)
{
    if (ARF_IS_SPECIAL(x))   /* x == 0 */
    {
        mag_fast_zero(y);
    }
    else
    {
        nn_srcptr xp;
        slong xn;

        ARF_GET_MPN_READONLY(xp, xn, x);

        MAG_MAN(y) = (xp[xn - 1] >> (FLINT_BITS - MAG_BITS)) + LIMB_ONE;
        MAG_EXP(y) = ARF_EXP(x);

        MAG_FAST_ADJUST_ONE_TOO_LARGE(y);
    }
}

ARF_INLINE void
arf_mag_fast_add_ulp(mag_t z, const mag_t x, const arf_t y, slong prec)
{
    mag_fast_add_2exp_si(z, x, ARF_EXP(y) - prec);
}

ARF_INLINE void
arf_mag_add_ulp(mag_t z, const mag_t x, const arf_t y, slong prec)
{
    if (ARF_IS_SPECIAL(y))
    {
        flint_throw(FLINT_ERROR, "error: ulp error not defined for special value!\n");
    }
    else if (MAG_IS_LAGOM(z) && MAG_IS_LAGOM(x) && ARF_IS_LAGOM(y))
    {
        arf_mag_fast_add_ulp(z, x, y, prec);
    }
    else
    {
        fmpz_t e;
        fmpz_init(e);
        fmpz_sub_ui(e, ARF_EXPREF(y), prec);
        mag_add_2exp_fmpz(z, x, e);
        fmpz_clear(e);
    }
}

ARF_INLINE void
arf_mag_set_ulp(mag_t z, const arf_t y, slong prec)
{
    if (ARF_IS_SPECIAL(y))
    {
        flint_throw(FLINT_ERROR, "error: ulp error not defined for special value!\n");
    }
    else
    {
        _fmpz_add_fast(MAG_EXPREF(z), ARF_EXPREF(y), 1 - prec);
        MAG_MAN(z) = MAG_ONE_HALF;
    }
}

void arf_get_fmpq(fmpq_t y, const arf_t x);

ARF_INLINE int
arf_set_fmpq(arf_t y, const fmpq_t x, slong prec, arf_rnd_t rnd)
{
    return arf_fmpz_div_fmpz(y, fmpq_numref(x), fmpq_denref(x), prec, rnd);
}

int arf_complex_mul(arf_t e, arf_t f, const arf_t a, const arf_t b,
                                      const arf_t c, const arf_t d,
                                      slong prec, arf_rnd_t rnd);

int arf_complex_mul_fallback(arf_t e, arf_t f,
        const arf_t a, const arf_t b,
        const arf_t c, const arf_t d,
        slong prec, arf_rnd_t rnd);

/* the exact complex product (flint_mpn_mul_complex) of finite nonzero
   LAGOM operands, and the approximate product of long balanced operands
   for the dot products (only the bits above 2^bottom matter: the
   complex high product kernel at the length they need, the error below
   2^bottom; returns 0 when the shape is unsuitable) */
void _arf_complex_mul_exact(arf_t e, arf_t f, const arf_t a, const arf_t b, const arf_t c, const arf_t d);
int _arf_complex_mul_high_approx(arf_t e, arf_t f, const arf_t a, const arf_t b, const arf_t c, const arf_t d, slong bottom);

int arf_complex_sqr(arf_t e, arf_t f, const arf_t a, const arf_t b,
                                      slong prec, arf_rnd_t rnd);

int arf_sum(arf_t s, arf_srcptr terms, slong len, slong prec, arf_rnd_t rnd);

double arf_get_d(const arf_t x, arf_rnd_t rnd);
void arf_set_d(arf_t x, double v);

ARF_INLINE slong
arf_allocated_bytes(const arf_t x)
{
    slong size = fmpz_allocated_bytes(ARF_EXPREF(x));

    if (ARF_HAS_PTR(x))
        size += ARF_PTR_ALLOC(x) * sizeof(ulong);

    return size;
}

int arf_load_str(arf_t res, const char * data);
char * arf_dump_str(const arf_t x);

#ifdef FLINT_HAVE_FILE
int arf_load_file(arf_t res, FILE *stream);
int arf_dump_file(FILE* stream, const arf_t x);
#endif

int arf_dot(arf_t res, const arf_t initial, int subtract,
    arf_srcptr x, slong xstep, arf_srcptr y, slong ystep, slong len, slong prec, arf_rnd_t rnd);

void arf_approx_dot(arf_t res, const arf_t initial, int subtract,
    arf_srcptr x, slong xstep, arf_srcptr y, slong ystep, slong len, slong prec, arf_rnd_t rnd);

#ifdef __cplusplus
}
#endif

#endif
