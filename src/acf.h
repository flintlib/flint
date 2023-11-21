/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef ACF_H
#define ACF_H

#ifdef ACF_INLINES_C
#define ACF_INLINE
#else
#define ACF_INLINE static inline
#endif

#include "arf.h"
#include "acf_types.h"

#ifdef __cplusplus
extern "C" {
#endif

#define acf_realref(x) (&(x)->real)
#define acf_imagref(x) (&(x)->imag)

ACF_INLINE void
acf_init(acf_t x)
{
    arf_init(acf_realref(x));
    arf_init(acf_imagref(x));
}

ACF_INLINE void
acf_clear(acf_t x)
{
    arf_clear(acf_realref(x));
    arf_clear(acf_imagref(x));
}

ACF_INLINE acf_ptr _acf_vec_init(slong n) { return (acf_ptr) _arf_vec_init(2 * n); }
ACF_INLINE void _acf_vec_clear(acf_ptr v, slong n) { _arf_vec_clear((arf_ptr) v, 2 * n); }

ACF_INLINE arf_ptr acf_real_ptr(acf_t z) { return acf_realref(z); }
ACF_INLINE arf_ptr acf_imag_ptr(acf_t z) { return acf_imagref(z); }

ACF_INLINE void
acf_set(acf_t z, const acf_t x)
{
    arf_set(acf_realref(z), acf_realref(x));
    arf_set(acf_imagref(z), acf_imagref(x));
}

ACF_INLINE void
acf_swap(acf_t z, acf_t x)
{
    arf_swap(acf_realref(z), acf_realref(x));
    arf_swap(acf_imagref(z), acf_imagref(x));
}

ACF_INLINE int
acf_equal(const acf_t x, const acf_t y)
{
    return arf_equal(acf_realref(x), acf_realref(y)) &&
           arf_equal(acf_imagref(x), acf_imagref(y));
}

/* todo: document */
ACF_INLINE void
acf_printd(const acf_t x, slong n)
{
    arf_printd(acf_realref(x), n);
    flint_printf(" + ");
    arf_printd(acf_imagref(x), n);
    flint_printf("*I");
}

/* todo: document */
ACF_INLINE slong
acf_bits(const acf_t x)
{
    slong b1, b2;
    b1 = arf_bits(acf_realref(x));
    b2 = arf_bits(acf_imagref(x));
    return FLINT_MAX(b1, b2);
}

ACF_INLINE slong
acf_allocated_bytes(const acf_t x)
{
    return arf_allocated_bytes(acf_realref(x)) + arf_allocated_bytes(acf_imagref(x));
}

/* todo: document */
ACF_INLINE void acf_randtest(acf_t x, flint_rand_t state, slong bits, slong mag_bits)
{
    arf_randtest(acf_realref(x), state, bits, mag_bits);
    arf_randtest(acf_imagref(x), state, bits, mag_bits);
}

/* todo: document */
ACF_INLINE void
acf_get_mag(mag_t res, const acf_t x)
{
    mag_t t, u;
    mag_init(t);
    mag_init(u);
    arf_get_mag(t, acf_realref(x));
    arf_get_mag(u, acf_imagref(x));
    mag_hypot(res, t, u);
    mag_clear(t);
    mag_clear(u);
}

/* todo: document */
ACF_INLINE void
acf_neg(acf_t z, const acf_t x)
{
    arf_neg(acf_realref(z), acf_realref(x));
    arf_neg(acf_imagref(z), acf_imagref(x));
}

/* todo: document */
ACF_INLINE int
acf_set_round(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)
{
    return arf_set_round(acf_realref(res), acf_realref(x), prec, rnd) |
           (arf_set_round(acf_imagref(res), acf_imagref(x), prec, rnd) << 1);
}

/* todo: document */
ACF_INLINE int
acf_neg_round(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)
{
    return arf_neg_round(acf_realref(res), acf_realref(x), prec, rnd) |
           (arf_neg_round(acf_imagref(res), acf_imagref(x), prec, rnd) << 1);
}

ACF_INLINE int
acf_add(acf_t res, const acf_t x, const acf_t y, slong prec, arf_rnd_t rnd)
{
    return arf_add(acf_realref(res), acf_realref(x), acf_realref(y), prec, rnd) |
           (arf_add(acf_imagref(res), acf_imagref(x), acf_imagref(y), prec, rnd) << 1);
}

ACF_INLINE int
acf_sub(acf_t res, const acf_t x, const acf_t y, slong prec, arf_rnd_t rnd)
{
    return arf_sub(acf_realref(res), acf_realref(x), acf_realref(y), prec, rnd) |
           (arf_sub(acf_imagref(res), acf_imagref(x), acf_imagref(y), prec, rnd) << 1);
}

ACF_INLINE int
acf_mul(acf_t res, const acf_t x, const acf_t y, slong prec, arf_rnd_t rnd)
{
    if (x == y)
        return arf_complex_sqr(acf_realref(res), acf_imagref(res),
            acf_realref(x), acf_imagref(x),
            prec, rnd);
    else
        return arf_complex_mul(acf_realref(res), acf_imagref(res),
            acf_realref(x), acf_imagref(x),
            acf_realref(y), acf_imagref(y),
            prec, rnd);
}

void acf_approx_inv(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd);
void acf_approx_div(acf_t res, const acf_t x, const acf_t y, slong prec, arf_rnd_t rnd);
void acf_approx_sqrt(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd);

void acf_approx_dot(acf_t res, const acf_t initial, int subtract, acf_srcptr x, slong xstep, acf_srcptr y, slong ystep, slong len, slong prec, arf_rnd_t rnd);

#ifdef __cplusplus
}
#endif

#endif
