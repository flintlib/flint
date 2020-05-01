/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#ifndef CA_QQBAR_H
#define CA_QQBAR_H

#ifdef CA_QQBAR_INLINES_C
#define CA_QQBAR_INLINE
#else
#define CA_QQBAR_INLINE static __inline__
#endif

#ifdef __cplusplus
extern "C" {
#endif

#include "flint/fmpz_poly.h"
#include "flint/fmpq_poly.h"
#include "acb.h"

typedef struct
{
    fmpz_poly_struct poly;
    acb_struct enclosure;
}
ca_qqbar_struct;

#define CA_QQBAR_POLY(x) (&((x)->poly))
#define CA_QQBAR_COEFFS(x) ((&((x)->poly))->coeffs)
#define CA_QQBAR_ENCLOSURE(x) (&((x)->enclosure))

#define CA_QQBAR_DEFAULT_PREC 128

typedef ca_qqbar_struct ca_qqbar_t[1];

/* Memory management */

void ca_qqbar_init(ca_qqbar_t res);

void ca_qqbar_clear(ca_qqbar_t res);

/* Assignment */

void ca_qqbar_swap(ca_qqbar_t x, ca_qqbar_t y);

void ca_qqbar_set(ca_qqbar_t res, const ca_qqbar_t x);

void ca_qqbar_set_si(ca_qqbar_t res, slong x);

void ca_qqbar_set_ui(ca_qqbar_t res, ulong x);

void ca_qqbar_set_fmpz(ca_qqbar_t res, const fmpz_t x);

void ca_qqbar_set_fmpq(ca_qqbar_t res, const fmpq_t x);

/* Properties */

CA_QQBAR_INLINE slong
ca_qqbar_degree(const ca_qqbar_t x)
{
    return fmpz_poly_degree(CA_QQBAR_POLY(x));
}

CA_QQBAR_INLINE int
ca_qqbar_is_rational(const ca_qqbar_t x)
{
    return (ca_qqbar_degree(x) == 1);
}

CA_QQBAR_INLINE int
ca_qqbar_is_integer(const ca_qqbar_t x)
{
    return ca_qqbar_is_rational(x) && fmpz_is_one(CA_QQBAR_COEFFS(x) + 1);
}

CA_QQBAR_INLINE int
ca_qqbar_is_algebraic_integer(const ca_qqbar_t x)
{
    return fmpz_is_one(CA_QQBAR_COEFFS(x) + ca_qqbar_degree(x));
}

CA_QQBAR_INLINE int
ca_qqbar_is_zero(const ca_qqbar_t x)
{
    return ca_qqbar_is_integer(x) && fmpz_is_zero(CA_QQBAR_COEFFS(x));
}

CA_QQBAR_INLINE int
ca_qqbar_is_one(const ca_qqbar_t x)
{
    return ca_qqbar_is_integer(x) && (fmpz_equal_si(CA_QQBAR_COEFFS(x), -1));
}

CA_QQBAR_INLINE int
ca_qqbar_is_neg_one(const ca_qqbar_t x)
{
    return ca_qqbar_is_integer(x) && fmpz_is_one(CA_QQBAR_COEFFS(x));
}

int ca_qqbar_real_sgn(const ca_qqbar_t x);

int ca_qqbar_imag_sgn(const ca_qqbar_t x);

CA_QQBAR_INLINE int
ca_qqbar_is_real(const ca_qqbar_t x)
{
    return ca_qqbar_imag_sgn(x) == 0;
}

/* Special values */

CA_QQBAR_INLINE void
ca_qqbar_zero(ca_qqbar_t res)
{
    ca_qqbar_set_ui(res, 0);
}

CA_QQBAR_INLINE void
ca_qqbar_one(ca_qqbar_t res)
{
    ca_qqbar_set_ui(res, 1);
}

void ca_qqbar_i(ca_qqbar_t res);

void ca_qqbar_phi(ca_qqbar_t res);

/* Random generation */

void ca_qqbar_randtest(ca_qqbar_t res, flint_rand_t state, slong deg, slong bits);

void ca_qqbar_randtest_real(ca_qqbar_t res, flint_rand_t state, slong deg, slong bits);

void ca_qqbar_randtest_nonreal(ca_qqbar_t res, flint_rand_t state, slong deg, slong bits);

/* Input and output */

void ca_qqbar_print(const ca_qqbar_t x);

/* Comparisons */

int ca_qqbar_equal(const ca_qqbar_t x, const ca_qqbar_t y);

int ca_qqbar_cmp_re(const ca_qqbar_t x, const ca_qqbar_t y);

int ca_qqbar_cmp_im(const ca_qqbar_t x, const ca_qqbar_t y);

int ca_qqbar_cmpabs_re(const ca_qqbar_t x, const ca_qqbar_t y);

int ca_qqbar_cmpabs_im(const ca_qqbar_t x, const ca_qqbar_t y);

int ca_qqbar_cmpabs(const ca_qqbar_t x, const ca_qqbar_t y);

/* Complex parts */

void ca_qqbar_conj(ca_qqbar_t res, const ca_qqbar_t x);

void ca_qqbar_re(ca_qqbar_t res, const ca_qqbar_t x);

void ca_qqbar_im(ca_qqbar_t res, const ca_qqbar_t x);

void ca_qqbar_re_im(ca_qqbar_t res1, ca_qqbar_t res2, const ca_qqbar_t x);

void ca_qqbar_abs(ca_qqbar_t res, const ca_qqbar_t x);

void ca_qqbar_abs2(ca_qqbar_t res, const ca_qqbar_t x);

void ca_qqbar_sgn(ca_qqbar_t res, const ca_qqbar_t x);

int ca_qqbar_csgn(const ca_qqbar_t x);

/* Integer parts */

void ca_qqbar_floor(fmpz_t res, const ca_qqbar_t x);

void ca_qqbar_ceil(fmpz_t res, const ca_qqbar_t x);

/* Arithmetic */

void ca_qqbar_neg(ca_qqbar_t res, const ca_qqbar_t x);

void ca_qqbar_add(ca_qqbar_t res, const ca_qqbar_t x, const ca_qqbar_t y);

void ca_qqbar_sub(ca_qqbar_t res, const ca_qqbar_t x, const ca_qqbar_t y);

void ca_qqbar_mul(ca_qqbar_t res, const ca_qqbar_t x, const ca_qqbar_t y);

void ca_qqbar_div(ca_qqbar_t res, const ca_qqbar_t x, const ca_qqbar_t y);

CA_QQBAR_INLINE void
ca_qqbar_sqr(ca_qqbar_t res, const ca_qqbar_t x)
{
    ca_qqbar_mul(res, x, x);
}

void ca_qqbar_inv(ca_qqbar_t res, const ca_qqbar_t x);

void ca_qqbar_mul_2exp_si(ca_qqbar_t res, const ca_qqbar_t x, slong exp);

void ca_qqbar_pow_ui(ca_qqbar_t res, const ca_qqbar_t x, ulong n);

void ca_qqbar_root_ui(ca_qqbar_t res, const ca_qqbar_t x, ulong n);

CA_QQBAR_INLINE void
ca_qqbar_sqrt(ca_qqbar_t res, const ca_qqbar_t x)
{
    ca_qqbar_root_ui(res, x, 2);
}

CA_QQBAR_INLINE void
ca_qqbar_rsqrt(ca_qqbar_t res, const ca_qqbar_t x)
{
    ca_qqbar_sqrt(res, x);
    ca_qqbar_inv(res, res);
}

/* Internal functions */

void ca_qqbar_scalar_op(ca_qqbar_t res, const ca_qqbar_t x, const fmpz_t a, const fmpz_t b, const fmpz_t c);

void ca_qqbar_fmpz_poly_composed_op(fmpz_poly_t res, const fmpz_poly_t A, const fmpz_poly_t B, int op);

void ca_qqbar_binary_op(ca_qqbar_t res, const ca_qqbar_t x, const ca_qqbar_t y, int op);

int _ca_qqbar_validate_enclosure(acb_t res, const fmpz_poly_t poly, const acb_t z, slong max_prec);

void _ca_qqbar_enclosure_raw(acb_t res, const fmpz_poly_t poly, const acb_t zin, slong prec);

void ca_qqbar_enclosure_raw(acb_t res, const ca_qqbar_t x, slong prec);

#ifdef __cplusplus
}
#endif

#endif

