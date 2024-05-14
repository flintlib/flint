/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef ARB_VEC_H
#define ARB_VEC_H

#ifdef ARB_VEC_INLINES_C
#define ARB_VEC_INLINE
#else
#define ARB_VEC_INLINE static inline error no inline file
#endif

#include "arb_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/* memory management *********************************************************/

arb_ptr _arb_vec_init(slong n);
void _arb_vec_clear(arb_ptr v, slong n);

arb_ptr _arb_vec_entry_ptr(arb_ptr vec, slong i);

void _arb_vec_swap(arb_ptr res, arb_ptr vec, slong len);

slong _arb_vec_allocated_bytes(arb_srcptr vec, slong len);
double _arb_vec_estimate_allocated_bytes(slong len, slong prec);

void _arb_vec_trim(arb_ptr res, arb_srcptr vec, slong len);

slong _arb_vec_bits(arb_srcptr x, slong len);

/* comparisons ***************************************************************/

int _arb_vec_is_zero(arb_srcptr vec, slong len);
int _arb_vec_is_finite(arb_srcptr x, slong len);

int _arb_vec_equal(arb_srcptr vec1, arb_srcptr vec2, slong len);
int _arb_vec_overlaps(arb_srcptr vec1, arb_srcptr vec2, slong len);
int _arb_vec_contains(arb_srcptr vec1, arb_srcptr vec2, slong len);

/* assignments and conversions ***********************************************/

void _arb_vec_set(arb_ptr res, arb_srcptr vec, slong len);
void _arb_vec_set_round(arb_ptr res, arb_srcptr vec, slong len, slong prec);

void _arb_vec_zero(arb_ptr A, slong n);

int _arb_vec_get_unique_fmpz_vec(fmpz * res,  arb_srcptr vec, slong len);

/* arithmetic ****************************************************************/

void _arb_vec_add(arb_ptr C, arb_srcptr A, arb_srcptr B, slong n, slong prec);
void _arb_vec_sub(arb_ptr C, arb_srcptr A, arb_srcptr B, slong n, slong prec);

void _arb_vec_neg(arb_ptr B, arb_srcptr A, slong n);

/* scalar arithmetic *********************************************************/

void _arb_vec_scalar_mul(arb_ptr res, arb_srcptr vec, slong len, const arb_t c, slong prec);
void _arb_vec_scalar_div(arb_ptr res, arb_srcptr vec, slong len, const arb_t c, slong prec);
void _arb_vec_scalar_mul_fmpz(arb_ptr res, arb_srcptr vec, slong len, const fmpz_t c, slong prec);
void _arb_vec_scalar_mul_2exp_si(arb_ptr res, arb_srcptr src, slong len, slong c);
void _arb_vec_scalar_addmul(arb_ptr res, arb_srcptr vec, slong len, const arb_t c, slong prec);

/* error arithmetic **********************************************************/

void _arb_vec_get_mag(mag_t bound, arb_srcptr vec, slong len);

void _arb_vec_add_error_arf_vec(arb_ptr res, arf_srcptr err, slong len);
void _arb_vec_add_error_mag_vec(arb_ptr res, mag_srcptr err, slong len);

void _arb_vec_indeterminate(arb_ptr vec, slong len);

/* miscellaneous *************************************************************/

void _arb_vec_set_powers(arb_ptr xs, const arb_t x, slong len, slong prec);

/* I/O ***********************************************************************/

void _arb_vec_printn(arb_srcptr vec, slong len, slong ndigits, ulong flags);
void _arb_vec_printd(arb_srcptr vec, slong len, slong ndigits);

#ifdef __cplusplus
}
#endif

#endif
