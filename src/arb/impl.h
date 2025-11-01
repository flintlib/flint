/*
    Copyright (C) 2025 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef ARB_IMPL_H
#define ARB_IMPL_H

#include "acb_types.h"

ulong euler_mod_p_powsum_noredc(ulong n, ulong p, const unsigned int * divtab);
ulong euler_mod_p_powsum(ulong n, ulong p, const unsigned int * divtab);
ulong euler_mod_p_powsum_1(ulong n, ulong p);
void divisor_table_odd(unsigned int * tab, slong len);
void arb_gamma_stirling_eval(arb_t s, const arb_t z, slong nterms, int digamma, slong prec);
int _arb_log_ui_smooth(arb_t res, ulong n, slong prec);
void arb_gamma_stirling_coeff(arb_t b, ulong k, int digamma, slong prec);
void _arb_fmpz_divapprox_newton(fmpz_t res, const fmpz_t x, const fmpz_t y, slong exp);
void _arb_dot_addmul_generic(nn_ptr sum, nn_ptr serr, nn_ptr tmp, slong sn, nn_srcptr xptr, slong xn, nn_srcptr yptr, slong yn, int negative, flint_bitcnt_t shift);
void _arb_dot_add_generic(nn_ptr sum, nn_ptr serr, nn_ptr tmp, slong sn, nn_srcptr xptr, slong xn, int negative, flint_bitcnt_t shift);
void mag_borwein_error(mag_t err, slong n);
void arb_exp_arf_huge(arb_t z, const arf_t x, slong mag, slong prec, int minus_one);
slong _arb_compute_bs_exponents(slong * tab, slong n);
slong _arb_get_exp_pos(const slong * tab, slong step);
void mag_set_ui_2exp_small(mag_t z, ulong x, slong e);
void mag_agm(mag_t res, const mag_t x, const mag_t y);
void arb_gamma_stirling_bound(mag_ptr err, const arb_t x, slong k0, slong knum, slong n);
void acb_gamma_stirling_bound(mag_ptr err, const acb_t z, slong k0, slong knum, slong n);
double d_lambertw_branch1(double x);
void arb_exp_taylor_sum_rs_generic(arb_t res, const arb_t x, slong N, slong prec);
void arb_sin_cos_fmpz_div_2exp_bsplit(arb_t wsin, arb_t wcos, const fmpz_t x, flint_bitcnt_t r, slong prec);
void arb_sin_cos_taylor_sum_rs(arb_t s, const arb_t x, slong N, int cosine, slong prec);

#endif
