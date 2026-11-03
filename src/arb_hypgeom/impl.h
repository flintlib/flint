/*
    Copyright (C) 2025 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef ARB_HYPGEOM_IMPL_H
#define ARB_HYPGEOM_IMPL_H

#include "arb_types.h"

void arb_hypgeom_gamma_stirling_choose_param(int * reflect, slong * r, slong * n, const arb_t x, int use_reflect, int digamma, slong prec);
void arb_hypgeom_gamma_stirling_inner(arb_t s, const arb_t z, slong N, slong prec);
int arb_hypgeom_gamma_exact(arb_t res, const arb_t x, int reciprocal, slong prec);
double arf_get_d_log2_abs_approx_clamped(const arf_t x);
int _arf_increment_fast(arf_t x, slong prec);
void _arb_increment_fast(arb_t x, slong prec);

/* roots of Legendre polynomials */
#define ARB_HYPGEOM_LEGENDRE_ROOT_MIN_N 64
/* n(n+1) and (j+1)(j+2), j <= n + 2, must fit in a slong */
#if FLINT_BITS == 64
#define ARB_HYPGEOM_LEGENDRE_ROOT_MAX_N (UWORD(1) << 28)
#else
#define ARB_HYPGEOM_LEGENDRE_ROOT_MAX_N (UWORD(1) << 14)
#endif
#define ARB_HYPGEOM_LEGENDRE_ROOT_MIN_PREC 128

void _arb_hypgeom_legendre_p_ui_root_petras(arb_t res, mag_t errc, mag_t errd, ulong n, ulong k, slong prec);
int _arb_hypgeom_legendre_p_ui_root_isolate(arb_t res, ulong n, ulong k, slong prec);
int _arb_hypgeom_legendre_p_ui_taylor_tail(mag_t T0, mag_t T1, ulong n, const arb_t m, const mag_t r, const mag_t tau, slong R);
void _arb_hypgeom_legendre_p_ui_disk_radius(mag_t r, const arb_t m, ulong n);
int _arb_hypgeom_legendre_p_ui_root_bitburst_param(arb_t res, arb_t weight, ulong n, ulong k, slong prec, slong limb_max_prec, double transfer_blocks_min, slong fft_bits);
void _arb_hypgeom_legendre_p_ui_root_enclosure(arb_t res, arb_t weight, ulong n, ulong k, slong prec);

#endif
