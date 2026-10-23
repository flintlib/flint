/*
    Copyright (C) 2016 Pascal Molin
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef ACB_DFT_H
#define ACB_DFT_H

#include "acb_types.h"
#include "gr_dft.h"

#ifdef __cplusplus
extern "C" {
#endif

/* This module is a thin wrapper around gr_dft, which computes DFTs
   over generic rings and provides the complex-ball functionality
   through fixed-point arithmetic with rigorous error bounds
   (gr_dft_acb), falling back to ball arithmetic when fixed point
   does not apply. Product DFTs, convolutions and access to specific
   algorithms are available in gr_dft. */

typedef gr_dft_acb_pre_struct acb_dft_pre_struct;
typedef acb_dft_pre_struct acb_dft_pre_t[1];

void acb_dft(acb_ptr w, acb_srcptr v, slong len, slong prec);
void acb_dft_inverse(acb_ptr w, acb_srcptr v, slong len, slong prec);

void acb_dft_precomp_init(acb_dft_pre_t pre, slong len, slong prec);
void acb_dft_precomp_clear(acb_dft_pre_t pre);
void acb_dft_precomp(acb_ptr w, acb_srcptr v, const acb_dft_pre_t pre, slong prec);
void acb_dft_inverse_precomp(acb_ptr w, acb_srcptr v, const acb_dft_pre_t pre, slong prec);

#ifdef __cplusplus
}
#endif

#endif
