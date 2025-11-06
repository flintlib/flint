/*
    Copyright (C) 2025 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef N_FFT_IMPL_H
#define N_FFT_IMPL_H

#include "longlong.h"  /* provides flint_clz */
#include "n_fft.h"

void dft_node_lazy_4_4(nn_ptr p, ulong depth, ulong node, n_fft_args_t F);
void dft_lazy_2_4(nn_ptr p, ulong depth, n_fft_args_t F);
void dft_lazy_1_4(nn_ptr p, ulong depth, n_fft_args_t F);

void idft_node_lazy_1_2(nn_ptr p, ulong depth, ulong node, n_fft_args_t F);
void idft_lazy_1_4(nn_ptr p, ulong depth, n_fft_args_t F);

void tft_node_lazy_4_4(nn_ptr p, ulong olen, ulong depth, ulong node, n_fft_args_t F);
void tft_lazy_1_4(nn_ptr p, ulong ilen, ulong olen, n_fft_args_t F);


/* exponent of next power of 2 for x > 2 */
FLINT_FORCE_INLINE
ulong n_clog2_gt2(ulong x)
{
    return FLINT_BITS - flint_clz(x - 1);
}

#endif  /* N_FFT_IMPL_H */

