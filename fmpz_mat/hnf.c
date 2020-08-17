/*
    Copyright (C) 2014 Alex J. Best

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mat.h"

void
fmpz_mat_hnf(fmpz_mat_t H, const fmpz_mat_t A)
{
    slong m = A->r, b = fmpz_mat_max_bits(A), cutoff = 2;

    if (b < 0)
        b = -b;

    if (b <= 2)
        cutoff = 52;
    else if (b <= 4)
        cutoff = 38;
    else if (b <= 8)
        cutoff = 30;
    else if (b <= 16)
        cutoff = 11;
    else if (b <= 32)
        cutoff = 11;
    else if (b <= 64)
        cutoff = 5;
    else if (b <= 128)
        cutoff = 4;
    else if (b <= 512)
        cutoff = 3;

    /* 
        TODO: we should call Micciancio-Warisnchi or Pauderis-Storjohann
        when implemented
    */
       
    if (m < cutoff)
        fmpz_mat_hnf_classical(H, A);
    else {
        flint_rand_t state;

        flint_randinit(state);

        fmpz_mat_hnf_pernet_stein(H, A, state);

        flint_randclear(state);
    }
}
