/*
    Copyright (C) 2015 William Hart
    Copyright (C) 2015 Vladimir Glazachev
    Copyright (C) 2024 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"

/* Returns a_precomp = floor(a * 2**FLINT_BITS / n) */
ulong
n_mulmod_precomp_shoup(ulong a, ulong n)
{
   ulong a_precomp, r;
   udiv_qrnnd(a_precomp, r, a, UWORD(0), n);  // requires a < n
   return a_precomp;
}

/*-------------------------------------------------------------*/
/* notes on Shoup's modular multiplication with precomputation */
/*-------------------------------------------------------------*/

// Sources:
// - NTL code (in particular file FFT.cpp, consulted in NTL v11.5.1)
// - Victor Shoup, "Arithmetic Software Libraries", 2021 (https://doi.org/10.1017/9781108854207.012)
//   (chapter 9 of the book https://doi.org/10.1017/9781108854207 )

// below, B == FLINT_BITS

// input: w, t, p
// output: (w*t) mod p
// constraints:
//    (Cn) n has <= B-1 bits (n < 2**(B-1))
//    (Ca) n has <= B-1 bits (n < 2**(B-1))
