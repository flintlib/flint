/*
    Copyright (C) 2025 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"

/* 8-bit inverses of the odd bytes, for n_binvert (inline in
   ulong_extras.h): flint_binvert_limb_tab[i] = (2 i + 1)^(-1) mod 2^8 */
const unsigned char flint_binvert_limb_tab[128] = {
  1, 171, 205, 183, 57, 163, 197, 239, 241, 27, 61, 167, 41, 19, 53, 223, 225,
  139, 173, 151, 25, 131, 165, 207, 209, 251, 29, 135, 9, 243, 21, 191, 193,
  107, 141, 119, 249, 99, 133, 175, 177, 219, 253, 103, 233, 211, 245, 159,
  161, 75, 109, 87, 217, 67, 101, 143, 145, 187, 221, 71, 201, 179, 213, 127,
  129, 43, 77, 55, 185, 35, 69, 111, 113, 155, 189, 39, 169, 147, 181, 95, 97,
  11, 45, 23, 153, 3, 37, 79, 81, 123, 157, 7, 137, 115, 149, 63, 65, 235, 13,
  247, 121, 227, 5, 47, 49, 91, 125, 231, 105, 83, 117, 31, 33, 203, 237, 215,
  89, 195, 229, 15, 17, 59, 93, 199, 73, 51, 85, 255
};
