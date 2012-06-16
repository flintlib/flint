/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2011 Fredrik Johansson

******************************************************************************/

#include <stdio.h>
#include <math.h>
#include "flint.h"
#include "arith.h"
#include "ulong_extras.h"

double arith_bernoulli_number_size(ulong n)
{
    double x;

    /* |B_n| < 2 */
    if (n <= 14)
        return 1.0;

    x = 2 + (n + 1) * log(n + 1) * 1.44269504088897;  /* 1/log(2) */
    x -= n * 4.0941911703612822; /* log2(2*pi*e) */

    return x + 2;
}
