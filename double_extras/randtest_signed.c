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

   Copyright (C) 2012 Fredrik Johansson
   Copyright (C) 2014 Abhinav Baid

******************************************************************************/

#include "double_extras.h"
#include "ulong_extras.h"

double
d_randtest_signed(flint_rand_t state, slong minexp, slong maxexp)
{
    double d, t;
    slong exp, kind;
    d = d_randtest(state);
    exp = minexp + n_randint(state, maxexp - minexp + 1);
    t = ldexp(d, exp);
    kind = n_randint(state, 3);
	if (kind == 2)
        return t;
    else if (kind == 1)
        return -t;
    else
        return 0;
}
