/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "thread_pool.h"

/*
    if f is continuous with
        f(0) = 0
        f'(x) = alpha for 0 < x < a
        f'(x) = beta for a < x < a + b
    return solution for x to f(x) = (a*alpha + b*beta)*yn/yd
    for 0 <= yn/yd <= 1
*/
ulong _thread_pool_find_work_2(
    ulong a, ulong alpha,
    ulong b, ulong beta,
    ulong yn, ulong yd)
{
    /* very low priority TODO: this can overflow only in very extreme cases */
    ulong y = yn*(a*alpha + b*beta)/yd;

    if (y <= a*alpha)
        return y/alpha;

    return a + (y - a*alpha)/beta;
}

