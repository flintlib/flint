/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

/* Assumes that the ch are all distinct. More efficient: compute ch from the
   Goepel relation, then check if it is even, etc. */

int
acb_theta_char_is_syzygous(ulong ch1, ulong ch2, ulong ch3, slong g)
{
    ulong n = 1 << (2 * g);
    ulong ch;

    if (ch1 == ch2 || ch2 == ch3 || ch1 == ch3)
    {
        return 0;
    }

    for (ch = 0; ch < n; ch++)
    {
        if (acb_theta_char_is_even(ch, g)
            && (ch != ch1)
            && (ch != ch2)
            && (ch != ch3)
            && acb_theta_char_is_goepel(ch, ch1, ch2, ch3, g))
        {
            return 1;
        }
    }
    return 0;
}
