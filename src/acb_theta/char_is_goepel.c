/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

int
acb_theta_char_is_goepel(ulong ch1, ulong ch2, ulong ch3, ulong ch4, slong g)
{
    if (ch1 == ch2 || ch1 == ch3 || ch1 == ch4
        || ch2 == ch3 || ch2 == ch4 || ch3 == ch4)
    {
        return 0;
    }

    return acb_theta_char_is_even(ch1, g)
        && acb_theta_char_is_even(ch2, g)
        && acb_theta_char_is_even(ch3, g)
        && acb_theta_char_is_even(ch4, g)
        && ((ch1 ^ ch2 ^ ch3 ^ ch4) == 0);
}
