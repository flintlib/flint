/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

slong
sp2gz_nb_fundamental(slong g)
{
    if (g == 1)
        return 1;
    if (g == 2)
        return 19;
    else
        return 19 * ((g * (g - 1)) / 2) + (1 << g);
}
