/*
    Copyright (C) 2017 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_calc.h"

void
acb_calc_integrate_opt_init(acb_calc_integrate_opt_t options)
{
    options->deg_limit = 0;
    options->eval_limit = 0;
    options->depth_limit = 0;
    options->use_heap = 0;
    options->verbose = 0;
}

