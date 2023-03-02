/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

slong
acb_theta_naive_fullprec(const acb_theta_eld_t E, slong prec)
{
    return prec + ceil(ACB_THETA_NAIVE_FULLPREC_ADDLOG
        * n_flog(1 + acb_theta_eld_nb_pts(E), 2));
}
