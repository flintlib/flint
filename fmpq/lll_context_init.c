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

    Copyright (C) 2011 Sebastian Pancratz
    Copyright (C) 2014 Abhinav Baid

******************************************************************************/

#include "fmpq.h"

void
fmpq_lll_context_init(fmpq_lll_t fl, slong delta_num, slong delta_den,
                      slong eta_num, slong eta_den)
{
    fmpq_init(fl->delta);
    fmpq_init(fl->eta);
    fmpq_set_si(fl->delta, delta_num, delta_den);
    fmpq_set_si(fl->eta, eta_num, eta_den);
}
