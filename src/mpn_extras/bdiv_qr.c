/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"

/*
    Hensel division q b == a (mod B^n) for odd b, dispatching to the block
    (classical) algorithm when the divisor is short relative to the
    precision, and to Karp-Markstein otherwise.
*/
void
flint_mpn_bdiv_qr(mp_ptr q, mp_ptr r, mp_srcptr a, mp_size_t an,
    mp_srcptr b, mp_size_t bn, mp_size_t n)
{
    if (n <= 2 || bn < n / 2)
        flint_mpn_bdiv_qr_classical(q, r, a, an, b, bn, n);
    else
        flint_mpn_bdiv_qr_karp_markstein(q, r, a, an, b, bn, n);
}
