/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "radix.h"

/*
    Hensel division q*b == a (mod B^n), dispatching to whichever low-level
    routine is faster for the shape of the inputs.

    The selection is intentionally simple for now.

    Note: "classical" is actually based on block multiplication, so it is
    fine asymptotically.
*/
int
radix_divmod_bn(nn_ptr q, nn_ptr rem, nn_srcptr a, slong an,
    nn_srcptr b, slong bn, slong n, const radix_t radix)
{
    if (n <= 2 || bn < n / 2)
        return radix_divmod_bn_classical(q, rem, a, an, b, bn, n, radix);
    else
        return radix_divmod_bn_karp_markstein(q, rem, a, an, b, bn, n, radix);
}
