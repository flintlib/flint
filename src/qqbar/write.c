/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "calcium.h"
#include "qqbar.h"

/* FIXME: Should this function be used somewhere? */

/* FIXME: Remove this guard against warnings. Best thing would probably be to
 * implement an *-impl.h to keep track of local functions. */
#ifdef __GNUC__
# pragma GCC diagnostic ignored "-Wmissing-prototypes"
# pragma message "qqbar_writen is currently unused/untested/undocumented!"
# pragma message "qqbar_writend is currently unused/untested/undocumented!"
# pragma message "qqbar_get_str_nd is currently unused/untested/undocumented!"
#endif

/* Todo: document, unify and generalize. This is currently only
   used for the Python interface, but other wrappers will want
   good printing code too. */

void
qqbar_writen(calcium_stream_t out, const qqbar_t x, slong n)
{
    acb_t t;
    slong prec;

    n = FLINT_MAX(1, n);
    prec = n * 3.333 + 10;

    acb_init(t);
    qqbar_get_acb(t, x, prec);

    calcium_write_acb(out, t, n, ARB_STR_NO_RADIUS);
    acb_clear(t);
}

void
qqbar_writend(calcium_stream_t out, const qqbar_t x, slong n)
{
    qqbar_writen(out, x, n);
    calcium_write(out, " (deg ");
    calcium_write_si(out, qqbar_degree(x));
    calcium_write(out, ")");
}

char * qqbar_get_str_nd(const qqbar_t x, slong n)
{
    calcium_stream_t out;
    calcium_stream_init_str(out);
    qqbar_writend(out, x, n);
    return out->s;
}
