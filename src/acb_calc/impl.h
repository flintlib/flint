/*
    Copyright (C) 2026 Fredrik Johansson
    Developed using Claude Fable 5.2

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef ACB_CALC_IMPL_H
#define ACB_CALC_IMPL_H

#include "arb_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/* The cache of Gauss-Legendre rules (per thread) used by
   acb_calc_integrate, and shared with other quadrature code: the rules
   with _acb_calc_gl_steps[i] nodes, 0 <= i < ACB_CALC_GL_STEPS, in
   steps of about 2^(1/2) (1, 2, 4, 6, 8, 12, 16, 22, 32, ...). */
#define ACB_CALC_GL_STEPS 38

extern const slong _acb_calc_gl_steps[ACB_CALC_GL_STEPS];

/* Sets *x and *w to the (n + 1) / 2 nonnegative nodes of the rule with
   n = _acb_calc_gl_steps[i] nodes, in decreasing order, and to their
   weights, with a precision of at least prec bits (the entries can have
   a higher precision and should be rounded by the user).  The pointers
   remain valid until the next call for the same i with a higher
   precision, or until flint_cleanup. */
void _acb_calc_gl_nodes(arb_srcptr * x, arb_srcptr * w, slong i, slong prec);

#ifdef __cplusplus
}
#endif

#endif
