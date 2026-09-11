/*
    Copyright (C) 2026 Maël Hostettler

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef GR_EC_IMPL_H
#define GR_EC_IMPL_H

#include "gr_ec.h"

/*
    Variants of the Jacobian doubling and mixed addition taking a caller
    supplied workspace of GR_EC_JAC_SCRATCH base ring elements. Allocating
    and initializing the temporaries costs a noticeable fraction of an
    operation over small moduli, so the scalar multiplication ladder hoists
    the workspace out of its loop and calls these directly.
*/

#define GR_EC_JAC_SCRATCH 9

WARN_UNUSED_RESULT int _gr_ec_jac_point_dbl_long_weierstrass_ws(gr_ec_jac_point_t res, const gr_ec_jac_point_t P, gr_ptr t, gr_ec_ctx_t ctx);
WARN_UNUSED_RESULT int _gr_ec_jac_point_dbl_short_weierstrass_ws(gr_ec_jac_point_t res, const gr_ec_jac_point_t P, gr_ptr t, gr_ec_ctx_t ctx);
WARN_UNUSED_RESULT int _gr_ec_jac_point_add_aff_point_long_weierstrass_ws(gr_ec_jac_point_t res, const gr_ec_jac_point_t P, const gr_ec_aff_point_t Q, gr_ptr t, gr_ec_ctx_t ctx);
WARN_UNUSED_RESULT int _gr_ec_jac_point_add_aff_point_short_weierstrass_ws(gr_ec_jac_point_t res, const gr_ec_jac_point_t P, const gr_ec_aff_point_t Q, gr_ptr t, gr_ec_ctx_t ctx);

/* Points the context at the method table and element size of a
   representation. Defined in generic.c; called when a context is created. */
void _gr_ec_ctx_init_methods(gr_ec_ctx_t ctx, gr_ec_repr_t repr);

#endif
