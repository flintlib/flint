/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Include functions *********************************************************/

#include "t-algebraic.c"
#include "t-api.c"
#include "t-map.c"
#include "t-lazy.c"
#include "t-lazy_api.c"
#include "t-flat.c"
#include "t-transcendental.c"
#include "t-richardson.c"
#include "t-elementary.c"
#include "t-roots.c"
#include "t-lazy_qqbar.c"
#include "t-lazy_numeric.c"
#include "t-modular.c"
#include "t-no_roots.c"
#include "t-gc.c"
#include "t-dense.c"
#include "t-dense_form.c"
#include "t-threads.c"
#include "t-subfields.c"
#include "t-poly.c"
#include "t-factor.c"
#include "t-trig.c"
#include "t-special.c"
#include "t-catalog.c"
#include "t-options.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(gr_tower_algebraic),
    TEST_FUNCTION(gr_tower_api),
    TEST_FUNCTION(gr_tower_map),
    TEST_FUNCTION(gr_tower_lazy),
    TEST_FUNCTION(gr_tower_lazy_api),
    TEST_FUNCTION(gr_tower_flat),
    TEST_FUNCTION(gr_tower_transcendental),
    TEST_FUNCTION(gr_tower_richardson),
    TEST_FUNCTION(gr_tower_elementary),
    TEST_FUNCTION(gr_tower_roots),
    TEST_FUNCTION(gr_tower_lazy_qqbar),
    TEST_FUNCTION(gr_tower_lazy_numeric),
    TEST_FUNCTION(gr_tower_modular),
    TEST_FUNCTION(gr_tower_no_roots),
    TEST_FUNCTION(gr_tower_gc),
    TEST_FUNCTION(gr_tower_dense),
    TEST_FUNCTION(gr_tower_dense_form),
    TEST_FUNCTION(gr_tower_threads),
    TEST_FUNCTION(gr_tower_subfields),
    TEST_FUNCTION(gr_tower_poly),
    TEST_FUNCTION(gr_tower_factor),
    TEST_FUNCTION(gr_tower_trig),
    TEST_FUNCTION(gr_tower_special),
    TEST_FUNCTION(gr_tower_catalog),
    TEST_FUNCTION(gr_tower_options)
};

/* main function *************************************************************/

TEST_MAIN(tests)
