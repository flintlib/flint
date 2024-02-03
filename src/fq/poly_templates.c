/*
    Copyright (C) 2007 David Howden
    Copyright (C) 2007-2012 William Hart
    Copyright (C) 2008 Richard Howell-Peak
    Copyright (C) 2010-2012 Sebastian Pancratz
    Copyright (C) 2011, 2012 Fredrik Johansson
    Copyright (C) 2012 Andres Goens
    Copyright (C) 2013 Martin Lee
    Copyright (C) 2013 Mike Hansen
    Copyright (C) 2017 Luca De Feo

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fq.h"
#include "fq_vec.h"
#include "fq_mat.h"
#include "fq_poly.h"

#ifdef T
#undef T
#endif

#define T fq
#define CAP_T FQ

#include "fq_poly_templates/add.c"
#include "fq_poly_templates/add_series.c"
#include "fq_poly_templates/add_si.c"
#include "fq_poly_templates/clear.c"
#include "fq_poly_templates/comparisons.c"
#include "fq_poly_templates/compose.c"
#include "fq_poly_templates/compose_mod.c"
#include "fq_poly_templates/compose_mod_brent_kung.c"
#include "fq_poly_templates/compose_mod_brent_kung_precomp_preinv.c"
#include "fq_poly_templates/compose_mod_brent_kung_preinv.c"
#include "fq_poly_templates/compose_mod_horner.c"
#include "fq_poly_templates/compose_mod_horner_preinv.c"
#include "fq_poly_templates/compose_mod_preinv.c"
#include "fq_poly_templates/deflate.c"
#include "fq_poly_templates/deflation.c"
#include "fq_poly_templates/derivative.c"
#include "fq_poly_templates/div.c"
#include "fq_poly_templates/div_newton_n_preinv.c"
#include "fq_poly_templates/div_series.c"
#include "fq_poly_templates/divides.c"
#include "fq_poly_templates/divrem.c"
#include "fq_poly_templates/divrem_f.c"
#include "fq_poly_templates/divrem_newton_n_preinv.c"
#include "fq_poly_templates/equal.c"
#include "fq_poly_templates/equal_trunc.c"
#include "fq_poly_templates/evaluate_fq.c"
#include "fq_poly_templates/evaluate_fq_vec.c"
#include "fq_poly_templates/evaluate_fq_vec_fast.c"
#include "fq_poly_templates/evaluate_fq_vec_iter.c"
#include "fq_poly_templates/fit_length.c"
#include "fq_poly_templates/gcd.c"
#include "fq_poly_templates/gcd_euclidean_f.c"
#include "fq_poly_templates/gen.c"
#include "fq_poly_templates/get_coeff.c"
#include "fq_poly_templates/get_str.c"
#include "fq_poly_templates/get_str_pretty.c"
#include "fq_poly_templates/hamming_weight.c"
#include "fq_poly_templates/inflate.c"
#include "fq_poly_templates/init.c"
#include "fq_poly_templates/inv_series_newton.c"
#include "fq_poly_templates/invsqrt_series.c"
#include "fq_poly_templates/io.c"
#include "fq_poly_templates/is_gen.c"
#include "fq_poly_templates/make_monic.c"
#define USE_MUL_REORDER 1
#include "fq_poly_templates/mul.c"
#undef USE_MUL_REORDER
/* #include "fq_poly_templates/mul_classical.c" */
#include "fq_poly_templates/mul_KS.c"
#include "fq_poly_templates/mul_reorder.c"
#include "fq_poly_templates/mulhigh.c"
#include "fq_poly_templates/mulhigh_classical.c"
#include "fq_poly_templates/mullow.c"
#include "fq_poly_templates/mullow_classical.c"
#include "fq_poly_templates/mullow_KS.c"
#include "fq_poly_templates/mulmod.c"
#include "fq_poly_templates/mulmod_preinv.c"
#include "fq_poly_templates/neg.c"
#include "fq_poly_templates/normalise.c"
#include "fq_poly_templates/one.c"
#include "fq_poly_templates/pow.c"
#include "fq_poly_templates/pow_trunc.c"
#include "fq_poly_templates/pow_trunc_binexp.c"
#include "fq_poly_templates/powmod_fmpz_binexp.c"
#include "fq_poly_templates/powmod_fmpz_binexp_preinv.c"
#include "fq_poly_templates/powmod_fmpz_sliding_preinv.c"
#include "fq_poly_templates/powmod_ui_binexp.c"
#include "fq_poly_templates/powmod_ui_binexp_preinv.c"
#include "fq_poly_templates/powmod_x_fmpz_preinv.c"
#include "fq_poly_templates/randtest.c"
#include "fq_poly_templates/randtest_irreducible.c"
#include "fq_poly_templates/randtest_monic.c"
#include "fq_poly_templates/realloc.c"
#include "fq_poly_templates/rem.c"
#include "fq_poly_templates/remove.c"
#include "fq_poly_templates/reverse.c"
#include "fq_poly_templates/scalar_addmul_fq.c"
#include "fq_poly_templates/scalar_div_fq.c"
#include "fq_poly_templates/scalar_mul_fq.c"
#include "fq_poly_templates/scalar_submul_fq.c"
#include "fq_poly_templates/set.c"
#include "fq_poly_templates/set_coeff.c"
#include "fq_poly_templates/set_fmpz_mod_poly.c"
#include "fq_poly_templates/set_fq.c"
#include "fq_poly_templates/set_length.c"
#include "fq_poly_templates/set_nmod_poly.c"
#include "fq_poly_templates/set_trunc.c"
#include "fq_poly_templates/shift_left.c"
#include "fq_poly_templates/shift_right.c"
#define USE_SQR_REORDER 1
#include "fq_poly_templates/sqr.c"
#undef USE_SQR_REORDER
#include "fq_poly_templates/sqr_classical.c"
#include "fq_poly_templates/sqr_KS.c"
#include "fq_poly_templates/sqr_reorder.c"
#include "fq_poly_templates/sqrt.c"
#include "fq_poly_templates/sqrt_series.c"
#include "fq_poly_templates/sub.c"
#include "fq_poly_templates/sub_series.c"
#include "fq_poly_templates/swap.c"
#include "fq_poly_templates/tree.c"
#include "fq_poly_templates/truncate.c"
#include "fq_poly_templates/xgcd.c"
#include "fq_poly_templates/xgcd_euclidean_f.c"
#include "fq_poly_templates/zero.c"

#undef CAP_T
#undef T
