/*
    Copyright (C) 2013 Mike Hansen
    Copyright (C) 2019 Edouard Rousseau
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "fmpz.h"
#include "fq_nmod.h"
#include "fq_zech.h"
#include "fq_zech_vec.h"
#include "fq_zech_mat.h"
#include "fq_zech_poly.h"
#include "fq_zech_poly_factor.h"

#ifdef T
#undef T
#endif

#define T fq_zech
#define CAP_T FQ_ZECH

#include "fq_poly_factor_templates/clear.c"
#include "fq_poly_factor_templates/concat.c"
#include "fq_poly_factor_templates/factor_berlekamp.c"
#include "fq_poly_factor_templates/factor.c"
#include "fq_poly_factor_templates/factor_cantor_zassenhaus.c"
#include "fq_poly_factor_templates/factor_distinct_deg.c"
#include "fq_poly_factor_templates/factor_equal_deg.c"
#include "fq_poly_factor_templates/factor_equal_deg_prob.c"
#include "fq_poly_factor_templates/factor_kaltofen_shoup.c"
#include "fq_poly_factor_templates/factor_split_single.c"
#include "fq_poly_factor_templates/factor_squarefree.c"
#include "fq_poly_factor_templates/fit_length.c"
#include "fq_poly_factor_templates/init.c"
#include "fq_poly_factor_templates/insert.c"
#include "fq_poly_factor_templates/is_irreducible_ben_or.c"
#include "fq_poly_factor_templates/is_irreducible.c"
#include "fq_poly_factor_templates/is_irreducible_ddf.c"
#include "fq_poly_factor_templates/is_squarefree.c"
#include "fq_poly_factor_templates/iterated_frobenius_preinv.c"
#include "fq_poly_factor_templates/pow.c"
#include "fq_poly_factor_templates/print.c"
#include "fq_poly_factor_templates/print_pretty.c"
#include "fq_poly_factor_templates/realloc.c"
#include "fq_poly_factor_templates/roots.c"
#include "fq_poly_factor_templates/set.c"

#undef CAP_T
#undef T
