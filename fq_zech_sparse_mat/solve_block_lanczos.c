/*
    Copyright (C) 2020 Kartik Venkatram

    Algorithm taken from P. Montgomery, "A Block Lanczos Algorithm for 
    Finding Dependencies over GF(2)", Advances in Cryptology - EUROCRYPT '95

    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by th e Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/


#include "fq_zech_sparse_mat.h"

#ifdef T
#undef T
#endif

#define T fq_zech
#define CAP_T FQ_ZECH
#include "fq_sparse_mat_templates/solve_block_lanczos.c"
#undef CAP_T
#undef T
