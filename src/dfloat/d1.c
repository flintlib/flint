/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fp_contract.h"
#include <float.h>
#include <string.h>
#include "dfloat.h"
#include "double_extras.h"
#include "machine_vectors.h"
#include "arf.h"
#include "arb.h"
#include "acb.h"
#include "gr.h"
#include "gr_generic.h"

#define N 1
#define DF(name) d1_##name
#define DB(name) d1b_##name
#define DFV(name) _d1_##name
#define DBV(name) _d1b_##name
#define DC(name) d1c_##name
#define DCB(name) d1cb_##name
#define DCV(name) _d1c_##name
#define DCBV(name) _d1cb_##name
#include "template.inc"
#include "ctemplate.inc"
