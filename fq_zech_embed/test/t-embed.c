/*
    Copyright (C) 2017 Luca De Feo

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_zech_embed.h"
#include "fq_zech_poly.h"


#ifdef T
#undef T
#endif
#ifdef B
#undef B
#endif

#define T fq_zech
#define CAP_T FQ_ZECH
#define B nmod
#include "fq_embed_templates/test/t-embed.c"
#undef B
#undef CAP_T
#undef T
