/*
    Copyright (C) 2017 Luca De Feo

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#ifndef FQ_NMOD_EMBED_H
#define FQ_NMOD_EMBED_H

#include "fq_nmod.h"

#define T fq_nmod
#define CAP_T FQ_NMOD
#include "fq_embed_templates.h"
#undef CAP_T
#undef T

#endif
