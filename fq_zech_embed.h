/*
    Copyright (C) 2017 Luca De Feo

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#ifndef FQ_ZECH_EMBED_H
#define FQ_ZECH_EMBED_H

#include "fq_zech.h"

#define T fq_zech
#define CAP_T FQ_ZECH
#include "fq_embed_templates.h"
#undef CAP_T
#undef T

#endif
