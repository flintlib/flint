/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2013 Andres Goens
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq.h"

#ifdef T
#undef T
#endif

#define T fq
#define CAP_T FQ
#include "fq_templates/test/t-mul_ui.c"
#undef CAP_T
#undef T
