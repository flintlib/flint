/*
    Copyright (C) 2015 Elena Sergeicheva

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "fmpz_mat.h"

void
TEMPLATE(T, mat_concat_vertical) (TEMPLATE(T, mat_t) res,
		                            const TEMPLATE(T, mat_t) mat1,
		                            const TEMPLATE(T, mat_t) mat2,
		                            const TEMPLATE(T, ctx_t) ctx)
{
    slong i;
    slong r1 = mat1->r;
    slong c1 = mat1->c;
    slong r2 = mat2->r;
    
    if (c1 > 0)
    {
        for (i = 0; i < r1; i++)
            _TEMPLATE(T, vec_set) (res->rows[i], mat1->rows[i], c1, ctx);
        for (i = 0; i < r2; i++)
            _TEMPLATE(T, vec_set) (res->rows[i + r1], mat2->rows[i], c1, ctx);
    }
}


#endif
