/*
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#ifdef T

#include "templates.h"

FLINT_DLL void TEMPLATE(T, gcdinv)(TEMPLATE(T, t) rop, TEMPLATE(T, t) inv,
                    const TEMPLATE(T, t) op,
                    const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL int TEMPLATE(T, is_invertible)(const TEMPLATE(T, t) op,
                           const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL int TEMPLATE(T, is_invertible_f)(TEMPLATE(T, t) rop, const TEMPLATE(T, t) op,
                             const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, div)(TEMPLATE(T, t) rop, const TEMPLATE(T, t) op1,
                 const TEMPLATE(T, t) op2, const TEMPLATE(T, ctx_t) ctx);

#endif
