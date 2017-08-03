/*
    Copyright (C) 2017 Luca De Feo

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#ifdef T

#include "templates.h"

#ifdef __cplusplus
// extern "C" {
#endif

FLINT_DLL void TEMPLATE(T, embed_gens)(TEMPLATE(T, t) gen_sub,
				       TEMPLATE(T, t) gen_sup,
				       const TEMPLATE(T, ctx_t) sub_ctx,
				       const TEMPLATE(T, ctx_t) sup_ctx);
FLINT_DLL void _TEMPLATE(T, embed_gens_naive)(TEMPLATE(T, t) gen_sub,
					      TEMPLATE(T, t) gen_sup,
					      const TEMPLATE(T, ctx_t) sub_ctx,
					      const TEMPLATE(T, ctx_t) sup_ctx);
FLINT_DLL void _TEMPLATE(T, embed_gens_allombert)(TEMPLATE(T, t) gen_sub,
						  TEMPLATE(T, t) gen_sup,
						  const TEMPLATE(T, ctx_t) sub_ctx,
						  const TEMPLATE(T, ctx_t) sup_ctx);

#ifdef __cplusplus
// }
#endif

#endif
