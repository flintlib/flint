/*
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "templates.h"

#ifdef __cplusplus
 extern "C" {
#endif

/*  Memory management  *******************************************************/

FLINT_DLL TEMPLATE(T, struct) * _TEMPLATE(T, vec_init)(slong len, const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void _TEMPLATE(T, vec_clear)(TEMPLATE(T, struct) * vec, slong len,
                             const TEMPLATE(T, ctx_t) ctx);

/*  Randomisation  ***********************************************************/
FLINT_DLL void _TEMPLATE(T, vec_randtest)(TEMPLATE(T, struct) * f, 
                           flint_rand_t state, 
                           slong len, 
                           const TEMPLATE(T, ctx_t) ctx);

/*  Norms  *******************************************************************/

/*  Input and output  ********************************************************/

#if defined (H_STDIO)               \
  || defined (_H_STDIO)             \
  || defined (_STDIO_H)             \
  || defined (_STDIO_H_)            \
  || defined (__STDIO_H)            \
  || defined (__STDIO_H__)          \
  || defined (_STDIO_INCLUDED)      \
  || defined (__dj_include_stdio_h_)\
  || defined (_FILE_DEFINED)        \
  || defined (__STDIO__)            \
  || defined (_MSL_STDIO_H)         \
  || defined (_STDIO_H_INCLUDED)    \
  || defined (_ISO_STDIO_ISO_H)     \
  || defined (__STDIO_LOADED)       \
  || defined (_STDIO)

FLINT_DLL int _TEMPLATE(T, vec_fprint)(FILE * file,
                             const TEMPLATE(T, struct) * vec,
                             slong len,
                             const TEMPLATE(T, ctx_t) ctx);

FQ_VEC_TEMPLATES_INLINE
int _TEMPLATE(T, vec_print)(const TEMPLATE(T, struct) * vec, slong len,
                            const TEMPLATE(T, ctx_t) ctx)
{
    return _TEMPLATE(T, vec_fprint)(stdout, vec, len, ctx);
}

#endif

/*  Conversions  *************************************************************/

/*  Assignment and basic manipulation  ***************************************/
FLINT_DLL void _TEMPLATE(T, vec_set)(TEMPLATE(T, struct) * v,
                      const TEMPLATE(T, struct) * f, 
                      slong len, const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void _TEMPLATE(T, vec_swap)(TEMPLATE(T, struct) * vec1,
                       TEMPLATE(T, struct) * vec2,
                       slong len2, const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void _TEMPLATE(T, vec_zero)(TEMPLATE(T, struct) * v, 
                       slong len,
                       const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void _TEMPLATE(T, vec_neg)(TEMPLATE(T, struct) * vec1,
                      const TEMPLATE(T, struct) * vec2,
                      slong len2,
                      const TEMPLATE(T, ctx_t) ctx);

/*  Comparison  **************************************************************/
FLINT_DLL int _TEMPLATE(T, vec_is_zero)(const TEMPLATE(T, struct) * vec, slong len,
                          const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL int _TEMPLATE(T, vec_equal)(const TEMPLATE(T, struct) * vec1,
                        const TEMPLATE(T, struct) * vec2, slong len,
                        const TEMPLATE(T, ctx_t) ctx);

/*  Sorting  *****************************************************************/

/*  Addition  ****************************************************************/

FLINT_DLL void _TEMPLATE(T, vec_add)(TEMPLATE(T, struct) * res,
                           const TEMPLATE(T, struct) * vec1, 
                           const TEMPLATE(T, struct) * vec2, 
                           slong len2,
                           const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void _TEMPLATE(T, vec_sub)(TEMPLATE(T, struct) * res,
                           const TEMPLATE(T, struct) * vec1, 
                           const TEMPLATE(T, struct) * vec2, 
                           slong len2,
                           const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void _TEMPLATE(T, TEMPLATE(vec_scalar_addmul, T))(TEMPLATE(T, struct) * poly1, 
                                             const TEMPLATE(T, struct) * poly2,
                                             slong len2,
                                             const TEMPLATE(T, t) x,
                                             const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void _TEMPLATE(T, TEMPLATE(vec_scalar_submul, T))(TEMPLATE(T, struct) * poly1, 
                                             const TEMPLATE(T, struct) * poly2,
                                             slong len2,
                                             const TEMPLATE(T, t) x,
                                             const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void _TEMPLATE3(T, vec_scalar_mul, T) (TEMPLATE(T, struct) * poly1,
                                  const TEMPLATE(T, struct) * poly2,
                                  slong len2, const TEMPLATE(T, t) x,
                                  const TEMPLATE(T, ctx_t) ctx);
     
/* ****************************************************************************/
FLINT_DLL void _TEMPLATE(T, vec_dot)(TEMPLATE(T, t) res,
                      const TEMPLATE(T, struct) * vec1,
                      const TEMPLATE(T, struct) * vec2,
                      slong len2,
                      const TEMPLATE(T, ctx_t) ctx);

#ifdef __cplusplus
 }
#endif

#endif
