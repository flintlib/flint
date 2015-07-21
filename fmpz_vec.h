/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2010 William Hart
    Copyright (C) 2010 Sebastian Pancratz

******************************************************************************/

#ifndef FMPZ_VEC_H
#define FMPZ_VEC_H

#ifdef FMPZ_VEC_INLINES_C
#define FMPZ_VEC_INLINE FLINT_DLL
#else
#define FMPZ_VEC_INLINE static __inline__
#endif

#include <gmp.h>
#include "fmpz.h"
#include "flint.h"
#include "mpf_vec.h"

#ifdef __cplusplus
 extern "C" {
#endif

#define FMPZ_VEC_NORM(vec, i)          \
do {                                   \
    while ((i) && vec[(i) - 1] == WORD(0))  \
        (i)--;                         \
} while (0)

#define FMPZ_VEC_SWAP(vec1, len1, vec2, len2) \
do {                                          \
    fmpz *__t;                                \
    slong __tn;                                \
    __t    = (vec1);                          \
    (vec1) = (vec2);                          \
    (vec2) = __t;                             \
    __tn   = (len1);                          \
    (len1) = (len2);                          \
    (len2) = __tn;                            \
} while (0);

/*  Memory management  *******************************************************/

FLINT_DLL fmpz * _fmpz_vec_init(slong len);

FLINT_DLL void _fmpz_vec_clear(fmpz * vec, slong len);

/*  Randomisation  ***********************************************************/

FLINT_DLL void _fmpz_vec_randtest(fmpz * f, flint_rand_t state, 
                        slong len, mp_bitcnt_t bits);

FLINT_DLL void _fmpz_vec_randtest_unsigned(fmpz * f, flint_rand_t state, 
                                 slong len, mp_bitcnt_t bits);

/*  Norms  *******************************************************************/

FLINT_DLL slong _fmpz_vec_max_bits(const fmpz * vec, slong len);

FLINT_DLL slong _fmpz_vec_max_bits_ref(const fmpz * vec, slong len);

FLINT_DLL mp_size_t _fmpz_vec_max_limbs(const fmpz * vec, slong len);

FLINT_DLL void _fmpz_vec_height(fmpz_t height, const fmpz * vec, slong len);

FLINT_DLL slong _fmpz_vec_height_index(const fmpz * vec, slong len);

/*  Input and output  ********************************************************/

FLINT_DLL int _fmpz_vec_fprint(FILE * file, const fmpz * vec, slong len);

FMPZ_VEC_INLINE
int _fmpz_vec_print(const fmpz * vec, slong len)
{
    return _fmpz_vec_fprint(stdout, vec, len);
}

FLINT_DLL int _fmpz_vec_fread(FILE * file, fmpz ** vec, slong * len);

FMPZ_VEC_INLINE
int _fmpz_vec_read(fmpz ** vec, slong * len)
{
    return _fmpz_vec_fread(stdin, vec, len);
}

/*  Conversions  *************************************************************/

FLINT_DLL void _fmpz_vec_set_nmod_vec(fmpz * res, 
                                       mp_srcptr poly, slong len, nmod_t mod);

FLINT_DLL void _fmpz_vec_get_nmod_vec(mp_ptr res, 
                                    const fmpz * poly, slong len, nmod_t mod);

FLINT_DLL slong _fmpz_vec_get_fft(mp_limb_t ** coeffs_f, 
                                 const fmpz * coeffs_m, slong l, slong length);

FLINT_DLL void _fmpz_vec_set_fft(fmpz * coeffs_m, slong length,
                               const mp_ptr * coeffs_f, slong limbs, slong sign);

FLINT_DLL slong _fmpz_vec_get_d_vec_2exp(double * appv, const fmpz * vec, slong len);

FLINT_DLL void _fmpz_vec_get_mpf_vec(mpf * appv, const fmpz * vec, slong len);

/*  Assignment and basic manipulation  ***************************************/

FLINT_DLL void _fmpz_vec_set(fmpz * vec1, const fmpz * vec2, slong len2);

FLINT_DLL void _fmpz_vec_swap(fmpz * vec1, fmpz * vec2, slong len2);

FLINT_DLL void _fmpz_vec_zero(fmpz * vec, slong len);

FLINT_DLL void _fmpz_vec_neg(fmpz * vec1, const fmpz * vec2, slong len2);

/*  Comparison  **************************************************************/

FLINT_DLL int _fmpz_vec_equal(const fmpz * vec1, const fmpz * vec2, slong len);

FLINT_DLL int _fmpz_vec_is_zero(const fmpz * vec, slong len);

/* Sorting  ******************************************************************/

FLINT_DLL void _fmpz_vec_sort(fmpz * vec, slong len);

/*  Addition  ****************************************************************/

FLINT_DLL void _fmpz_vec_add(fmpz * res, const fmpz * vec1, 
                                               const fmpz * vec2, slong len2);

FLINT_DLL void _fmpz_vec_sub(fmpz * res, const fmpz * vec1, 
                                               const fmpz * vec2, slong len2);

/*  Scalar multiplication and division  **************************************/

FLINT_DLL void _fmpz_vec_scalar_mul_si(fmpz * vec1, 
                                       const fmpz * vec2, slong len2, slong c);

FLINT_DLL void _fmpz_vec_scalar_mul_ui(fmpz * vec1, 
                                      const fmpz * vec2, slong len2, ulong c);

FLINT_DLL void _fmpz_vec_scalar_mul_fmpz(fmpz * vec1, 
                               const fmpz * vec2, slong len2, const fmpz_t x);

FLINT_DLL void _fmpz_vec_scalar_mul_2exp(fmpz * vec1, 
                                    const fmpz * vec2, slong len2, ulong exp);

FLINT_DLL void _fmpz_vec_scalar_divexact_fmpz(fmpz * vec1, const fmpz * vec2, 
                                                  slong len2, const fmpz_t x);

FLINT_DLL void _fmpz_vec_scalar_divexact_si(fmpz * vec1, 
                                       const fmpz * vec2, slong len2, slong c);

FLINT_DLL void _fmpz_vec_scalar_divexact_ui(fmpz * vec1, 
                                      const fmpz * vec2, slong len2, ulong c);

FLINT_DLL void _fmpz_vec_scalar_fdiv_q_fmpz(fmpz * vec1, 
                               const fmpz * vec2, slong len2, const fmpz_t c);

FLINT_DLL void _fmpz_vec_scalar_fdiv_q_si(fmpz * vec1, 
                                       const fmpz * vec2, slong len2, slong c);

FLINT_DLL void _fmpz_vec_scalar_fdiv_q_ui(fmpz * vec1, 
                                      const fmpz * vec2, slong len2, ulong c);

FLINT_DLL void _fmpz_vec_scalar_fdiv_q_2exp(fmpz * vec1, const fmpz * vec2, 
                                                       slong len2, ulong exp);

FLINT_DLL void _fmpz_vec_scalar_fdiv_r_2exp(fmpz * vec1, const fmpz * vec2, 
                                                       slong len2, ulong exp);

FLINT_DLL void _fmpz_vec_scalar_tdiv_q_fmpz(fmpz * vec1, 
                               const fmpz * vec2, slong len2, const fmpz_t c);

FLINT_DLL void _fmpz_vec_scalar_tdiv_q_si(fmpz * vec1, 
                                       const fmpz * vec2, slong len2, slong c);

FLINT_DLL void _fmpz_vec_scalar_tdiv_q_ui(fmpz * vec1, 
                                      const fmpz * vec2, slong len2, ulong c);

FLINT_DLL void _fmpz_vec_scalar_tdiv_q_2exp(fmpz * vec1, const fmpz * vec2, 
                                                       slong len2, ulong exp);

FLINT_DLL void _fmpz_vec_scalar_addmul_si(fmpz * vec1, 
                                       const fmpz * vec2, slong len2, slong c);

FLINT_DLL void _fmpz_vec_scalar_addmul_fmpz(fmpz * poly1, const fmpz * poly2, 
                                                  slong len2, const fmpz_t x);

FLINT_DLL void _fmpz_vec_scalar_addmul_si_2exp(fmpz * vec1, const fmpz * vec2, 
                                               slong len2, slong c, ulong exp);

FLINT_DLL void _fmpz_vec_scalar_submul_si(fmpz * vec1, 
                                       const fmpz * vec2, slong len2, slong c);

FLINT_DLL void _fmpz_vec_scalar_submul_fmpz(fmpz * vec1, const fmpz * vec2, 
                                                  slong len2, const fmpz_t x);

FLINT_DLL void _fmpz_vec_scalar_submul_si_2exp(fmpz * vec1, const fmpz * vec2, 
                                               slong len2, slong c, ulong exp);

/*  Vector sum and product  **************************************************/

FLINT_DLL void _fmpz_vec_sum(fmpz_t res, const fmpz * vec, slong len);

FLINT_DLL void _fmpz_vec_prod(fmpz_t res, const fmpz * vec, slong len);

/*  Reduction mod p **********************************************************/

FLINT_DLL void _fmpz_vec_scalar_mod_fmpz(fmpz *res, const fmpz *vec, slong len, const fmpz_t p);

FLINT_DLL void _fmpz_vec_scalar_smod_fmpz(fmpz *res, const fmpz *vec, slong len, const fmpz_t p);

/*  Gaussian content  ********************************************************/

FLINT_DLL void _fmpz_vec_content(fmpz_t res, const fmpz * vec, slong len);

FLINT_DLL void _fmpz_vec_lcm(fmpz_t res, const fmpz * vec, slong len);

/*  Dot product  *************************************************************/

FLINT_DLL void _fmpz_vec_dot(fmpz_t res, const fmpz * vec1, const fmpz * vec2, slong len2);

#ifdef __cplusplus
}
#endif

#endif

