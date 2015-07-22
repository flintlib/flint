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
    Copyright (C) 2010 Fredrik Johansson

******************************************************************************/

#ifndef FMPZ_MAT_H
#define FMPZ_MAT_H

#ifdef FMPZ_MAT_INLINES_C
#define FMPZ_MAT_INLINE FLINT_DLL
#else
#define FMPZ_MAT_INLINE static __inline__
#endif

#undef ulong
#define ulong ulongxx /* interferes with system includes */
#include <stdio.h>
#undef ulong

#include <gmp.h>
#define ulong mp_limb_t
#include "flint.h"
#include "fmpz.h"
#include "nmod_mat.h"
#include "d_mat.h"
#include "mpf_mat.h"
#include "fmpz_poly.h"

#ifdef __cplusplus
 extern "C" {
#endif

typedef struct
{
    fmpz * entries;
    slong r;
    slong c;
    fmpz ** rows;
} fmpz_mat_struct;

typedef fmpz_mat_struct fmpz_mat_t[1];

/* Memory management  ********************************************************/

FMPZ_MAT_INLINE
fmpz * fmpz_mat_entry(const fmpz_mat_t mat, slong i, slong j)
{
   return mat->rows[i] + j;
}

FMPZ_MAT_INLINE
slong fmpz_mat_nrows(const fmpz_mat_t mat)
{
   return mat->r;
}

FMPZ_MAT_INLINE
slong fmpz_mat_ncols(const fmpz_mat_t mat)
{
   return mat->c;
}

FLINT_DLL void fmpz_mat_init(fmpz_mat_t mat, slong rows, slong cols);
FLINT_DLL void fmpz_mat_init_set(fmpz_mat_t mat, const fmpz_mat_t src);
FLINT_DLL void fmpz_mat_swap(fmpz_mat_t mat1, fmpz_mat_t mat2);
FLINT_DLL void fmpz_mat_set(fmpz_mat_t mat1, const fmpz_mat_t mat2);
FLINT_DLL void fmpz_mat_clear(fmpz_mat_t mat);

FLINT_DLL int fmpz_mat_equal(const fmpz_mat_t mat1, const fmpz_mat_t mat2);
FLINT_DLL int fmpz_mat_is_zero(const fmpz_mat_t mat);
FLINT_DLL int fmpz_mat_is_one(const fmpz_mat_t mat);

FMPZ_MAT_INLINE
int fmpz_mat_is_empty(const fmpz_mat_t mat)
{
    return (mat->r == 0) || (mat->c == 0);
}

FMPZ_MAT_INLINE
int fmpz_mat_is_square(const fmpz_mat_t mat)
{
    return (mat->r == mat->c);
}

FLINT_DLL void fmpz_mat_zero(fmpz_mat_t mat);
FLINT_DLL void fmpz_mat_one(fmpz_mat_t mat);

/* Windows and concatenation */

FLINT_DLL void fmpz_mat_window_init(fmpz_mat_t window, const fmpz_mat_t mat, slong r1,
    slong c1, slong r2, slong c2);

FLINT_DLL void fmpz_mat_window_clear(fmpz_mat_t window);

FLINT_DLL void fmpz_mat_concat_horizontal(fmpz_mat_t res,
                           const fmpz_mat_t mat1,  const fmpz_mat_t mat2);

FLINT_DLL void fmpz_mat_concat_vertical(fmpz_mat_t res,
                           const fmpz_mat_t mat1,  const fmpz_mat_t mat2);

/* Input and output  *********************************************************/

FLINT_DLL int fmpz_mat_fprint(FILE * file, const fmpz_mat_t mat);

FLINT_DLL int fmpz_mat_fprint_pretty(FILE * file, const fmpz_mat_t mat);

FMPZ_MAT_INLINE
int fmpz_mat_print(const fmpz_mat_t mat)
{
    return fmpz_mat_fprint(stdout, mat);
}

FMPZ_MAT_INLINE
int fmpz_mat_print_pretty(const fmpz_mat_t mat)
{
    return fmpz_mat_fprint_pretty(stdout, mat);
}

FLINT_DLL int fmpz_mat_fread(FILE* file, fmpz_mat_t mat);

FMPZ_MAT_INLINE
int fmpz_mat_read(fmpz_mat_t mat)
{
    return fmpz_mat_fread(stdin, mat);
}

/* Random matrix generation  *************************************************/

FLINT_DLL void fmpz_mat_randbits(fmpz_mat_t mat, flint_rand_t state, mp_bitcnt_t bits);
FLINT_DLL void fmpz_mat_randtest(fmpz_mat_t mat, flint_rand_t state, mp_bitcnt_t bits);
FLINT_DLL void fmpz_mat_randtest_unsigned(fmpz_mat_t mat, flint_rand_t state, mp_bitcnt_t bits);
FLINT_DLL void fmpz_mat_randintrel(fmpz_mat_t mat, flint_rand_t state, mp_bitcnt_t bits);
FLINT_DLL void fmpz_mat_randsimdioph(fmpz_mat_t mat, flint_rand_t state, mp_bitcnt_t bits, mp_bitcnt_t bits2);
FLINT_DLL void fmpz_mat_randntrulike(fmpz_mat_t mat, flint_rand_t state, mp_bitcnt_t bits, ulong q);
FLINT_DLL void fmpz_mat_randntrulike2(fmpz_mat_t mat, flint_rand_t state, mp_bitcnt_t bits, ulong q);
FLINT_DLL void fmpz_mat_randajtai(fmpz_mat_t mat, flint_rand_t state, double alpha);
FLINT_DLL void fmpz_mat_randrank(fmpz_mat_t mat, flint_rand_t state, slong rank, mp_bitcnt_t bits);
FLINT_DLL void fmpz_mat_randdet(fmpz_mat_t mat, flint_rand_t state, const fmpz_t det);
FLINT_DLL void fmpz_mat_randops(fmpz_mat_t mat, flint_rand_t state, slong count);
FLINT_DLL int fmpz_mat_randpermdiag(fmpz_mat_t mat, flint_rand_t state, const fmpz * diag, slong n);

/* Norms */

FLINT_DLL slong fmpz_mat_max_bits(const fmpz_mat_t mat);

/* Transpose */

FLINT_DLL void fmpz_mat_transpose(fmpz_mat_t B, const fmpz_mat_t A);

/* Addition and subtraction */

FLINT_DLL void fmpz_mat_add(fmpz_mat_t C, const fmpz_mat_t A, const fmpz_mat_t B);
FLINT_DLL void fmpz_mat_sub(fmpz_mat_t C, const fmpz_mat_t A, const fmpz_mat_t B);
FLINT_DLL void fmpz_mat_neg(fmpz_mat_t B, const fmpz_mat_t A);

/* Scalar operations */
FLINT_DLL void fmpz_mat_scalar_mul_fmpz(fmpz_mat_t B, const fmpz_mat_t A, const fmpz_t c);
FLINT_DLL void fmpz_mat_scalar_mul_si(fmpz_mat_t B, const fmpz_mat_t A, slong c);
FLINT_DLL void fmpz_mat_scalar_mul_ui(fmpz_mat_t B, const fmpz_mat_t A, ulong c);

FLINT_DLL void fmpz_mat_scalar_addmul_fmpz(fmpz_mat_t B, const fmpz_mat_t A, const fmpz_t c);
FLINT_DLL void fmpz_mat_scalar_addmul_si(fmpz_mat_t B, const fmpz_mat_t A, slong c);
FLINT_DLL void fmpz_mat_scalar_addmul_ui(fmpz_mat_t B, const fmpz_mat_t A, ulong c);

FLINT_DLL void fmpz_mat_scalar_submul_fmpz(fmpz_mat_t B, const fmpz_mat_t A, const fmpz_t c);
FLINT_DLL void fmpz_mat_scalar_submul_si(fmpz_mat_t B, const fmpz_mat_t A, slong c);
FLINT_DLL void fmpz_mat_scalar_submul_ui(fmpz_mat_t B, const fmpz_mat_t A, ulong c);

FLINT_DLL void fmpz_mat_scalar_addmul_nmod_mat_fmpz(fmpz_mat_t B, const nmod_mat_t A, const fmpz_t c);
FLINT_DLL void fmpz_mat_scalar_addmul_nmod_mat_ui(fmpz_mat_t B, const nmod_mat_t A, ulong c);

FLINT_DLL void fmpz_mat_scalar_divexact_fmpz(fmpz_mat_t B, const fmpz_mat_t A, const fmpz_t c);
FLINT_DLL void fmpz_mat_scalar_divexact_si(fmpz_mat_t B, const fmpz_mat_t A, slong c);
FLINT_DLL void fmpz_mat_scalar_divexact_ui(fmpz_mat_t B, const fmpz_mat_t A, ulong c);

FLINT_DLL void fmpz_mat_scalar_mul_2exp(fmpz_mat_t B, const fmpz_mat_t A, ulong exp);
FLINT_DLL void fmpz_mat_scalar_tdiv_q_2exp(fmpz_mat_t B, const fmpz_mat_t A, ulong exp);

FLINT_DLL void fmpz_mat_scalar_mod_fmpz(fmpz_mat_t B, const fmpz_mat_t A, const fmpz_t m);

/* Multiplication */

FLINT_DLL void fmpz_mat_mul(fmpz_mat_t C, const fmpz_mat_t A, const fmpz_mat_t B);

FLINT_DLL void fmpz_mat_mul_classical(fmpz_mat_t C, const fmpz_mat_t A,
    const fmpz_mat_t B);

FLINT_DLL void fmpz_mat_mul_classical_inline(fmpz_mat_t C, const fmpz_mat_t A,
    const fmpz_mat_t B);

FLINT_DLL void _fmpz_mat_mul_multi_mod(fmpz_mat_t C, const fmpz_mat_t A,
    const fmpz_mat_t B, mp_bitcnt_t bits);

FLINT_DLL void fmpz_mat_mul_multi_mod(fmpz_mat_t C, const fmpz_mat_t A,
    const fmpz_mat_t B);

FLINT_DLL void fmpz_mat_sqr_bodrato(fmpz_mat_t B, const fmpz_mat_t A);

FLINT_DLL void fmpz_mat_sqr(fmpz_mat_t B, const fmpz_mat_t A);

FLINT_DLL void fmpz_mat_pow(fmpz_mat_t B, const fmpz_mat_t A, ulong exp);

/* Content */

FLINT_DLL void fmpz_mat_content(fmpz_t ret, const fmpz_mat_t A);

/* Permutations */

FMPZ_MAT_INLINE
void fmpz_mat_swap_rows(fmpz_mat_t mat, slong * perm, slong r, slong s)
{
    if (r != s)
    {
        fmpz * u;
        slong t;

        if (perm)
        {
            t = perm[s];
            perm[s] = perm[r];
            perm[r] = t;
        }

        u = mat->rows[s];
        mat->rows[s] = mat->rows[r];
        mat->rows[r] = u;
    }
}

/* Gaussian elimination *****************************************************/

FLINT_DLL slong fmpz_mat_find_pivot_any(const fmpz_mat_t mat,
                                    slong start_row, slong end_row, slong c);

FLINT_DLL slong fmpz_mat_fflu(fmpz_mat_t B, fmpz_t den, slong * perm,
                            const fmpz_mat_t A, int rank_check);

FLINT_DLL slong fmpz_mat_rref(fmpz_mat_t B, fmpz_t den, const fmpz_mat_t A);
FLINT_DLL slong fmpz_mat_rref_fflu(fmpz_mat_t B, fmpz_t den, const fmpz_mat_t A);
FLINT_DLL slong fmpz_mat_rref_mul(fmpz_mat_t B, fmpz_t den, const fmpz_mat_t A);
FLINT_DLL int fmpz_mat_is_in_rref_with_rank(const fmpz_mat_t A, const fmpz_t den,
        slong rank);

/* Modular gaussian elimination *********************************************/

FLINT_DLL slong fmpz_mat_rref_mod(slong * perm, fmpz_mat_t A, const fmpz_t p);

/* Trace ********************************************************************/

FLINT_DLL void fmpz_mat_trace(fmpz_t trace, const fmpz_mat_t mat);

/* Determinant **************************************************************/

FLINT_DLL void fmpz_mat_det(fmpz_t det, const fmpz_mat_t A);

FLINT_DLL void fmpz_mat_det_cofactor(fmpz_t det, const fmpz_mat_t A);
FLINT_DLL void _fmpz_mat_det_cofactor_2x2(fmpz_t det, fmpz ** const x);
FLINT_DLL void _fmpz_mat_det_cofactor_3x3(fmpz_t det, fmpz ** const x);
FLINT_DLL void _fmpz_mat_det_cofactor_4x4(fmpz_t det, fmpz ** const x);

FLINT_DLL void fmpz_mat_det_bareiss(fmpz_t det, const fmpz_mat_t A);

FLINT_DLL void fmpz_mat_det_modular(fmpz_t det, const fmpz_mat_t A, int proved);

FLINT_DLL void fmpz_mat_det_modular_accelerated(fmpz_t det,
    const fmpz_mat_t A, int proved);

FLINT_DLL void fmpz_mat_det_modular_given_divisor(fmpz_t det, const fmpz_mat_t A,
        const fmpz_t d, int proved);

FLINT_DLL void fmpz_mat_det_bound(fmpz_t bound, const fmpz_mat_t A);
FLINT_DLL void fmpz_mat_det_divisor(fmpz_t d, const fmpz_mat_t A);

/* Characteristic polynomial ************************************************/

FLINT_DLL void _fmpz_mat_charpoly(fmpz *cp, const fmpz_mat_t mat);
FLINT_DLL void fmpz_mat_charpoly(fmpz_poly_t cp, const fmpz_mat_t mat);

/* Rank *********************************************************************/

FLINT_DLL slong fmpz_mat_rank(const fmpz_mat_t A);

/* Nonsingular solving ******************************************************/

FLINT_DLL void fmpz_mat_solve_bound(fmpz_t N, fmpz_t D,
        const fmpz_mat_t A, const fmpz_mat_t B);

FLINT_DLL int fmpz_mat_solve(fmpz_mat_t X, fmpz_t den,
        const fmpz_mat_t A, const fmpz_mat_t B);

FLINT_DLL int fmpz_mat_solve_cramer(fmpz_mat_t X, fmpz_t den,
        const fmpz_mat_t A, const fmpz_mat_t B);

FLINT_DLL int fmpz_mat_solve_fflu(fmpz_mat_t X, fmpz_t den,
        const fmpz_mat_t A, const fmpz_mat_t B);

FLINT_DLL void fmpz_mat_solve_fflu_precomp(fmpz_mat_t X, const slong * perm,
        const fmpz_mat_t FFLU, const fmpz_mat_t B);

FLINT_DLL int fmpz_mat_solve_dixon(fmpz_mat_t X, fmpz_t mod,
        const fmpz_mat_t A, const fmpz_mat_t B);

/* Nullspace ****************************************************************/

FLINT_DLL slong fmpz_mat_nullspace(fmpz_mat_t res, const fmpz_mat_t mat);

/* Inverse ******************************************************************/

FLINT_DLL int fmpz_mat_inv(fmpz_mat_t B, fmpz_t den, const fmpz_mat_t A);

/* Modular reduction and reconstruction *************************************/

FLINT_DLL void fmpz_mat_set_nmod_mat(fmpz_mat_t A, const nmod_mat_t Amod);

FLINT_DLL void fmpz_mat_set_nmod_mat_unsigned(fmpz_mat_t A, const nmod_mat_t Amod);

FLINT_DLL void fmpz_mat_get_nmod_mat(nmod_mat_t Amod, const fmpz_mat_t A);

FLINT_DLL void fmpz_mat_CRT_ui(fmpz_mat_t res, const fmpz_mat_t mat1,
                        const fmpz_t m1, const nmod_mat_t mat2, int sign);

FLINT_DLL void fmpz_mat_multi_mod_ui_precomp(nmod_mat_t * residues, slong nres, 
    const fmpz_mat_t mat, const fmpz_comb_t comb, fmpz_comb_temp_t temp);

FLINT_DLL void fmpz_mat_multi_mod_ui(nmod_mat_t * residues, slong nres, const fmpz_mat_t mat);


FLINT_DLL void fmpz_mat_multi_CRT_ui_precomp(fmpz_mat_t mat,
    nmod_mat_t * const residues, slong nres,
    const fmpz_comb_t comb, fmpz_comb_temp_t temp, int sign);

FLINT_DLL void fmpz_mat_multi_CRT_ui(fmpz_mat_t mat, nmod_mat_t * const residues,
    slong nres, int sign);

/* HNF and SNF **************************************************************/

FLINT_DLL void fmpz_mat_hnf(fmpz_mat_t H, const fmpz_mat_t A);
FLINT_DLL void fmpz_mat_hnf_transform(fmpz_mat_t H, fmpz_mat_t U, const fmpz_mat_t A);
FLINT_DLL void fmpz_mat_hnf_classical(fmpz_mat_t H, const fmpz_mat_t A);
FLINT_DLL void fmpz_mat_hnf_xgcd(fmpz_mat_t H, const fmpz_mat_t A);
FLINT_DLL void fmpz_mat_hnf_minors(fmpz_mat_t H, const fmpz_mat_t A);
FLINT_DLL void fmpz_mat_hnf_modular(fmpz_mat_t H, const fmpz_mat_t A, const fmpz_t D);
FLINT_DLL int fmpz_mat_hnf_pernet_stein(fmpz_mat_t H, const fmpz_mat_t A, flint_rand_t state);
FLINT_DLL int fmpz_mat_is_in_hnf(const fmpz_mat_t A);

FLINT_DLL void fmpz_mat_snf(fmpz_mat_t S, const fmpz_mat_t A);
FLINT_DLL void fmpz_mat_snf_diagonal(fmpz_mat_t S, const fmpz_mat_t A);
FLINT_DLL void fmpz_mat_snf_kannan_bachem(fmpz_mat_t S, const fmpz_mat_t A);
FLINT_DLL void fmpz_mat_snf_iliopoulos(fmpz_mat_t S, const fmpz_mat_t A,
        const fmpz_t mod);
FLINT_DLL int fmpz_mat_is_in_snf(const fmpz_mat_t A);

/* Special matrices **********************************************************/

FLINT_DLL int fmpz_mat_is_hadamard(const fmpz_mat_t A);

FLINT_DLL int fmpz_mat_hadamard(fmpz_mat_t A);

/* Gram matrix **************************************************************/

FLINT_DLL void fmpz_mat_gram(fmpz_mat_t B, const fmpz_mat_t A);

/* Conversions **************************************************************/

FLINT_DLL int fmpz_mat_get_d_mat(d_mat_t B, const fmpz_mat_t A);

FLINT_DLL int fmpz_mat_get_d_mat_transpose(d_mat_t B, const fmpz_mat_t A);

FLINT_DLL void fmpz_mat_get_mpf_mat(mpf_mat_t B, const fmpz_mat_t A);

/* Cholesky Decomposition ****************************************************/

FLINT_DLL void fmpz_mat_chol_d(d_mat_t R, const fmpz_mat_t A);

/* LLL ***********************************************************************/

FLINT_DLL int fmpz_mat_is_reduced(const fmpz_mat_t A, double delta, double eta);

FLINT_DLL int fmpz_mat_is_reduced_gram(const fmpz_mat_t A, double delta, double eta);

FLINT_DLL int fmpz_mat_is_reduced_with_removal(const fmpz_mat_t A, double delta, double eta, const fmpz_t gs_B, int newd);

FLINT_DLL int fmpz_mat_is_reduced_gram_with_removal(const fmpz_mat_t A, double delta, double eta, const fmpz_t gs_B, int newd);

/* Classical LLL *************************************************************/

FLINT_DLL void fmpz_mat_lll_original(fmpz_mat_t A, const fmpq_t delta, const fmpq_t eta);

/* Modified LLL **************************************************************/

FLINT_DLL void fmpz_mat_lll_storjohann(fmpz_mat_t A, const fmpq_t delta, const fmpq_t eta);

#ifdef __cplusplus
}
#endif

#endif

