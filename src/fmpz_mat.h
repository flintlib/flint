/*
    Copyright (C) 2010 William Hart
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef FMPZ_MAT_H
#define FMPZ_MAT_H

#ifdef FMPZ_MAT_INLINES_C
#define FMPZ_MAT_INLINE
#else
#define FMPZ_MAT_INLINE static inline
#endif

#include "fmpz_types.h"
#include "nmod_types.h"
#include "d_mat.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Element access  ********************************************************/

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

/* Memory management  ********************************************************/

void fmpz_mat_init(fmpz_mat_t mat, slong rows, slong cols);
void fmpz_mat_init_set(fmpz_mat_t mat, const fmpz_mat_t src);

void fmpz_mat_swap(fmpz_mat_t mat1, fmpz_mat_t mat2);
void fmpz_mat_swap_entrywise(fmpz_mat_t mat1, fmpz_mat_t mat2);
void fmpz_mat_swap_cols(fmpz_mat_t mat, slong * perm, slong r, slong s);

void fmpz_mat_set(fmpz_mat_t mat1, const fmpz_mat_t mat2);
void fmpz_mat_clear(fmpz_mat_t mat);

int fmpz_mat_equal(const fmpz_mat_t mat1, const fmpz_mat_t mat2);
int fmpz_mat_is_zero(const fmpz_mat_t mat);
int fmpz_mat_is_one(const fmpz_mat_t mat);

int fmpz_mat_is_zero_row(const fmpz_mat_t mat, slong i);

#define fmpz_mat_col_equal(M, m, n) _Pragma("GCC warning \"'fmpz_mat_col_equal' is deprecated in favor for 'fmpz_mat_equal_col'\"") fmpz_mat_equal_col(M, m, n)
#define fmpz_mat_row_equal(M, m, n) _Pragma("GCC warning \"'fmpz_mat_row_equal' is deprecated in favor for 'fmpz_mat_equal_row'\"") fmpz_mat_equal_row(M, m, n)
int fmpz_mat_equal_col(fmpz_mat_t M, slong m, slong n);
int fmpz_mat_equal_row(fmpz_mat_t M, slong m, slong n);

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

void fmpz_mat_zero(fmpz_mat_t mat);
void fmpz_mat_one(fmpz_mat_t mat);

/* Windows and concatenation *************************************************/

void fmpz_mat_window_init(fmpz_mat_t window, const fmpz_mat_t mat, slong r1, slong c1, slong r2, slong c2);
void fmpz_mat_window_clear(fmpz_mat_t window);

void _fmpz_mat_window_readonly_init_strip_initial_zero_rows(fmpz_mat_t A, const fmpz_mat_t B);
#define _fmpz_mat_window_readonly_clear(A) /* Do nothing */

void fmpz_mat_concat_horizontal(fmpz_mat_t res, const fmpz_mat_t mat1, const fmpz_mat_t mat2);
void fmpz_mat_concat_vertical(fmpz_mat_t res, const fmpz_mat_t mat1, const fmpz_mat_t mat2);

/* Input and output  *********************************************************/

#ifdef FLINT_HAVE_FILE
int fmpz_mat_fprint(FILE * file, const fmpz_mat_t mat);
int fmpz_mat_fprint_pretty(FILE * file, const fmpz_mat_t mat);

int fmpz_mat_fread(FILE* file, fmpz_mat_t mat);
#endif

int fmpz_mat_print(const fmpz_mat_t mat);
int fmpz_mat_print_pretty(const fmpz_mat_t mat);

int fmpz_mat_read(fmpz_mat_t mat);

/* Random matrix generation  *************************************************/

void fmpz_mat_randbits(fmpz_mat_t mat, flint_rand_t state, flint_bitcnt_t bits);
void fmpz_mat_randtest(fmpz_mat_t mat, flint_rand_t state, flint_bitcnt_t bits);
void fmpz_mat_randtest_unsigned(fmpz_mat_t mat, flint_rand_t state, flint_bitcnt_t bits);
void fmpz_mat_randintrel(fmpz_mat_t mat, flint_rand_t state, flint_bitcnt_t bits);
void fmpz_mat_randsimdioph(fmpz_mat_t mat, flint_rand_t state, flint_bitcnt_t bits, flint_bitcnt_t bits2);
void fmpz_mat_randntrulike(fmpz_mat_t mat, flint_rand_t state, flint_bitcnt_t bits, ulong q);
void fmpz_mat_randntrulike2(fmpz_mat_t mat, flint_rand_t state, flint_bitcnt_t bits, ulong q);
void fmpz_mat_randajtai(fmpz_mat_t mat, flint_rand_t state, double alpha);
void fmpz_mat_randrank(fmpz_mat_t mat, flint_rand_t state, slong rank, flint_bitcnt_t bits);
void fmpz_mat_randdet(fmpz_mat_t mat, flint_rand_t state, const fmpz_t det);
void fmpz_mat_randops(fmpz_mat_t mat, flint_rand_t state, slong count);
int fmpz_mat_randpermdiag(fmpz_mat_t mat, flint_rand_t state, const fmpz * diag, slong n);

/* Norms */

slong fmpz_mat_max_bits(const fmpz_mat_t mat);

/* Transpose */

void fmpz_mat_transpose(fmpz_mat_t B, const fmpz_mat_t A);

/* Addition and subtraction */

void fmpz_mat_add(fmpz_mat_t C, const fmpz_mat_t A, const fmpz_mat_t B);
void fmpz_mat_sub(fmpz_mat_t C, const fmpz_mat_t A, const fmpz_mat_t B);
void fmpz_mat_neg(fmpz_mat_t B, const fmpz_mat_t A);

/* Scalar operations */
void fmpz_mat_scalar_mul_fmpz(fmpz_mat_t B, const fmpz_mat_t A, const fmpz_t c);
void fmpz_mat_scalar_mul_si(fmpz_mat_t B, const fmpz_mat_t A, slong c);
void fmpz_mat_scalar_mul_ui(fmpz_mat_t B, const fmpz_mat_t A, ulong c);

void fmpz_mat_scalar_addmul_fmpz(fmpz_mat_t B, const fmpz_mat_t A, const fmpz_t c);
void fmpz_mat_scalar_addmul_si(fmpz_mat_t B, const fmpz_mat_t A, slong c);
void fmpz_mat_scalar_addmul_ui(fmpz_mat_t B, const fmpz_mat_t A, ulong c);

void fmpz_mat_scalar_submul_fmpz(fmpz_mat_t B, const fmpz_mat_t A, const fmpz_t c);
void fmpz_mat_scalar_submul_si(fmpz_mat_t B, const fmpz_mat_t A, slong c);
void fmpz_mat_scalar_submul_ui(fmpz_mat_t B, const fmpz_mat_t A, ulong c);

void fmpz_mat_scalar_addmul_nmod_mat_fmpz(fmpz_mat_t B, const nmod_mat_t A, const fmpz_t c);
void fmpz_mat_scalar_addmul_nmod_mat_ui(fmpz_mat_t B, const nmod_mat_t A, ulong c);

void fmpz_mat_scalar_divexact_fmpz(fmpz_mat_t B, const fmpz_mat_t A, const fmpz_t c);
void fmpz_mat_scalar_divexact_si(fmpz_mat_t B, const fmpz_mat_t A, slong c);
void fmpz_mat_scalar_divexact_ui(fmpz_mat_t B, const fmpz_mat_t A, ulong c);

void fmpz_mat_scalar_mul_2exp(fmpz_mat_t B, const fmpz_mat_t A, ulong exp);
void fmpz_mat_scalar_tdiv_q_2exp(fmpz_mat_t B, const fmpz_mat_t A, ulong exp);
void fmpz_mat_scalar_smod(fmpz_mat_t B, const fmpz_mat_t A, const fmpz_t P);

void fmpz_mat_scalar_mod_fmpz(fmpz_mat_t B, const fmpz_mat_t A, const fmpz_t m);

/* Multiplication */

void fmpz_mat_mul(fmpz_mat_t C, const fmpz_mat_t A, const fmpz_mat_t B);

void fmpz_mat_mul_classical(fmpz_mat_t C, const fmpz_mat_t A,
    const fmpz_mat_t B);

void fmpz_mat_mul_strassen(fmpz_mat_t C, const fmpz_mat_t A, const fmpz_mat_t B);

void fmpz_mat_mul_classical_inline(fmpz_mat_t C, const fmpz_mat_t A,
    const fmpz_mat_t B);

void _fmpz_mat_mul_fft(fmpz_mat_t C,
                                    const fmpz_mat_t A, slong abits,
                                    const fmpz_mat_t B, slong bbits, int sign);

void fmpz_mat_mul_fft(fmpz_mat_t C, const fmpz_mat_t A, const fmpz_mat_t B);

void _fmpz_mat_mul_multi_mod(fmpz_mat_t C, const fmpz_mat_t A,
                           const fmpz_mat_t B, int sign, flint_bitcnt_t Cbits);

void fmpz_mat_mul_multi_mod(fmpz_mat_t C, const fmpz_mat_t A,
                                                           const fmpz_mat_t B);

int _fmpz_mat_mul_blas(fmpz_mat_t C,
                                    const fmpz_mat_t A, flint_bitcnt_t Abits,
                                    const fmpz_mat_t B, flint_bitcnt_t Bbits,
                                               int sign, flint_bitcnt_t Cbits);

int fmpz_mat_mul_blas(fmpz_mat_t C, const fmpz_mat_t A,
                                                           const fmpz_mat_t B);

void _fmpz_mat_mul_small_1(fmpz_mat_t C, const fmpz_mat_t A,
                                                           const fmpz_mat_t B);

void _fmpz_mat_mul_small_2a(fmpz_mat_t C, const fmpz_mat_t A,
                                                           const fmpz_mat_t B);

void _fmpz_mat_mul_small_2b(fmpz_mat_t C, const fmpz_mat_t A,
                                                           const fmpz_mat_t B);

void _fmpz_mat_mul_small_internal(fmpz_mat_t C, const fmpz_mat_t A,
                                     const fmpz_mat_t B, flint_bitcnt_t Cbits);

void _fmpz_mat_mul_small(fmpz_mat_t C, const fmpz_mat_t A,
                                                           const fmpz_mat_t B);

void _fmpz_mat_mul_double_word(fmpz_mat_t C, const fmpz_mat_t A,
                                                           const fmpz_mat_t B);

void _fmpz_mat_mul_double_word_internal(fmpz_mat_t C,
        const fmpz_mat_t A, const fmpz_mat_t B, int sign, flint_bitcnt_t bits);

void fmpz_mat_sqr_bodrato(fmpz_mat_t B, const fmpz_mat_t A);

void fmpz_mat_sqr(fmpz_mat_t B, const fmpz_mat_t A);

void fmpz_mat_pow(fmpz_mat_t B, const fmpz_mat_t A, ulong exp);

void fmpz_mat_mul_fmpz_vec(fmpz * c, const fmpz_mat_t A,
                                                   const fmpz * b, slong blen);

void fmpz_mat_mul_fmpz_vec_ptr(fmpz * const * c, const fmpz_mat_t A,
                                           const fmpz * const * b, slong blen);

void fmpz_mat_fmpz_vec_mul(fmpz * c, const fmpz * a, slong alen,
                                                           const fmpz_mat_t B);

void fmpz_mat_fmpz_vec_mul_ptr(fmpz * const * c,
                       const fmpz * const * a, slong alen, const fmpz_mat_t B);

/* Kronecker product */

void fmpz_mat_kronecker_product(fmpz_mat_t C, const fmpz_mat_t A, const fmpz_mat_t B);

/* Content */

void fmpz_mat_content(fmpz_t ret, const fmpz_mat_t A);

/* Permutations */

FMPZ_MAT_INLINE
void fmpz_mat_swap_rows(fmpz_mat_t mat, slong * perm, slong r, slong s)
{
    if (r != s && !fmpz_mat_is_empty(mat))
    {
        if (perm != NULL)
            FLINT_SWAP(slong, perm[r], perm[s]);

        FLINT_SWAP(fmpz *, mat->rows[r], mat->rows[s]);
    }
}

void fmpz_mat_invert_rows(fmpz_mat_t mat, slong * perm);
void fmpz_mat_invert_cols(fmpz_mat_t mat, slong * perm);

/* Gaussian elimination *****************************************************/

slong fmpz_mat_find_pivot_any(const fmpz_mat_t mat,
                                      slong start_row, slong end_row, slong c);

slong fmpz_mat_find_pivot_smallest(const fmpz_mat_t mat,
                                      slong start_row, slong end_row, slong c);

slong fmpz_mat_fflu(fmpz_mat_t B, fmpz_t den, slong * perm,
                            const fmpz_mat_t A, int rank_check);

slong fmpz_mat_rank_small_inplace(fmpz_mat_t B);
slong fmpz_mat_rref(fmpz_mat_t B, fmpz_t den, const fmpz_mat_t A);
slong fmpz_mat_rref_fflu(fmpz_mat_t B, fmpz_t den, const fmpz_mat_t A);
slong fmpz_mat_rref_mul(fmpz_mat_t B, fmpz_t den, const fmpz_mat_t A);
int fmpz_mat_is_in_rref_with_rank(const fmpz_mat_t A, const fmpz_t den,
        slong rank);

/* Modular gaussian elimination *********************************************/

slong fmpz_mat_rref_mod(slong * perm, fmpz_mat_t A, const fmpz_t p);

/* Modular Howell and strong echelon form ***********************************/

slong fmpz_mat_howell_form_mod(fmpz_mat_t A, const fmpz_t mod);

void fmpz_mat_strong_echelon_form_mod(fmpz_mat_t A, const fmpz_t mod);

/* Trace ********************************************************************/

void fmpz_mat_trace(fmpz_t trace, const fmpz_mat_t mat);

/* Determinant **************************************************************/

void fmpz_mat_det(fmpz_t det, const fmpz_mat_t A);

void fmpz_mat_det_cofactor(fmpz_t det, const fmpz_mat_t A);

void fmpz_mat_det_bareiss(fmpz_t det, const fmpz_mat_t A);

void fmpz_mat_det_modular(fmpz_t det, const fmpz_mat_t A, int proved);

void fmpz_mat_det_modular_accelerated(fmpz_t det,
    const fmpz_mat_t A, int proved);

void fmpz_mat_det_modular_given_divisor(fmpz_t det, const fmpz_mat_t A,
        const fmpz_t d, int proved);

void fmpz_mat_det_bound(fmpz_t bound, const fmpz_mat_t A);
void fmpz_mat_det_bound_nonzero(fmpz_t bound, const fmpz_mat_t A);

void fmpz_mat_det_divisor(fmpz_t d, const fmpz_mat_t A);

/* Transforms */

void fmpz_mat_similarity(fmpz_mat_t A, slong r, fmpz_t d);

/* Characteristic polynomial ************************************************/

void _fmpz_mat_charpoly_berkowitz(fmpz * rop, const fmpz_mat_t op);

void fmpz_mat_charpoly_berkowitz(fmpz_poly_t cp,
                                                          const fmpz_mat_t mat);

void _fmpz_mat_charpoly_modular(fmpz * rop, const fmpz_mat_t op);

void fmpz_mat_charpoly_modular(fmpz_poly_t cp, const fmpz_mat_t mat);

FMPZ_MAT_INLINE
void _fmpz_mat_charpoly(fmpz * cp, const fmpz_mat_t mat)
{
   _fmpz_mat_charpoly_modular(cp, mat);
}

FMPZ_MAT_INLINE
void fmpz_mat_charpoly(fmpz_poly_t cp, const fmpz_mat_t mat)
{
   if (mat->r != mat->c)
   {
       flint_throw(FLINT_ERROR, "Exception (nmod_mat_charpoly).  Non-square matrix.\n");
   }

   fmpz_mat_charpoly_modular(cp, mat);
}

/* Characteristic polynomial ************************************************/

slong _fmpz_mat_minpoly_modular(fmpz * rop, const fmpz_mat_t op);

void fmpz_mat_minpoly_modular(fmpz_poly_t cp, const fmpz_mat_t mat);

FMPZ_MAT_INLINE
slong _fmpz_mat_minpoly(fmpz * cp, const fmpz_mat_t mat)
{
   return _fmpz_mat_minpoly_modular(cp, mat);
}

FMPZ_MAT_INLINE
void fmpz_mat_minpoly(fmpz_poly_t cp, const fmpz_mat_t mat)
{
   if (mat->r != mat->c)
   {
       flint_throw(FLINT_ERROR, "Exception (fmpz_mat_minpoly).  Non-square matrix.\n");
   }

   fmpz_mat_minpoly_modular(cp, mat);
}

/* Rank **********************************************************************/

slong fmpz_mat_rank(const fmpz_mat_t A);

/* Nonsingular solving *******************************************************/

void fmpz_mat_solve_bound(fmpz_t N, fmpz_t D,
                                       const fmpz_mat_t A, const fmpz_mat_t B);

int fmpz_mat_solve(fmpz_mat_t X, fmpz_t den,
                                       const fmpz_mat_t A, const fmpz_mat_t B);

int fmpz_mat_solve_cramer(fmpz_mat_t X, fmpz_t den,
                                       const fmpz_mat_t A, const fmpz_mat_t B);

int fmpz_mat_solve_fflu(fmpz_mat_t X, fmpz_t den,
                                       const fmpz_mat_t A, const fmpz_mat_t B);

int fmpz_mat_solve_fflu_precomp(fmpz_mat_t X, const slong * perm,
                                    const fmpz_mat_t FFLU, const fmpz_mat_t B);

mp_limb_t
fmpz_mat_find_good_prime_and_invert(nmod_mat_t Ainv,
		                   const fmpz_mat_t A, const fmpz_t det_bound);

mp_limb_t *
fmpz_mat_dixon_get_crt_primes(slong * num_primes,
		                              const fmpz_mat_t A, mp_limb_t p);

void
_fmpz_mat_solve_dixon(fmpz_mat_t X, fmpz_t mod,
		  const fmpz_mat_t A, const fmpz_mat_t B,
 			       const nmod_mat_t Ainv, mp_limb_t p,
                                               const fmpz_t N, const fmpz_t D);

int fmpz_mat_solve_dixon(fmpz_mat_t X, fmpz_t mod,
        const fmpz_mat_t A, const fmpz_mat_t B);

void
_fmpz_mat_solve_dixon_den(fmpz_mat_t X, fmpz_t den,
                     const fmpz_mat_t A, const fmpz_mat_t B,
                                 const nmod_mat_t Ainv, mp_limb_t p,
                                               const fmpz_t N, const fmpz_t D);

int
fmpz_mat_solve_dixon_den(fmpz_mat_t X, fmpz_t den,
		                       const fmpz_mat_t A, const fmpz_mat_t B);

int
fmpz_mat_solve_multi_mod_den(fmpz_mat_t X, fmpz_t den,
	                               const fmpz_mat_t A, const fmpz_mat_t B);

int
fmpz_mat_can_solve_multi_mod_den(fmpz_mat_t X, fmpz_t den,
                                       const fmpz_mat_t A, const fmpz_mat_t B);

int
fmpz_mat_can_solve_fflu(fmpz_mat_t X, fmpz_t den,
 		                       const fmpz_mat_t A, const fmpz_mat_t B);

int
fmpz_mat_can_solve(fmpz_mat_t X, fmpz_t den,
                                       const fmpz_mat_t A, const fmpz_mat_t B);

/* Nullspace *****************************************************************/

slong fmpz_mat_nullspace(fmpz_mat_t res, const fmpz_mat_t mat);

/* Inverse *******************************************************************/

int fmpz_mat_inv(fmpz_mat_t B, fmpz_t den, const fmpz_mat_t A);

/* Modular reduction and reconstruction **************************************/

void fmpz_mat_set_nmod_mat(fmpz_mat_t A, const nmod_mat_t Amod);
void fmpz_mat_set_nmod_mat_unsigned(fmpz_mat_t A, const nmod_mat_t Amod);
void fmpz_mat_get_nmod_mat(nmod_mat_t Amod, const fmpz_mat_t A);

void fmpz_mat_CRT_ui(fmpz_mat_t res, const fmpz_mat_t mat1, const fmpz_t m1, const nmod_mat_t mat2, int sign);

#ifdef FMPZ_H
void fmpz_mat_multi_mod_ui_precomp(nmod_mat_t * residues, slong nres, const fmpz_mat_t mat, const fmpz_comb_t comb, fmpz_comb_temp_t temp);
void fmpz_mat_multi_CRT_ui_precomp(fmpz_mat_t mat, nmod_mat_t * const residues, slong nres, const fmpz_comb_t comb, fmpz_comb_temp_t temp, int sign);
#endif

void fmpz_mat_multi_mod_ui(nmod_mat_t * residues, slong nres, const fmpz_mat_t mat);
void fmpz_mat_multi_CRT_ui(fmpz_mat_t mat, nmod_mat_t * const residues, slong nres, int sign);

/* HNF and SNF **************************************************************/

void fmpz_mat_hnf(fmpz_mat_t H, const fmpz_mat_t A);
void fmpz_mat_hnf_transform(fmpz_mat_t H, fmpz_mat_t U, const fmpz_mat_t A);
void fmpz_mat_hnf_classical(fmpz_mat_t H, const fmpz_mat_t A);
void fmpz_mat_hnf_xgcd(fmpz_mat_t H, const fmpz_mat_t A);
void fmpz_mat_hnf_minors(fmpz_mat_t H, const fmpz_mat_t A);
void fmpz_mat_hnf_minors_transform(fmpz_mat_t H, fmpz_mat_t U, const fmpz_mat_t A);
void fmpz_mat_hnf_modular(fmpz_mat_t H, const fmpz_mat_t A, const fmpz_t D);
void fmpz_mat_hnf_modular_eldiv(fmpz_mat_t A, const fmpz_t D);
void fmpz_mat_hnf_pernet_stein(fmpz_mat_t H, const fmpz_mat_t A, flint_rand_t state);
int fmpz_mat_is_in_hnf(const fmpz_mat_t A);

void fmpz_mat_snf(fmpz_mat_t S, const fmpz_mat_t A);
void fmpz_mat_snf_diagonal(fmpz_mat_t S, const fmpz_mat_t A);
void fmpz_mat_snf_kannan_bachem(fmpz_mat_t S, const fmpz_mat_t A);
void fmpz_mat_snf_iliopoulos(fmpz_mat_t S, const fmpz_mat_t A,
        const fmpz_t mod);
int fmpz_mat_is_in_snf(const fmpz_mat_t A);

/* Special matrices **********************************************************/

int fmpz_mat_is_hadamard(const fmpz_mat_t A);

int fmpz_mat_hadamard(fmpz_mat_t A);

/* Gram matrix **************************************************************/

void fmpz_mat_gram(fmpz_mat_t B, const fmpz_mat_t A);

int fmpz_mat_is_spd(const fmpz_mat_t A);

/* Conversions **************************************************************/

int fmpz_mat_get_d_mat(d_mat_t B, const fmpz_mat_t A);

int fmpz_mat_get_d_mat_transpose(d_mat_t B, const fmpz_mat_t A);

/* Cholesky Decomposition ****************************************************/

void fmpz_mat_chol_d(d_mat_t R, const fmpz_mat_t A);

/* LLL ***********************************************************************/

int fmpz_mat_is_reduced(const fmpz_mat_t A,
                                                     double delta, double eta);

int fmpz_mat_is_reduced_gram(const fmpz_mat_t A,
                                                     double delta, double eta);

int fmpz_mat_is_reduced_with_removal(const fmpz_mat_t A,
                        double delta, double eta, const fmpz_t gs_B, int newd);

int fmpz_mat_is_reduced_gram_with_removal(const fmpz_mat_t A,
                        double delta, double eta, const fmpz_t gs_B, int newd);

/* Classical LLL *************************************************************/

void fmpz_mat_lll_original(fmpz_mat_t A,
                                         const fmpq_t delta, const fmpq_t eta);

/* Modified LLL **************************************************************/

void fmpz_mat_lll_storjohann(fmpz_mat_t A,
                                         const fmpq_t delta, const fmpq_t eta);

/* Column partitioning *******************************************************/

int fmpz_mat_col_partition(slong * part,
                                              fmpz_mat_t M, int short_circuit);

/* Van Hoeij helper function *************************************************/

int fmpz_mat_next_col_van_hoeij(fmpz_mat_t M, fmpz_t P,
                                       fmpz_mat_t col, slong exp, slong U_exp);

#ifdef __cplusplus
}
#endif

#endif
