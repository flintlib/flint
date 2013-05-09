/*============================================================================
    Copyright 2006 Jason Papadopoulos.    
    Copyright 2006, 2011 William Hart.

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

===============================================================================

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	
       				   --jasonp@boo.net 9/8/06
       				   
The following modifications were made by William Hart:
    -added the utility function get_null_entry
    -reformatted original code so it would operate as a standalone 
     filter and block Lanczos module
--------------------------------------------------------------------*/


#undef ulong /* avoid clash with stdlib */
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#define ulong unsigned long 

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "qsieve.h"

#define BIT(x) (((uint64_t)(1)) << (x))

static const uint64_t bitmask[64] = {
	BIT( 0), BIT( 1), BIT( 2), BIT( 3), BIT( 4), BIT( 5), BIT( 6), BIT( 7),
	BIT( 8), BIT( 9), BIT(10), BIT(11), BIT(12), BIT(13), BIT(14), BIT(15),
	BIT(16), BIT(17), BIT(18), BIT(19), BIT(20), BIT(21), BIT(22), BIT(23),
	BIT(24), BIT(25), BIT(26), BIT(27), BIT(28), BIT(29), BIT(30), BIT(31),
	BIT(32), BIT(33), BIT(34), BIT(35), BIT(36), BIT(37), BIT(38), BIT(39),
	BIT(40), BIT(41), BIT(42), BIT(43), BIT(44), BIT(45), BIT(46), BIT(47),
	BIT(48), BIT(49), BIT(50), BIT(51), BIT(52), BIT(53), BIT(54), BIT(55),
	BIT(56), BIT(57), BIT(58), BIT(59), BIT(60), BIT(61), BIT(62), BIT(63),
};

/*--------------------------------------------------------------------*/
uint64_t get_null_entry(uint64_t * nullrows, len_t i, len_t l) {
   
   /* Returns true if the entry with indices i,l is 1 in the
      supplied 64xN matrix. This is used to read the nullspace
      vectors which are output by the Lanczos routine */
      
    return nullrows[i]&bitmask[l];
}

/*--------------------------------------------------------------------*/
void reduce_matrix(qs_t qs_inf, len_t *nrows, len_t *ncols, la_col_t *cols) {

	/* Perform light filtering on the nrows x ncols
	   matrix specified by cols[]. The processing here is
	   limited to deleting columns that contain a singleton
	   row, then resizing the matrix to have a few more
	   columns than rows. Because deleting a column reduces
	   the counts in several different rows, the process
	   must iterate to convergence.
	   
	   Note that this step is not intended to make the Lanczos
	   iteration run any faster (though it will); it's just
	   that if we don't go to this trouble then there are 
	   factorizations for which the matrix step will fail 
	   outright  */

	len_t r, c, i, j, k;
	len_t passes;
	len_t *counts;
	len_t reduced_rows;
	len_t reduced_cols;

	/* count the number of nonzero entries in each row */

	counts = (len_t *)flint_calloc((size_t)*nrows, sizeof(len_t));
	for (i = 0; i < *ncols; i++) {
		for (j = 0; j < cols[i].weight; j++)
			counts[cols[i].data[j]]++;
	}

	reduced_rows = *nrows;
	reduced_cols = *ncols;
	passes = 0;

	do {
		r = reduced_rows;

		/* remove any columns that contain the only entry
		   in one or more rows, then update the row counts
		   to reflect the missing column. Iterate until
		   no more columns can be deleted */

		do {
			c = reduced_cols;
			for (i = j = 0; i < reduced_cols; i++) {
				la_col_t *col = cols + i;
				for (k = 0; k < col->weight; k++) {
					if (counts[col->data[k]] < 2)
						break;
				}
	
				if (k < col->weight) {
					for (k = 0; k < col->weight; k++) {
						counts[col->data[k]]--;
					}
					free_col(col);
				   clear_col(col);
				}
				else {
					cols[j++] = cols[i];
					if (j-1 != i) clear_col(col);
				}
			}
			reduced_cols = j;
		} while (c != reduced_cols);
	
		/* count the number of rows that contain a
		   nonzero entry */

		for (i = reduced_rows = 0; i < *nrows; i++) {
			if (counts[i])
				reduced_rows++;
		}

		/* Because deleting a column reduces the weight
		   of many rows, the number of nonzero rows may
		   be much less than the number of columns. Delete
		   more columns until the matrix has the correct
		   aspect ratio. Columns at the end of cols[] are
		   the heaviest, so delete those (and update the
		   row counts again) */

		if (reduced_cols > reduced_rows + qs_inf->extra_rels) {
			for (i = reduced_rows + qs_inf->extra_rels;
					i < reduced_cols; i++) {

				la_col_t *col = cols + i;
				for (j = 0; j < col->weight; j++) {
					counts[col->data[j]]--;
				}
				free_col(col);
				clear_col(col);
			}
			reduced_cols = reduced_rows + qs_inf->extra_rels;
		}

		/* if any columns were deleted in the previous step,
		   then the matrix is less dense and more columns
		   can be deleted; iterate until no further deletions
		   are possible */

		passes++;

	} while (r != reduced_rows);

#if (QS_DEBUG & 128)
	printf("reduce to %ld x %ld in %ld passes\n", 
			reduced_rows, reduced_cols, passes);
#endif

	flint_free(counts);
    
	/* record the final matrix size. Note that we can't touch
	   nrows because all the column data (and the sieving relations
	   that produced it) would have to be updated */

	*ncols = reduced_cols;
}

/*-------------------------------------------------------------------*/
static void mul_64x64_64x64(uint64_t *a, uint64_t *b, uint64_t *c ) {

	/* c[][] = x[][] * y[][], where all operands are 64 x 64
	   (i.e. contain 64 words of 64 bits each). The result
	   may overwrite a or b. */

	uint64_t ai, bj, accum;
	uint64_t tmp[64];
	unsigned long i, j;

	for (i = 0; i < 64; i++) {
		j = 0;
		accum = 0;
		ai = a[i];

		while (ai) {
			bj = b[j];
			if( ai & 1 )
				accum ^= bj;
			ai >>= 1;
			j++;
		}

		tmp[i] = accum;
	}
	memcpy(c, tmp, sizeof(tmp));
}

/*-----------------------------------------------------------------------*/
static void precompute_Nx64_64x64(uint64_t *x, uint64_t *c) {

	/* Let x[][] be a 64 x 64 matrix in GF(2), represented
	   as 64 words of 64 bits each. Let c[][] be an 8 x 256
	   matrix of 64-bit words. This code fills c[][] with
	   a bunch of "partial matrix multiplies". For 0<=i<256,
	   the j_th row of c[][] contains the matrix product

	   	( i << (8*j) ) * x[][]

	   where the quantity in parentheses is considered a 
	   1 x 64 vector of elements in GF(2). The resulting
	   table can dramatically speed up matrix multiplies
	   by x[][]. */

	uint64_t accum, xk;
	unsigned long i, j, k, index;

	for (j = 0; j < 8; j++) {
		for (i = 0; i < 256; i++) {
			k = 0;
			index = i;
			accum = 0;
			while (index) {
				xk = x[k];
				if (index & 1)
					accum ^= xk;
				index >>= 1;
				k++;
			}
			c[i] = accum;
		}

		x += 8;
		c += 256;
	}
}

/*-------------------------------------------------------------------*/
static void mul_Nx64_64x64_acc(uint64_t *v, uint64_t *x, uint64_t *c, 
				uint64_t *y, len_t n) {

	/* let v[][] be a n x 64 matrix with elements in GF(2), 
	   represented as an array of n 64-bit words. Let c[][]
	   be an 8 x 256 scratch matrix of 64-bit words.
	   This code multiplies v[][] by the 64x64 matrix 
	   x[][], then XORs the n x 64 result into y[][] */

    len_t i;
	uint64_t word;

	precompute_Nx64_64x64(x, c);

	for (i = 0; i < n; i++) {
		word = v[i];
		y[i] ^=  c[ 0*256 + ((word>> 0) & 0xff) ]
		       ^ c[ 1*256 + ((word>> 8) & 0xff) ]
		       ^ c[ 2*256 + ((word>>16) & 0xff) ]
		       ^ c[ 3*256 + ((word>>24) & 0xff) ]
		       ^ c[ 4*256 + ((word>>32) & 0xff) ]
		       ^ c[ 5*256 + ((word>>40) & 0xff) ]
		       ^ c[ 6*256 + ((word>>48) & 0xff) ]
		       ^ c[ 7*256 + ((word>>56)       ) ];
	}
}

/*-------------------------------------------------------------------*/
static void mul_64xN_Nx64(uint64_t *x, uint64_t *y,
			   uint64_t *c, uint64_t *xy, len_t n) {

	/* Let x and y be n x 64 matrices. This routine computes
	   the 64 x 64 matrix xy[][] given by transpose(x) * y.
	   c[][] is a 256 x 8 scratch matrix of 64-bit words. */

	len_t i;

	memset(c, 0, 256 * 8 * sizeof(uint64_t));
	memset(xy, 0, 64 * sizeof(uint64_t));

	for (i = 0; i < n; i++) {
		uint64_t xi = x[i];
		uint64_t yi = y[i];
		c[ 0*256 + ( xi        & 0xff) ] ^= yi;
		c[ 1*256 + ((xi >>  8) & 0xff) ] ^= yi;
		c[ 2*256 + ((xi >> 16) & 0xff) ] ^= yi;
		c[ 3*256 + ((xi >> 24) & 0xff) ] ^= yi;
		c[ 4*256 + ((xi >> 32) & 0xff) ] ^= yi;
		c[ 5*256 + ((xi >> 40) & 0xff) ] ^= yi;
		c[ 6*256 + ((xi >> 48) & 0xff) ] ^= yi;
		c[ 7*256 + ((xi >> 56)       ) ] ^= yi;
	}


	for(i = 0; i < 8; i++) {

		unsigned long j;
		uint64_t a0, a1, a2, a3, a4, a5, a6, a7;

		a0 = a1 = a2 = a3 = 0;
		a4 = a5 = a6 = a7 = 0;

		for (j = 0; j < 256; j++) {
			if ((j >> i) & 1) {
				a0 ^= c[0*256 + j];
				a1 ^= c[1*256 + j];
				a2 ^= c[2*256 + j];
				a3 ^= c[3*256 + j];
				a4 ^= c[4*256 + j];
				a5 ^= c[5*256 + j];
				a6 ^= c[6*256 + j];
				a7 ^= c[7*256 + j];
			}
		}

		xy[ 0] = a0; xy[ 8] = a1; xy[16] = a2; xy[24] = a3;
		xy[32] = a4; xy[40] = a5; xy[48] = a6; xy[56] = a7;
		xy++;
	}
}

/*-------------------------------------------------------------------*/
static len_t find_nonsingular_sub(uint64_t *t, len_t *s, 
				len_t *last_s, len_t last_dim, 
				uint64_t *w) {

	/* given a 64x64 matrix t[][] (i.e. sixty-four
	   64-bit words) and a list of 'last_dim' column 
	   indices enumerated in last_s[]: 
	   
	     - find a submatrix of t that is invertible 
	     - invert it and copy to w[][]
	     - enumerate in s[] the columns represented in w[][] */

	len_t i, j;
	len_t dim;
	len_t cols[64];
	uint64_t M[64][2];
	uint64_t mask, *row_i, *row_j;
	uint64_t m0, m1;

	/* M = [t | I] for I the 64x64 identity matrix */

	for (i = 0; i < 64; i++) {
		M[i][0] = t[i]; 
		M[i][1] = bitmask[i];
	}

	/* put the column indices from last_s[] into the
	   back of cols[], and copy to the beginning of cols[]
	   any column indices not in last_s[] */

	mask = 0;
	for (i = 0; i < last_dim; i++) {
		cols[63 - i] = last_s[i];
		mask |= bitmask[last_s[i]];
	}
	for (i = j = 0; i < 64; i++) {
		if (!(mask & bitmask[i]))
			cols[j++] = i;
	}

	/* compute the inverse of t[][] */

	for (i = dim = 0; i < 64; i++) {
	
		/* find the next pivot row and put in row i */

		mask = bitmask[cols[i]];
		row_i = M[cols[i]];

		for (j = i; j < 64; j++) {
			row_j = M[cols[j]];
			if (row_j[0] & mask) {
				m0 = row_j[0];
				m1 = row_j[1];
				row_j[0] = row_i[0];
				row_j[1] = row_i[1];
				row_i[0] = m0; 
				row_i[1] = m1;
				break;
			}
		}
				
		/* if a pivot row was found, eliminate the pivot
		   column from all other rows */

		if (j < 64) {
			for (j = 0; j < 64; j++) {
				row_j = M[cols[j]];
				if ((row_i != row_j) && (row_j[0] & mask)) {
					row_j[0] ^= row_i[0];
					row_j[1] ^= row_i[1];
				}
			}

			/* add the pivot column to the list of 
			   accepted columns */

			s[dim++] = cols[i];
			continue;
		}

		/* otherwise, use the right-hand half of M[]
		   to compensate for the absence of a pivot column */

		for (j = i; j < 64; j++) {
			row_j = M[cols[j]];
			if (row_j[1] & mask) {
				m0 = row_j[0];
				m1 = row_j[1];
				row_j[0] = row_i[0];
				row_j[1] = row_i[1];
				row_i[0] = m0; 
				row_i[1] = m1;
				break;
			}
		}
				
		if (j == 64) {
#if (QS_DEBUG & 128)
			printf("lanczos error: submatrix "
					"is not invertible\n");
#endif
			return 0;
		}
			
		/* eliminate the pivot column from the other rows
		   of the inverse */

		for (j = 0; j < 64; j++) {
			row_j = M[cols[j]];
			if ((row_i != row_j) && (row_j[1] & mask)) {
				row_j[0] ^= row_i[0];
				row_j[1] ^= row_i[1];
			}
		}

		/* wipe out the pivot row */

		row_i[0] = row_i[1] = 0;
	}

	/* the right-hand half of M[] is the desired inverse */
	
	for (i = 0; i < 64; i++) 
		w[i] = M[i][1];

	/* The block Lanczos recurrence depends on all columns
	   of t[][] appearing in s[] and/or last_s[]. 
	   Verify that condition here */

	mask = 0;
	for (i = 0; i < dim; i++)
		mask |= bitmask[s[i]];
	for (i = 0; i < last_dim; i++)
		mask |= bitmask[last_s[i]];

	if (mask != (uint64_t)(-1)) {
#if (QS_DEBUG & 128)
		printf("lanczos error: not all columns used\n");
#endif
		return 0;
	}

	return dim;
}

/*-------------------------------------------------------------------*/
void mul_MxN_Nx64(len_t vsize, len_t dense_rows,
		len_t ncols, la_col_t *A,
		uint64_t *x, uint64_t *b) {

	/* Multiply the vector x[] by the matrix A (stored
	   columnwise) and put the result in b[]. vsize
	   refers to the number of uint64_t's allocated for
	   x[] and b[]; vsize is probably different from ncols */

	len_t i, j;

	memset(b, 0, vsize * sizeof(uint64_t));
	
	for (i = 0; i < ncols; i++) {
		la_col_t *col = A + i;
		len_t *row_entries = col->data;
		uint64_t tmp = x[i];

		for (j = 0; j < col->weight; j++) {
			b[row_entries[j]] ^= tmp;
		}
	}

	if (dense_rows) {
		for (i = 0; i < ncols; i++) {
			la_col_t *col = A + i;
			len_t *row_entries = col->data + col->weight;
			uint64_t tmp = x[i];
	
			for (j = 0; j < dense_rows; j++) {
				if (row_entries[j / 32] & 
						((len_t)1 << (j % 32))) {
					b[j] ^= tmp;
				}
			}
		}
	}
}

/*-------------------------------------------------------------------*/
void mul_trans_MxN_Nx64(len_t dense_rows, len_t ncols,
			la_col_t *A, uint64_t *x, uint64_t *b) {

	/* Multiply the vector x[] by the transpose of the
	   matrix A and put the result in b[]. Since A is stored
	   by columns, this is just a matrix-vector product */

	len_t i, j;

	for (i = 0; i < ncols; i++) {
		la_col_t *col = A + i;
		len_t *row_entries = col->data;
		uint64_t accum = 0;

		for (j = 0; j < col->weight; j++) {
			accum ^= x[row_entries[j]];
		}
		b[i] = accum;
	}

	if (dense_rows) {
		for (i = 0; i < ncols; i++) {
			la_col_t *col = A + i;
			len_t *row_entries = col->data + col->weight;
			uint64_t accum = b[i];
	
			for (j = 0; j < dense_rows; j++) {
				if (row_entries[j / 32] &
						((len_t)1 << (j % 32))) {
					accum ^= x[j];
				}
			}
			b[i] = accum;
		}
	}
}

/*-----------------------------------------------------------------------*/
static void transpose_vector(len_t ncols, uint64_t *v, uint64_t **trans) {

	/* Hideously inefficent routine to transpose a
	   vector v[] of 64-bit words into a 2-D array
	   trans[][] of 64-bit words */

	len_t i, j;
	len_t col;
	uint64_t mask, word;

	for (i = 0; i < ncols; i++) {
		col = i / 64;
		mask = bitmask[i % 64];
		word = v[i];
		j = 0;
		while (word) {
			if (word & 1)
				trans[j][col] |= mask;
			word = word >> 1;
			j++;
		}
	}
}

/*-----------------------------------------------------------------------*/
void combine_cols(len_t ncols, 
		uint64_t *x, uint64_t *v, 
		uint64_t *ax, uint64_t *av) {

	/* Once the block Lanczos iteration has finished, 
	   x[] and v[] will contain mostly nullspace vectors
	   between them, as well as possibly some columns
	   that are linear combinations of nullspace vectors.
	   Given vectors ax[] and av[] that are the result of
	   multiplying x[] and v[] by the matrix, this routine 
	   will use Gauss elimination on the columns of [ax | av] 
	   to find all of the linearly dependent columns. The
	   column operations needed to accomplish this are mir-
	   rored in [x | v] and the columns that are independent
	   are skipped. Finally, the dependent columns are copied
	   back into x[] and represent the nullspace vector output
	   of the block Lanczos code.
	   
	   v[] and av[] can be NULL, in which case the elimination
	   process assumes 64 dependencies instead of 128 */

	len_t i, j, k, bitpos, col, col_words, num_deps;
	uint64_t mask;
	uint64_t *matrix[128], *amatrix[128], *tmp;

	num_deps = 128;
	if (v == NULL || av == NULL)
		num_deps = 64;

	col_words = (ncols + 63) / 64;

	for (i = 0; i < num_deps; i++) {
		matrix[i] = (uint64_t *)flint_calloc((size_t)col_words, 
					     sizeof(uint64_t));
		amatrix[i] = (uint64_t *)flint_calloc((size_t)col_words, 
					      sizeof(uint64_t));
	}

	/* operations on columns can more conveniently become 
	   operations on rows if all the vectors are first
	   transposed */

	transpose_vector(ncols, x, matrix);
	transpose_vector(ncols, ax, amatrix);
	if (num_deps == 128) {
		transpose_vector(ncols, v, matrix + 64);
		transpose_vector(ncols, av, amatrix + 64);
	}

	/* Keep eliminating rows until the unprocessed part
	   of amatrix[][] is all zero. The rows where this
	   happens correspond to linearly dependent vectors
	   in the nullspace */

	for (i = bitpos = 0; i < num_deps && bitpos < ncols; bitpos++) {

		/* find the next pivot row */

		mask = bitmask[bitpos % 64];
		col = bitpos / 64;
		for (j = i; j < num_deps; j++) {
			if (amatrix[j][col] & mask) {
				tmp = matrix[i];
				matrix[i] = matrix[j];
				matrix[j] = tmp;
				tmp = amatrix[i];
				amatrix[i] = amatrix[j];
				amatrix[j] = tmp;
				break;
			}
		}
		if (j == num_deps)
			continue;

		/* a pivot was found; eliminate it from the
		   remaining rows */

		for (j++; j < num_deps; j++) {
			if (amatrix[j][col] & mask) {

				/* Note that the entire row, *not*
				   just the nonzero part of it, must
				   be eliminated; this is because the
				   corresponding (dense) row of matrix[][]
				   must have the same operation applied */

				for (k = 0; k < col_words; k++) {
					amatrix[j][k] ^= amatrix[i][k];
					matrix[j][k] ^= matrix[i][k];
				}
			}
		}
		i++;
	}

	/* transpose rows i to 64 back into x[] */

	for (j = 0; j < ncols; j++) {
		uint64_t word = 0;

		col = j / 64;
		mask = bitmask[j % 64];

		for (k = i; k < 64; k++) {
			if (matrix[k][col] & mask)
				word |= bitmask[k];
		}
		x[j] = word;
	}

	for (i = 0; i < num_deps; i++) {
		flint_free(matrix[i]);
		flint_free(amatrix[i]);
	}
}

/*-----------------------------------------------------------------------*/
uint64_t * block_lanczos(flint_rand_t state, len_t nrows, 
			len_t dense_rows, len_t ncols, la_col_t *B) {
	
	/* Solve Bx = 0 for some nonzero x; the computed
	   solution, containing up to 64 of these nullspace
	   vectors, is returned */

	uint64_t *vnext, *v[3], *x, *v0;
	uint64_t *winv[3];
	uint64_t *vt_a_v[2], *vt_a2_v[2];
	uint64_t *scratch;
	uint64_t *d, *e, *f, *f2;
	uint64_t *tmp;
	len_t s[2][64];
	len_t i, iter;
	len_t n = ncols;
	len_t dim0, dim1;
	uint64_t mask0, mask1;
	len_t vsize;

	/* allocate all of the size-n variables. Note that because
	   B has been preprocessed to ignore singleton rows, the
	   number of rows may really be less than nrows and may
	   be greater than ncols. vsize is the maximum of these
	   two numbers  */

	vsize = FLINT_MAX(nrows, ncols);
	v[0] = (uint64_t *)flint_malloc(vsize * sizeof(uint64_t));
	v[1] = (uint64_t *)flint_malloc(vsize * sizeof(uint64_t));
	v[2] = (uint64_t *)flint_malloc(vsize * sizeof(uint64_t));
	vnext = (uint64_t *)flint_malloc(vsize * sizeof(uint64_t));
	x = (uint64_t *)flint_malloc(vsize * sizeof(uint64_t));
	v0 = (uint64_t *)flint_malloc(vsize * sizeof(uint64_t));
	scratch = (uint64_t *)flint_malloc(FLINT_MAX(vsize, 256 * 8) * sizeof(uint64_t));

	/* allocate all the 64x64 variables */

	winv[0] = (uint64_t *)flint_malloc(64 * sizeof(uint64_t));
	winv[1] = (uint64_t *)flint_malloc(64 * sizeof(uint64_t));
	winv[2] = (uint64_t *)flint_malloc(64 * sizeof(uint64_t));
	vt_a_v[0] = (uint64_t *)flint_malloc(64 * sizeof(uint64_t));
	vt_a_v[1] = (uint64_t *)flint_malloc(64 * sizeof(uint64_t));
	vt_a2_v[0] = (uint64_t *)flint_malloc(64 * sizeof(uint64_t));
	vt_a2_v[1] = (uint64_t *)flint_malloc(64 * sizeof(uint64_t));
	d = (uint64_t *)flint_malloc(64 * sizeof(uint64_t));
	e = (uint64_t *)flint_malloc(64 * sizeof(uint64_t));
	f = (uint64_t *)flint_malloc(64 * sizeof(uint64_t));
	f2 = (uint64_t *)flint_malloc(64 * sizeof(uint64_t));

	/* The iterations computes v[0], vt_a_v[0],
	   vt_a2_v[0], s[0] and winv[0]. Subscripts larger
	   than zero represent past versions of these
	   quantities, which start off empty (except for
	   the past version of s[], which contains all
	   the column indices */
	   
	memset(v[1], 0, vsize * sizeof(uint64_t));
	memset(v[2], 0, vsize * sizeof(uint64_t));
	for (i = 0; i < 64; i++) {
		s[1][i] = i;
		vt_a_v[1][i] = 0;
		vt_a2_v[1][i] = 0;
		winv[1][i] = 0;
		winv[2][i] = 0;
	}
	dim0 = 0;
	dim1 = 64;
	mask1 = (uint64_t)(-1);
	iter = 0;

	/* The computed solution 'x' starts off random,
	   and v[0] starts off as B*x. This initial copy
	   of v[0] must be saved off separately */

	for (i = 0; i < n; i++)
#if FLINT_BITS==64
		v[0][i] = (uint64_t) n_randlimb(state);
#else
		v[0][i] = (uint64_t) n_randlimb(state) + ((uint64_t) n_randlimb(state) << 32);
#endif

	memcpy(x, v[0], vsize * sizeof(uint64_t));
	mul_MxN_Nx64(vsize, dense_rows, ncols, B, v[0], scratch);
	mul_trans_MxN_Nx64(dense_rows, ncols, B, scratch, v[0]);
	memcpy(v0, v[0], vsize * sizeof(uint64_t));

	/* perform the iteration */

	while (1) {
		iter++;

		/* multiply the current v[0] by a symmetrized
		   version of B, or B'B (apostrophe means 
		   transpose). Use "A" to refer to B'B  */

		mul_MxN_Nx64(vsize, dense_rows, ncols, B, v[0], scratch);
		mul_trans_MxN_Nx64(dense_rows, ncols, B, scratch, vnext);

		/* compute v0'*A*v0 and (A*v0)'(A*v0) */

		mul_64xN_Nx64(v[0], vnext, scratch, vt_a_v[0], n);
		mul_64xN_Nx64(vnext, vnext, scratch, vt_a2_v[0], n);

		/* if the former is orthogonal to itself, then
		   the iteration has finished */

		for (i = 0; i < 64; i++) {
			if (vt_a_v[0][i] != 0)
				break;
		}
		if (i == 64) {
			break;
		}

		/* Find the size-'dim0' nonsingular submatrix
		   of v0'*A*v0, invert it, and list the column
		   indices present in the submatrix */

		dim0 = find_nonsingular_sub(vt_a_v[0], s[0], 
					    s[1], dim1, winv[0]);
		if (dim0 == 0)
			break;

		/* mask0 contains one set bit for every column
		   that participates in the inverted submatrix
		   computed above */

		mask0 = 0;
		for (i = 0; i < dim0; i++)
			mask0 |= bitmask[s[0][i]];

		/* compute d */

		for (i = 0; i < 64; i++)
			d[i] = (vt_a2_v[0][i] & mask0) ^ vt_a_v[0][i];

		mul_64x64_64x64(winv[0], d, d);

		for (i = 0; i < 64; i++)
			d[i] = d[i] ^ bitmask[i];

		/* compute e */

		mul_64x64_64x64(winv[1], vt_a_v[0], e);

		for (i = 0; i < 64; i++)
			e[i] = e[i] & mask0;

		/* compute f */

		mul_64x64_64x64(vt_a_v[1], winv[1], f);

		for (i = 0; i < 64; i++)
			f[i] = f[i] ^ bitmask[i];

		mul_64x64_64x64(winv[2], f, f);

		for (i = 0; i < 64; i++)
			f2[i] = ((vt_a2_v[1][i] & mask1) ^ 
				   vt_a_v[1][i]) & mask0;

		mul_64x64_64x64(f, f2, f);

		/* compute the next v */

		for (i = 0; i < n; i++)
			vnext[i] = vnext[i] & mask0;

		mul_Nx64_64x64_acc(v[0], d, scratch, vnext, n);
		mul_Nx64_64x64_acc(v[1], e, scratch, vnext, n);
		mul_Nx64_64x64_acc(v[2], f, scratch, vnext, n);
		
		/* update the computed solution 'x' */

		mul_64xN_Nx64(v[0], v0, scratch, d, n);
		mul_64x64_64x64(winv[0], d, d);
		mul_Nx64_64x64_acc(v[0], d, scratch, x, n);

		/* rotate all the variables */

		tmp = v[2]; 
		v[2] = v[1]; 
		v[1] = v[0]; 
		v[0] = vnext; 
		vnext = tmp;
		
		tmp = winv[2]; 
		winv[2] = winv[1]; 
		winv[1] = winv[0]; 
		winv[0] = tmp;
		
		tmp = vt_a_v[1]; vt_a_v[1] = vt_a_v[0]; vt_a_v[0] = tmp;
		
		tmp = vt_a2_v[1]; vt_a2_v[1] = vt_a2_v[0]; vt_a2_v[0] = tmp;

		memcpy(s[1], s[0], 64 * sizeof(len_t));
		mask1 = mask0;
		dim1 = dim0;
	}

#if (QS_DEBUG & 128)
	printf("lanczos halted after %ld iterations\n", iter);
#endif

	/* free unneeded storage */

    
    flint_free(vnext);
	flint_free(scratch);
	flint_free(v0);
	flint_free(vt_a_v[0]);
	flint_free(vt_a_v[1]);
	flint_free(vt_a2_v[0]);
	flint_free(vt_a2_v[1]);
	flint_free(winv[0]);
	flint_free(winv[1]);
	flint_free(winv[2]);
	flint_free(d);
	flint_free(e);
	flint_free(f);
	flint_free(f2);

	/* if a recoverable failure occurred, start everything
	   over again */

	if (dim0 == 0) {
#if (QS_DEBUG & 128)
		printf("linear algebra failed; retrying...\n");
#endif
		flint_free(x);
		flint_free(v[0]);
		flint_free(v[1]);
		flint_free(v[2]);
		return NULL;
	}

	/* convert the output of the iteration to an actual
	   collection of nullspace vectors */

	mul_MxN_Nx64(vsize, dense_rows, ncols, B, x, v[1]);
	mul_MxN_Nx64(vsize, dense_rows, ncols, B, v[0], v[2]);

	combine_cols(ncols, x, v[0], v[1], v[2]);

	/* verify that these really are linear dependencies of B */

	mul_MxN_Nx64(vsize, dense_rows, ncols, B, x, v[0]);
	
	for (i = 0; i < ncols; i++) {
		if (v[0][i] != 0)
			break;
	}
	if (i < ncols) {
		printf("lanczos error: dependencies don't work %ld\n",i);
		abort();
	}
	
	flint_free(v[0]);
	flint_free(v[1]);
	flint_free(v[2]);
	return x;
}
