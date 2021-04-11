#ifndef FMPZ_POLY_ROOTS_H_
#define FMPZ_POLY_ROOTS_H_

#include "fmpz_poly.h"

#include "fq.h"
#include "padic.h"
#include "qadic.h"

typedef struct
{
  fq_struct *x0;
  slong *multiplicity;
  slong num;
  slong alloc;
} fmpz_poly_roots_fq_struct;

typedef fmpz_poly_roots_fq_struct fmpz_poly_roots_fq_t[1];

void fmpz_poly_roots_fq_init2 (fmpz_poly_roots_fq_t roots, slong n, fq_ctx_t fctx);
void fmpz_poly_roots_fq_clear (fmpz_poly_roots_fq_t roots, fq_ctx_t fctx);
void fmpz_poly_roots_fq_print (fmpz_poly_roots_fq_t roots, fq_ctx_t fctx);
char* fmpz_poly_roots_fq_get_str_pretty (fmpz_poly_roots_fq_t roots, fq_ctx_t fctx);
int fmpz_poly_roots_fq_fprint_pretty (FILE *file, fmpz_poly_roots_fq_t roots, fq_ctx_t fctx);
int fmpz_poly_roots_fq_print_pretty (fmpz_poly_roots_fq_t roots, fq_ctx_t fctx);
void fmpz_poly_roots_fq(fmpz_poly_roots_fq_t roots, fmpz_poly_t poly, fq_ctx_t fctx);

typedef struct
{
  padic_struct *x0;
  slong *multiplicity;
  slong num;
  slong alloc;
} fmpz_poly_roots_padic_struct;

typedef fmpz_poly_roots_padic_struct fmpz_poly_roots_padic_t[1];

void fmpz_poly_roots_padic_init2 (fmpz_poly_roots_padic_t roots, slong n);
void fmpz_poly_roots_padic_clear (fmpz_poly_roots_padic_t roots);
char* fmpz_poly_roots_padic_get_str (fmpz_poly_roots_padic_t roots, padic_ctx_t fctx);
int fmpz_poly_roots_padic_fprint (FILE* file, fmpz_poly_roots_padic_t roots, padic_ctx_t fctx);
int fmpz_poly_roots_padic_print (fmpz_poly_roots_padic_t roots, padic_ctx_t fctx);
void fmpz_poly_roots_padic (fmpz_poly_roots_padic_t roots, fmpz_poly_t poly,
			    padic_ctx_t fctx);

typedef struct
{
  qadic_struct *x0;
  slong *multiplicity;
  slong num;
  slong alloc;
} fmpz_poly_roots_qadic_struct;

typedef fmpz_poly_roots_qadic_struct fmpz_poly_roots_qadic_t[1];

void fmpz_poly_roots_qadic_init2 (fmpz_poly_roots_qadic_t roots, slong n);
void fmpz_poly_roots_qadic_clear (fmpz_poly_roots_qadic_t roots);
char* fmpz_poly_roots_qadic_get_str_pretty (fmpz_poly_roots_qadic_t roots, qadic_ctx_t fctx);
int fmpz_poly_roots_qadic_fprint_pretty (FILE * file, fmpz_poly_roots_qadic_t roots, qadic_ctx_t fctx);
int fmpz_poly_roots_qadic_print_pretty (fmpz_poly_roots_qadic_t roots, qadic_ctx_t fctx);
void fmpz_poly_roots_qadic (fmpz_poly_roots_qadic_t roots, fmpz_poly_t poly,
			    qadic_ctx_t fctx);

#endif
