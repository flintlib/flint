/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2013 Qingwen GUAN

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include "fmpz.h"

/* printing *******************************************************************/

int fmpz_fprint(FILE * file, const fmpz_t x)
{
	if (!COEFF_IS_MPZ(*x))
        return fprintf(file, WORD_FMT "d", *x);
	else
        return (int) mpz_out_str(file, 10, COEFF_TO_PTR(*x));
}

int fmpz_print(const fmpz_t x) { return fmpz_fprint(stdout, x); }

/* reading ********************************************************************/

int
fmpz_fread(FILE * file, fmpz_t f)
{
    mpz_t t;
    size_t r;

    mpz_init(t);
    r = mpz_inp_str(t, file, 10);
    fmpz_set_mpz(f, t);
    mpz_clear(t);

    return (r > 0) ? 1 : 0;
}

int fmpz_read(fmpz_t f) { return fmpz_fread(stdout, f); }

/* file I/O ********************************************************************/

size_t fmpz_inp_raw(fmpz_t x, FILE * fin)
{
    mpz_t v;
    size_t size;

    mpz_init( v );
    size = mpz_inp_raw( v, fin );
    fmpz_set_mpz( x, v );
    mpz_clear( v );

    return size;
}

size_t fmpz_out_raw(FILE * fout, const fmpz_t x)
{
    mpz_t v;
    size_t size;

    mpz_init( v );
    fmpz_get_mpz( v, x );
    size = mpz_out_raw( fout, v );
    mpz_clear( v );

    return size;
}
