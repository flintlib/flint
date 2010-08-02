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

    Copyright (C) 2009 William Hart
    Copyright (C) 2010 Sebastian Pancratz

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <mpir.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result;
    printf("content....");
    fflush(stdout);

    _fmpz_vec_randinit();

    // Check that content(a f) = abs(a) content(f)
    for (i = 0; i < 10000; i++)
    {
        fmpz_t a, c, d;
        fmpz *f;
        long len = n_randint(100);

        fmpz_init(a);
        fmpz_init(c);
        fmpz_init(d);
        f = _fmpz_vec_init(len);
        _fmpz_vec_randtest(f, len, n_randint(200));
        fmpz_randtest(a, n_randint(100));

        _fmpz_vec_content(c, f, len);
        _fmpz_vec_scalar_mul_fmpz(f, f, len, a);
        fmpz_abs(a, a);
        fmpz_mul(c, a, c);
        _fmpz_vec_content(d, f, len);

        result = (fmpz_equal(c, d));
        if (!result)
        {
            printf("FAIL:\n");
            fmpz_print(c), printf("\n\n");
            fmpz_print(d), printf("\n\n");
            abort();
        }

        fmpz_clear(a);
        fmpz_clear(c);
        fmpz_clear(d);
        _fmpz_vec_clear(f, len);
    }

    _fmpz_vec_randclear();
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
