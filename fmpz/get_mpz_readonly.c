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

    Copyright (C) 2011 Sebastian Pancratz

******************************************************************************/

#include <mpir.h>
#include "flint.h"
#include "fmpz.h"

__mpz_struct * fmpz_get_mpz_readonly(const fmpz_t f)
{
    if (COEFF_IS_MPZ(*f))
    {
        return COEFF_TO_PTR(*f);
    }
    else
    {
        __mpz_struct *z; 

        if (*f > 0)
        {
            z            = malloc(sizeof(__mpz_struct) + sizeof(mp_limb_t));
            z->_mp_alloc = 1;
            z->_mp_size  = 1;
            z->_mp_d     = (mp_ptr) ((char *) z + sizeof(__mpz_struct));
            z->_mp_d[0]  = *f;
        }
        else if (*f < 0)
        {
            z            = malloc(sizeof(__mpz_struct) + sizeof(mp_limb_t));
            z->_mp_alloc = 1;
            z->_mp_size  = -1;
            z->_mp_d     = (mp_ptr) ((char *) z + sizeof(__mpz_struct));
            z->_mp_d[0]  = -(*f);
        }
        else
        {
            z            = malloc(sizeof(__mpz_struct));
            z->_mp_alloc = 1;
            z->_mp_size  = 0;
        }

        return z;
    }
}

