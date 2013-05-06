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

    Copyright (C) 2011 Fredrik Johansson

******************************************************************************/

#include <math.h>
#include <gmp.h>
#include <mpfr.h>
#include "flint.h"
#include "fmpz.h"
#include "arith.h"

/* This nice round number precisely fits on 32 bits */
#define NUMBER_OF_SMALL_PARTITIONS 128

const unsigned int
partitions_lookup[NUMBER_OF_SMALL_PARTITIONS] = 
{
    1UL,1UL,2UL,3UL,5UL,7UL,11UL,15UL,22UL,30UL,42UL,56UL,77UL,101UL,135UL,
    176UL,231UL,297UL,385UL,490UL,627UL,792UL,1002UL,1255UL,1575UL,1958UL,
    2436UL,3010UL,3718UL,4565UL,5604UL,6842UL,8349UL,10143UL,12310UL,14883UL,
    17977UL,21637UL,26015UL,31185UL,37338UL,44583UL,53174UL,63261UL,75175UL,
    89134UL,105558UL,124754UL,147273UL,173525UL,204226UL,239943UL,281589UL,
    329931UL,386155UL,451276UL,526823UL,614154UL,715220UL,831820UL,966467UL,
    1121505UL,1300156UL,1505499UL,1741630UL,2012558UL,2323520UL,2679689UL,
    3087735UL,3554345UL,4087968UL,4697205UL,5392783UL,6185689UL,7089500UL,
    8118264UL,9289091UL,10619863UL,12132164UL,13848650UL,15796476UL,18004327UL,
    20506255UL,23338469UL,26543660UL,30167357UL,34262962UL,38887673UL,
    44108109UL,49995925UL,56634173UL,64112359UL,72533807UL,82010177UL,
    92669720UL,104651419UL,118114304UL,133230930UL,150198136UL,169229875UL,
    190569292UL,214481126UL,241265379UL,271248950UL,304801365UL,342325709UL,
    384276336UL,431149389UL,483502844UL,541946240UL,607163746UL,679903203UL,
    761002156UL,851376628UL,952050665UL,1064144451UL,1188908248UL,1327710076UL,
    1482074143UL,1653668665UL,1844349560UL,2056148051UL,2291320912UL,
    2552338241UL,2841940500UL,3163127352UL,3519222692UL,3913864295UL
};

void
arith_number_of_partitions(fmpz_t x, ulong n)
{
    if (n < NUMBER_OF_SMALL_PARTITIONS)
    {
        fmpz_set_ui(x, partitions_lookup[n]);
    }
    else
    {
        mpfr_t t;
        mpfr_init(t);
        arith_number_of_partitions_mpfr(t, n);
        mpfr_get_z(_fmpz_promote(x), t, MPFR_RNDN);
        _fmpz_demote_val(x);
        mpfr_clear(t);
    }
}
