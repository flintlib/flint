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

    Copyright (C) 2013 Tom Bachmann

******************************************************************************/
#include <iostream>
#include <string>
#include <stdexcept>
#include <stdio.h>

#include "cxx.h"
#include "cxx/test/helpers.h"

int
main()
{
    std::cout << "mpz....";
    //flint::mpz a('a');
    flint::mpz a((short)-4);
    flint::mpz b; b = "8";
    std::cout << a + b << std::endl;
    std::cout << "PASS" << std::endl;
}
