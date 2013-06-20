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
#ifndef CXX_TEST_HELPERS_H
#define CXX_TEST_HELPERS_H

#include <iostream>
#include <cstdlib>

#define tassert(expr) do \
{ \
    if (!(expr)) \
    { \
        std::cout << "FAIL\n" __FILE__ ":" << __LINE__ << ": assertion " #expr " failed\n"; \
        std::exit(1); \
    } \
} while (0)

#endif
