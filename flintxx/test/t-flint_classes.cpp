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
#include "flintxx/test/helpers.h"
#include "fmpzxx.h"

using namespace flint;

void
test_traits()
{
    using namespace flint_classes;
    tassert((mp::equal_types<to_nonref<fmpzxx>::type, fmpzxx>::val));
    tassert((mp::equal_types<to_nonref<fmpzxx_ref>::type, fmpzxx>::val));
    tassert((mp::equal_types<to_nonref<fmpzxx_srcref>::type, fmpzxx>::val));
    tassert((mp::equal_types<to_ref<fmpzxx>::type, fmpzxx_ref>::val));
    tassert((mp::equal_types<to_srcref<fmpzxx>::type, fmpzxx_srcref>::val));
    tassert((mp::equal_types<to_ref<fmpzxx_ref>::type, fmpzxx_ref>::val));
    tassert((is_ref<fmpzxx, fmpzxx_ref>::val));
    tassert((is_ref<fmpzxx_ref, fmpzxx_ref>::val));
    tassert((!is_ref<fmpzxx, fmpzxx>::val));
    tassert((!is_ref<fmpzxx, fmpzxx_srcref>::val));
    tassert((is_srcref<fmpzxx, fmpzxx_srcref>::val));
    tassert((!is_srcref<fmpzxx, fmpzxx>::val));
    tassert((!is_srcref<fmpzxx, fmpzxx_ref>::val));
    tassert(is_nonref<fmpzxx>::val);
    tassert(!is_nonref<fmpzxx_srcref>::val);
    tassert(!is_nonref<fmpzxx_ref>::val);
}


fmpzxx_ref make_ref(fmpzxx& a) {return fmpzxx_ref(a);}
void
test_ref_lvalues()
{
    fmpzxx a(1);
    fmpzxx_ref b(a);
    b *= 2;
    make_ref(a) *= 2;
    tassert(a == 4);
}

int
main()
{
    std::cout << "flint_classes....";

    test_traits();
    test_ref_lvalues();

    std::cout << "PASS" << std::endl;
}
