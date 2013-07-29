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
#include <sstream>
#include <string>

#include "fmpz_matxx.h"
#include "flintxx/test/helpers.h"

using namespace flint;

void
test_init()
{
    // TODO
}

void
test_assignment()
{
    // TODO
}

void
test_conversion()
{
    // TODO
}

template<class Expr>
bool has_explicit_temporaries(const Expr&)
{
    return Expr::ev_traits_t::rule_t::temporaries_t::len != 0;
}
template<class T, class Expr>
bool compare_temporaries(const Expr&)
{
    return mp::equal_types<T,
             typename Expr::ev_traits_t::rule_t::temporaries_t>::val;
}
void
test_arithmetic()
{
    fmpz_matxx A(10, 10);
    fmpz_matxx v(10, 1);
    for(unsigned i = 0;i < 10;++i)
        v.at(i, 0) = i;

    tassert(transpose(v).rows() == 1);
    tassert(transpose(v).cols() == 10);

    tassert(!has_explicit_temporaries(trace(transpose(v))));
    tassert(!has_explicit_temporaries(trace(A + v*transpose(v))));
    tassert(!has_explicit_temporaries(trace((v*transpose(v) + A))));
    tassert(!has_explicit_temporaries(trace(v*transpose(v) + v*transpose(v))));
    tassert((compare_temporaries<tuple<fmpzxx*, empty_tuple> >(
                    ((A+A)*(fmpzxx(1)+fmpzxx(1))))));

    tassert(trace(transpose(v)) == 0);
    tassert(trace(A + v*transpose(v)) == 285);
    tassert(trace(v*transpose(v) + A) == 285);
    tassert(trace(v*transpose(v) + v*transpose(v)) == 2*285);
    tassert(trace((A+A)*(fmpzxx(1) + fmpzxx(1))) == 0);
}

void
test_functions()
{
    // TODO
}

int
main()
{
    std::cout << "fmpz_matxx....";

    test_init();
    test_assignment();
    test_conversion();
    test_arithmetic();
    test_functions();

    std::cout << "PASS" << std::endl;
    return 0;
}
