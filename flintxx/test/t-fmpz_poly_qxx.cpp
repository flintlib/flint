/*
    Copyright (C) 2013 Tom Bachmann

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <iostream>
#include <sstream>
#include <string>

#include "fmpz_poly_qxx.h"

#include "flintxx/test/helpers.h"

using namespace flint;

void
test_manipulation()
{
    fmpz_poly_qxx f;
    f.num() = "3  0 1 1";
    f.den() = "3  0 -1 1";
    tassert(!f.is_canonical());
    f.canonicalise();
    tassert(f.num().to_string() == "2  1 1");
    tassert(f.den().to_string() == "2  -1 1");

    tassert(fmpz_poly_qxx::zero().is_zero());
    tassert(fmpz_poly_qxx::one().is_one());
}

void
test_assignment_conversion()
{
    fmpz_poly_qxx f, g;
    f = 1;
    tassert(f.is_one());
    g = 0;
    f = g;
    tassert(f.is_zero());
    tassert(f.to_string() == "0");

    f = "4  1 0 0 1";
    tassert(f.num().to_string() == "4  1 0 0 1");
    tassert(f.den().to_string() == "1  1");

    f.den() = "2  -1 1";
    tassert(f.to_string() == "4  1 0 0 1/2  -1 1");
    g = "4  1 0 0 1/2  -1 1";
    tassert(f == g);
    tassert(f == fmpz_poly_qxx("4  1 0 0 1/2  -1 1"));

    tassert(f.pretty("x") == "(x^3+1)/(x-1)");
}

void
test_arithmetic()
{
    fmpz_poly_qxx f, g;
    g = "4  1 0 0 1/2  -1 1";
    f = "1  1";

    tassert((f + g).to_string() == "4  0 1 0 1/2  -1 1");
    tassert((g - f).to_string() == "4  2 -1 0 1/2  -1 1");
    tassert(g - f == g + (-f));
    tassert(inv(g).to_string() == "2  -1 1/4  1 0 0 1");

    tassert(2 * g == g * 2);
    f = 2*g;
    tassert(f.num() == 2*g.num() && f.den() == g.den());
    f /= 2;
    tassert(f == g);

    f = "1  1";
    tassert((f*g).num() == f.num()*g.num());
    tassert((f*g).den() == f.den()*g.den());
    tassert((f/g).num() == f.num()*g.den());
    tassert((f/g).den() == f.den()*g.num());
}

// Won't compile if the expression is not done using addmul
template<class T>
bool is_ternary(const T&)
{
    return T::ev_traits_t::temp_rule_t::TERNARY_OP_MARKER + 1;
}

// test stuff which we should get automatically - addmul, references etc
void
test_extras()
{
    // TODO
}

void
test_functions()
{
    fmpz_poly_qxx f;
    tassert(f.is_zero() && !f.is_one());
    f.num() = 1;
    tassert(f.is_one());

    f = "4  1 0 0 1/2  -1 1";
    tassert(pow(f, 4u) == f*f*f*f);

    tassert(derivative(f).to_string() == "4  -1 0 -3 2/3  1 -2 1");

    // test static methods
    frandxx rand;
    tassert(fmpz_poly_qxx::randtest(rand, 10, 8, 10, 8).num().degree() < 10);
    tassert(fmpz_poly_qxx::randtest(rand, 10, 8, 10, 8).den().degree() < 10);
    tassert(flog(height(fmpz_poly_qxx::randtest(
                        rand, 10, 8, 10, 8).num()), 2u) < 8);
    tassert(flog(height(fmpz_poly_qxx::randtest(
                        rand, 10, 8, 10, 8).den()), 2u) < 8);
    tassert(!fmpz_poly_qxx::randtest_not_zero(rand, 10, 8, 10, 8).is_zero());

    tassert(f.derivative() == derivative(f));
    tassert(f.inv() == inv(f));
    tassert(f.pow(7u) == pow(f, 7u));
}

void
test_printing()
{
    if(0)
    {
        // make sure these compile
        fmpz_poly_qxx f;
        print(f);
        print_pretty(f, "x");
    }
}

void
test_unified_access()
{
    fmpz_poly_qxx a;
    const fmpz_poly_qxx& b = a;
    tassert(b.num().is_zero() && b.den().is_one());
}

int
main()
{
    std::cout << "fmpz_poly_qxx....";

    test_manipulation();
    test_assignment_conversion();
    test_arithmetic();
    test_functions();
    test_extras();
    test_printing();
    test_unified_access();

    std::cout << "PASS" << std::endl;
    return 0;
}

