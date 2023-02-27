/*
    Copyright (C) 2013 Tom Bachmann

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <sstream>

#include "flintxx/vector.h"

#include "flintxx/test/helpers.h"
#include "flintxx/test/myint.h"
#include "fmpz_vecxx.h"

using namespace flint;

template<class Vec, class Tuple>
struct tuple_has
    : mp::or_<
          mp::equal_types<Vec, typename Tuple::head_t>,
          tuple_has<Vec, typename Tuple::tail_t> > { };
template<class Vec>
struct tuple_has<Vec, empty_tuple> : mp::false_ { };

template<class Vec, class Expr>
bool has_vector_temporaries(const Expr& e)
{
    return tuple_has<Vec*, typename Expr::ev_traits_t::rule_t::temporaries_t>::val;
}

template<class Vec>
void test(const Vec& original, const char* str = "(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)")
{
    Vec v(original);
    Vec w(original);
    Vec u(original);
    Vec x(original);
    for(slong i = 0;i < (slong) v.size();++i)
    {
        v[i] = i;
        w[i] = WORD(0);
        u[i] = 2*i;
        x[i] = 8*i;
    }
    tassert(v.to_string() == str);
    tassert(v == v);
    tassert(v != w);

    tassert(u == v + v);
    tassert(v == v + w);
    tassert(x == ((v + v) + (v + v)) + ((v + v) + (v + v)));
    tassert((v + w)[0] == v[0] + w[0]);

    if(!mp::equal_types<Vec, fmpz_vecxx>::val)
        tassert(!has_vector_temporaries<Vec>(
                    ((v + v) + (v + v)) + ((v + v) + (v + v))));
}

// b/c generic vector has no assignment implementation ...
void
test_fmpz_vecxx()
{
    fmpz_vecxx a(10), b(10);
    a[0] = 1;
    tassert(a != b);
    a = b;
    tassert(a == b);
}

int
main()
{
    std::cout << "vector....";

    typedef make_vector<myint>::type intvec;
    typedef make_vector_n<myint, 10>::type intvec10;

    test(intvec(10));
    test(intvec10());
    test(fmpz_vecxx(10), "10  0 1 2 3 4 5 6 7 8 9");
    test_fmpz_vecxx();

    std::cout << "PASS" << std::endl;

    return 0;
}
