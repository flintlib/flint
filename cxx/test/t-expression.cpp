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

#include<sstream>

#include "cxx/expression.h"
#include "cxx/test/helpers.h"

using namespace flint;

template<class Operation, class Data>
class my_expression
    : public expression<derived_wrapper<my_expression>, Operation, Data>
{
public:
    my_expression() {};
    template<class T>
    explicit my_expression(const T& t) : my_expression::expression(t) {}

    template<class T>
    my_expression& operator=(const T& t)
    {
        this->set(t);
        return *this;
    }

protected:
    explicit my_expression(const Data& d) : my_expression::expression(d) {}

    template<class D, class O, class Da>
    friend class flint::expression;
};

struct data
{
    int payload;
    bool* destroyed;
    int extra;

    data() : payload(-1), destroyed(0), extra(0) {}

    ~data()
    {
        if(destroyed)
            *destroyed = true;
    }

    data(const data& d)
    {
        payload = d.payload;
        destroyed = 0;
        extra = d.extra | 1 << 30;
    }
};

typedef my_expression<operations::immediate, data> myint;

namespace flint {
namespace rules {
template<>
struct empty_initialization<myint>
{
    static void doit(myint& v)
    {
        v._data().extra = 42;
    }
};

template<>
struct initialization<myint, myint>
{
    static void doit(myint& to, const myint& from)
    {
        to._data().payload = from._data().payload;
        to._data().extra = 1;
    }
};

template<>
struct initialization<myint, int>
{
    static void doit(myint& to, int i)
    {
        to._data().payload = i;
        to._data().extra = 2;
    }
};

template<>
struct initialization<myint, char>
{
    static void doit(myint& to, char i)
    {
        to._data().payload = i;
        to._data().extra = 3;
    }
};

template<>
struct print<myint>
{
    static void doit(const myint& v, std::ostream& o)
    {
        o << v._data().payload;
    }
};

template<>
struct assignment<myint, myint>
{
    static void doit(myint& to, const myint& from)
    {
        to._data().payload = from._data().payload;
        to._data().extra = 4;
    }
};

template<>
struct assignment<myint, long>
{
    static void doit(myint& to, long from)
    {
        to._data().payload = from;
        to._data().extra = 5;
    }
};
} // rules
} // flint

void
test_initialization()
{
    myint a;
    tassert(a._data().payload == -1 && a._data().extra == 42);

    a._data().payload = 17;
    myint b = a;
    tassert(b._data().payload == 17 && b._data().extra == (1 | 1 << 30));

    myint c(8);
    tassert(c._data().payload == 8 && c._data().extra == 2);

    myint d('a');
    tassert(d._data().payload == 'a' && d._data().extra == 3);
}

void
test_destruction()
{
    bool destroyed = false;
    {
        myint a;
        a._data().destroyed = &destroyed;
    }
    tassert(destroyed);
}

void
test_printing()
{
    std::ostringstream o;
    myint a(4);
    o << a;
    tassert(o.str() == "4");
}

void
test_assignment()
{
    myint a;
    myint b(43);
    a = b;
    tassert(a._data().payload == 43);
    // tassert(a._data().extra == 4); // this won't happen because of copy init

    a = 44l;
    tassert(a._data().payload == 44 && a._data().extra == 5);
}

int
main()
{
    test_initialization();
    test_destruction();
    test_printing();
    test_assignment();

    // TODO test that certain things *don't* compile?

    return 0;
}
