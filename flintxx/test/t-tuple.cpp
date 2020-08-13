/*
    Copyright (C) 2013 Tom Bachmann

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <set>

#include "flintxx/test/helpers.h"
#include "flintxx/tuple.h"
#include "flintxx/mp.h"

using namespace flint;
using namespace mp;
using namespace traits;

void
test_make()
{
    typedef make_tuple<int, char, int> make3;
    typedef make_tuple<int, char> make2;
    typedef make_tuple<int> make1;

    tassert((equal_types<
                make3::type,
                tuple<int, tuple<char, tuple<int, empty_tuple> > >
              >::val));
    tassert((equal_types<
                make2::type,
                tuple<int, tuple<char, empty_tuple> >
              >::val));
    tassert((equal_types<make1::type, tuple<int, empty_tuple> >::val));
    tassert(make3::type::len == 3);
    tassert(make2::type::len == 2);
    tassert(make1::type::len == 1);

    make3::type tup = make3::make(1, 'b', 3);
    tassert(tup.head == 1);
    tassert(tup.tail.head == 'b');
    tassert(tup.tail.tail.head == 3);
    tup.first() = 2;
    tassert(tup.head == 2);
    tassert(tup.second() == 'b');

    make2::type pair = make2::make(1, 'b');
    tassert(pair.first() == 1);
    tassert(pair.second() == 'b');

    make1::type singleton = make1::make(1);
    tassert(singleton.first() == 1);

    typedef make_tuple<const int&, const int&> makeref;
    tassert((equal_types<makeref::type::head_t, const int&>::val));
    int a = 4;
    makeref::type reftup = makeref::make(a, a);
    tassert(reftup.head == 4);
    a = 5;
    tassert(reftup.second() == 5);
}

template<int i, bool default_ctor = false>
struct type_n
{
    int payload;
    type_n(int p) : payload(p) {};
};
template<int i>
struct type_n<i, true>
{
    int payload;
    type_n(int p) : payload(p) {};
    type_n() {};
};

void
test_back()
{
    typedef type_n<0> A;
    typedef type_n<1> B;
    typedef type_n<2> C;
    typedef type_n<3> D;

    typedef make_tuple<A, B, C> make3;
    typedef make_tuple<A*, B*, C*> make3p;
    typedef make_tuple<B, C> make2A;
    typedef make_tuple<A, C> make2B;

    typedef back_tuple<make3p::type, void> backvoid;
    typedef back_tuple<make3p::type, D> backD;
    typedef back_tuple<make3p::type, A> backA;
    typedef back_tuple<make3p::type, B> backB;

    tassert((equal_types<backvoid::type, make3::type>::val));
    tassert((equal_types<backD::type, make3::type>::val));
    tassert((equal_types<backA::type, make2A::type>::val));
    tassert((equal_types<backB::type, make2B::type>::val));

    make3::type backing = make3::make(1, 2, 3);
    make3p::type pointers;
    backvoid::init(pointers, backing, 0);
    tassert(pointers.first()->payload == 1);
    backing.second() = 4;
    tassert(pointers.second()->payload == 4);

    A ret = 4;
    make2A::type backing2 = make2A::make(1, 2);
    backA::init(pointers, backing2, &ret);
    tassert(pointers.first()->payload == 4);
    tassert(pointers.second()->payload == 1);
    ret.payload = 5;
    backing2.head.payload = 6;
    tassert(pointers.first()->payload == 5);
    tassert(pointers.second()->payload == 6);

    // Test some corner cases.
    tassert((equal_types<back_tuple<empty_tuple>::type, empty_tuple>::val));
    typedef tuple<int*, empty_tuple> tup1_t;
    tassert((equal_types<back_tuple<tup1_t, int>::type, empty_tuple>::val));
}

void
test_concat()
{
    typedef type_n<0> A;
    typedef type_n<1> B;
    typedef type_n<2> C;
    typedef type_n<3> D;
    typedef type_n<4> E;

    typedef make_tuple<A, B, C> make1st;
    typedef make_tuple<D, E> make2nd;

    typedef concat_tuple<make1st::type, make2nd::type> concater;
    tassert((equal_types<
                concater::type,
                tuple<A, tuple<B, tuple<C, make2nd::type> > >
              >::val));

    concater::type concated = concater::type(0,
            tuple<B, tuple<C, make2nd::type> >(1,
                tuple<C, make2nd::type>(2, make2nd::make(3, 4))));
    make2nd::type second = concater::get_second(concated);
    tassert(second.first().payload == 3);
    tassert(second.second().payload == 4);
    make1st::type first = concater::get_first(concated);
    tassert(first.first().payload == 0);
    tassert(first.second().payload == 1);
    tassert(first.tail.second().payload == 2);

    first.first().payload = 17;
    second.first().payload = 18;
    concated = concater::doit(first, second);
    tassert(concated.first().payload == 17);
    tassert(concated.tail.tail.tail.head.payload == 18);

    // Test a few corner cases.
    tassert((equal_types<
                concat_tuple<empty_tuple, empty_tuple>::type,
                empty_tuple
              >::val));
    typedef tuple<int, empty_tuple> tup1_t;
    tassert((equal_types<
                concat_tuple<tup1_t, empty_tuple>::type,
                tup1_t
              >::val));
    tassert((equal_types<
                concat_tuple<empty_tuple, tup1_t>::type,
                tup1_t
              >::val));
}

namespace flint{
namespace detail
{
template<class T>
struct fillit
{
    static void doit(T& t, int start)
    {
        t.head.payload = start;
        detail::fillit<typename T::tail_t>::doit(t.tail, ++start);
    }
};

template<>
struct fillit<empty_tuple>
{
    static void doit(empty_tuple e, int start)
    {
    }
};

template<class T>
struct values
{
    static void doit(const T& t, std::set<int>& values)
    {
        values.insert(t.head.payload);
        detail::values<typename T::tail_t>::doit(t.tail, values);
    }
};

template<>
struct values<empty_tuple>
{
    static void doit(empty_tuple e, std::set<int>& values)
    {
    }
};
}
}

template<class T>
void fillit(T& t, int start = 1)
{
    detail::fillit<T>::doit(t, start);
}

template<class T>
std::set<int> values(const T& t)
{
    std::set<int> ret;
    detail::values<T>::doit(t, ret);
    return ret;
}

void
test_merge()
{
    typedef type_n<0, true> A;
    typedef type_n<1, true> B;
    typedef type_n<2, true> C;
    typedef type_n<3, true> D;
    typedef type_n<4, true> E;

    // Test that merging tuples <T1, T2, T3> and <U1, U2>
    // yields something of length n, and that the constitutents can be
    // recovered. Also do other order.
#define TEST32n(T1, T2, T3, U1, U2, n) \
    { \
        typedef make_tuple<T1, T2, T3> m1; \
        typedef make_tuple<U1, U2> m2; \
        typedef merge_tuple<m1::type, m2::type> merge1; \
        typedef merge_tuple<m2::type, m1::type> merge2; \
        tassert(merge1::type::len == n); \
        tassert(merge2::type::len == n); \
        merge1::type merged1; \
        { \
            fillit(merged1); \
            m1::type n1 = merge1::get_first(merged1); \
            m2::type n2 = merge1::get_second(merged1); \
            std::set<int> vals1 = values(n1), vals2 = values(n2); \
            vals1.insert(vals2.begin(), vals2.end()); \
            tassert(vals1.size() == n); \
        } \
        merge2::type merged2; \
        { \
            fillit(merged2); \
            m2::type n2 = merge2::get_first(merged2); \
            m1::type n1 = merge2::get_second(merged2); \
            std::set<int> vals1 = values(n1), vals2 = values(n2); \
            vals1.insert(vals2.begin(), vals2.end()); \
            tassert(vals1.size() == n); \
        } \
    }

    TEST32n(A, B, C, A, B, 3);
    TEST32n(A, B, C, B, A, 3);
    TEST32n(A, B, C, A, C, 3);
    TEST32n(A, B, C, C, A, 3);
    TEST32n(A, B, C, C, B, 3);
    TEST32n(A, B, C, B, C, 3);
    TEST32n(A, A, A, A, A, 3);

    TEST32n(A, B, C, A, D, 4);
    TEST32n(A, B, C, B, D, 4);
    TEST32n(A, B, C, C, D, 4);
    TEST32n(A, B, C, D, A, 4);
    TEST32n(A, B, C, D, B, 4);
    TEST32n(A, B, C, D, C, 4);
    TEST32n(A, A, A, A, B, 4);

    TEST32n(A, B, C, D, E, 5);
    TEST32n(A, A, A, B, B, 5);

    typedef merge_tuple<empty_tuple, empty_tuple> emptymerge;
    tassert((equal_types<emptymerge::type, empty_tuple>::val));

    typedef merge_tuple<make_tuple<A, A>::type, empty_tuple> anothermerge;
    tassert(anothermerge::type::len == 2);
}

void
test_traits()
{
    tassert((is_homogeneous_tuple<make_tuple<int, int, int>::type, int>::val
                == true));
    tassert((is_homogeneous_tuple<make_tuple<int, char, int>::type, int>::val
                == false));

    tassert(is_tuple<empty_tuple>::val == true);
    tassert(is_tuple<make_tuple<int>::type>::val == true);
    tassert((is_tuple<make_tuple<int, slong>::type>::val == true));
    tassert(is_tuple<int>::val == false);
}

struct filler
{
    template<class T>
    T create() const
    {
        if(mp::equal_types<T, slong>::val)
            return 1;
        return 0;
    }
};
void
test_htuples()
{
    tassert((equal_types<make_tuple<int, int, int>::type,
                         make_homogeneous_tuple<int, 3>::type>::val));

    typedef make_tuple<int, int, int> maker;
    maker::type tpl = maker::make(1, 2, 3);

    tassert(htuples::extract<1>(tpl) == make_tuple<int>::make(1));
    tassert((htuples::extract<2>(tpl) == make_tuple<int, int>::make(1, 2)));
    tassert((htuples::extract<3>(tpl) == make_tuple<int, int, int>::make(1, 2, 3)));

    typedef make_tuple<int, int> make2;
    tassert(htuples::removeres(tpl, 1) == make2::make(2, 3));
    tassert(htuples::removeres(tpl, 2) == make2::make(1, 3));
    tassert(htuples::removeres(tpl, 3) == make2::make(1, 2));
    tassert(htuples::removeres(tpl, 4) == make2::make(1, 2));

    typedef make_tuple<int, slong, int> maker3;
    tassert((maker3::make(0, 1, 0) == htuples::fill<maker3::type>(filler())));
}

void
test_get()
{
    typedef make_tuple<int, int&, const int&>::type tuple_t;
    int a, b;
    tuple_t t(17, make_tuple<int&, const int&>::type(a,
                make_tuple<const int&>::type(b, empty_tuple())));

    tassert((equal_types<int, tuple_get<tuple_t, 0>::type>::val));
    tassert((equal_types<int&, tuple_get<tuple_t, 1>::type>::val));
    tassert((equal_types<const int&, tuple_get<tuple_t, 2>::type>::val));

    tassert((tuple_get<tuple_t, 0>::get(t) == 17));
    tuple_get<tuple_t, 1>::get(t) = 42;
    b = 3;
    tassert((tuple_get<tuple_t, 2>::get(t) == 3));
    tassert(a == 42);

    const tuple_t ct(t);
    tuple_get<tuple_t, 1>::get(ct) = 17;
    tassert((tuple_get<tuple_t, 1>::get(ct) == 17));
}

void
test_set()
{
    make_tuple<int, int, int>::type tup;
    tup.set(make_tuple<int, char, short>::make(1, 2, 3));
    tassert((tup == make_tuple<int, int, int>::make(1, 2, 3)));
}

void
test_equals_elementwise()
{
    tassert((make_tuple<int, int, int>::make(1, 2, 3).equals_elementwise(
                    make_tuple<int, char, slong>::make(1, 2, 3))));
    tassert((!make_tuple<int, int, int>::make(0, 2, 3).equals_elementwise(
                    make_tuple<int, char, slong>::make(1, 2, 3))));
}

int
main()
{
    std::cout << "tuple....";

    test_make();
    test_back();
    test_concat();
    test_merge();
    test_traits();
    test_htuples();
    test_get();
    test_set();
    test_equals_elementwise();

    std::cout << "PASS" << std::endl;
    return 0;
}

