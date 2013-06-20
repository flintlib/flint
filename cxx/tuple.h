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

#ifndef CXX_TUPLE_H
#define CXX_TUPLE_H

#include "cxx/traits.h"

// TODO
// * tuple merging

namespace flint {
///////////////////////////
// A general tuple class.
///////////////////////////
// This simply stores a head and a tail. Conventionally, tuples a built by
// nesting the tail. The last entry should be an empty_tuple (see below).
// So e.g. a pair of integers would be tuple<int, tuple<int, empty_tuple> >.
//
// There are some helpers in the mp namespace below.
template<class Head, class Tail>
struct tuple
{
    Head head;
    Tail tail;

    typedef Head head_t;

    // Obtain a reference to head (convenience name for pairs).
    typename traits::reference<head_t>::type first() {return head;}
    typename traits::forwarding<head_t>::type first() const {return head;}

    // Obtain a reference to the head of the tail (convenience name for pairs).
    typename traits::reference<typename Tail::head_t>::type
        second() {return tail.head;}
    typename traits::forwarding<typename Tail::head_t>::type
        second() const {return tail.head;}

    tuple() {};
    tuple(typename traits::forwarding<Head>::type h,
          typename traits::forwarding<Tail>::type t)
        : head(h), tail(t)
    {
    }
};

template<>
struct tuple<void, void>
{
    struct empty { };
    typedef empty head_t;
    typedef empty tail_t;
    empty head;
    empty tail;
};
typedef tuple<void, void> empty_tuple;

namespace mp {

// Helper to conveniently define tuple types, and marshall objects into
// tuples.
// Typical usage:
// typedef make_tuple<int, char, long> maker;
// maker::type my_tuple = maker::make(1, 'a', 2l);
template<class T1, class T2 = void, class T3 = void>
struct make_tuple
{
    typedef tuple<T1, tuple<T2, tuple<T3, empty_tuple> > > type;
    static type make(typename traits::forwarding<T1>::type t1,
              typename traits::forwarding<T2>::type t2,
              typename traits::forwarding<T3>::type t3)
    {
        return type(t1, make_tuple<T2, T3>::make(t2, t3));
    }
};
template<class T1, class T2>
struct make_tuple<T1, T2, void>
{
    typedef tuple<T1, tuple<T2, empty_tuple> > type;
    static type make(typename traits::forwarding<T1>::type t1,
              typename traits::forwarding<T2>::type t2)
    {
        return type(t1, make_tuple<T2>::make(t2));
    }
};
template<class T1>
struct make_tuple<T1, void, void>
{
    typedef tuple<T1, empty_tuple> type;
    static type make(typename traits::forwarding<T1>::type t1)
    {
        return type(t1, empty_tuple());
    }
};

// Create a tuple backing a tuple of points.
// That is to say, given a tuple like <A*, B*, C*, D*>,
// compute a backing tuple type <A, B, C, D>.
// 
// If Tuple::head_t is the same as Return*, do not back the head
// and instead feed it in separately. I.e. if Return is A*, then type
// will be just <B*, C*, D*>.
//
// The static member init(to, from, ret) can be used to initalize the tuple
// of pointers "to" to point to its backing in "from" and "ret".
template<class Tuple, class Return = void>
struct back_tuple;

// General case: non-empty tuple <Head, Tail>, and Return type cannot be
// merged in.
template<class Head, class Tail, class Return>
struct back_tuple<tuple<Head*, Tail>, Return>
{
    typedef tuple<Head, typename back_tuple<Tail, empty_tuple>::type> type;
    static void init(tuple<Head*, Tail>& to, type& from, Return* ret /* unused */)
    {
        back_tuple<Tail, Return>::init(to.tail, from.tail, ret);
        to.head = &from.head;
    }
};

// Merging case: non-empty tuple <Head, Tail>, and Return is Head*
template<class Head, class Tail>
struct back_tuple<tuple<Head*, Tail>, Head>
{
    typedef typename back_tuple<Tail, void>::type type;
    static void init(tuple<Head*, Tail>& to, type& from, Head* ret)
    {
        to.head = ret;
        back_tuple<Tail, void>::init(to.tail, from, ret /* unused */ );
    }
};

// Base case: empty tuple; nothing to do.
template<class Return>
struct back_tuple<empty_tuple, Return>
{
    typedef empty_tuple type;
    static void init(empty_tuple& to, type& from, Return* ret)
    {
    }
};

// Helper to concatenate two tuples.
//
// This has one member type, and two static member functions get_second and
// get_first.
// "type" is a tuple type which can store the data of both Tuple1 and Tuple2.
// Then, given an element of type "type", get_first and get_second can be used
// to fill types Tuple1 and Tuple2. Note that get_second can return a constant
// reference, whereas get_first has to do copying.
// (But these copies should usually be inlined and optimized away by the
// compiler.)
//
// Example: Tuple1 = <A, B, C>, Tuple2 = <D, E>.
// Then type = <A, B, C, D, E>
//   get_second(t) = t.tail.tail.tail
//   get_first(t) = tuple(t.head, tuple(t.tail.head, tuple(
//                                   t.tail.tail.head, empty_tuple())));
//
template<class Tuple1, class Tuple2>
struct concat_tuple;

// Degenerate case: Tuple1 is empty.
template<class Tuple2>
struct concat_tuple<empty_tuple, Tuple2>
{
    typedef Tuple2 type;
    static const Tuple2& get_second(const Tuple2& t) {return t;}
    static empty_tuple get_first(const Tuple2& t) {return empty_tuple();}
};
// General case: non-empty Tuple1.
template<class Head, class Tail, class Tuple2>
struct concat_tuple<tuple<Head, Tail>, Tuple2>
{
    typedef tuple<Head, typename concat_tuple<Tail, Tuple2>::type> type;
    static const Tuple2& get_second(const type& t)
    {
        return concat_tuple<Tail, Tuple2>::get_second(t.tail);
    }
    static tuple<Head, Tail> get_first(const type& o)
    {
        return tuple<Head, Tail>(
                o.head,
                concat_tuple<Tail, Tuple2>::get_first(o.tail));
    }
};
} // mp
} // flint

#endif
