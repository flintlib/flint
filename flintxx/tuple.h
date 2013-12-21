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

#include "traits.h"

namespace flint {
///////////////////////////
// A general tuple class.
///////////////////////////
// This simply stores a head and a tail. Conventionally, tuples a built by
// nesting the tail. The last entry should be an empty_tuple (see below).
// So e.g. a pair of integers would be tuple<int, tuple<int, empty_tuple> >.
//
// There are some helpers in the mp namespace below.
namespace detail {
struct FILLIT { };
} // detail
namespace mp {
namespace htuples {
// Filling of tuples
template<class Tuple, class Filler>
Tuple fill(const Filler& f)
{
    return Tuple(detail::FILLIT(), f);
}
} // htuples
} // mp

template<class Head, class Tail>
struct tuple
{
    Head head;
    Tail tail;

    typedef Head head_t;
    typedef Tail tail_t;
    static const unsigned len = 1 + Tail::len;

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

    bool operator==(const tuple& other)
    {
        return head == other.head && tail == other.tail;
    }

    template<class T>
    void set(const T& t)
    {
        head = t.head;
        tail.set(t.tail);
    }
    template<class T>
    void set(T& t)
    {
        head = t.head;
        tail.set(t.tail);
    }

    template<class T>
    bool equals_elementwise(const T& t) const
    {
        return head == t.head && tail.equals_elementwise(t.tail);
    }

private:
    template<class Filler>
    tuple(detail::FILLIT fillit, const Filler& f)
        : head(f.template create<head_t>()), tail(fillit, f)
    {
    }

    template<class Tuple, class Filler>
    friend Tuple mp::htuples::fill(const Filler&);
    template<class H, class T>
    friend struct tuple;
};

struct empty_tuple
{
    struct empty { };
    typedef empty head_t;
    typedef empty tail_t;
    empty head;
    empty tail;
    static const unsigned len = 0;

    bool operator==(const empty_tuple&) {return true;}

    empty_tuple() {}

    void set(empty_tuple) {}
    bool equals_elementwise(empty_tuple) const {return true;}

private:
    template<class Filler>
    empty_tuple(detail::FILLIT fillit, const Filler& f)
    {
    }

    template<class Tuple, class Filler>
    friend Tuple mp::htuples::fill(const Filler&);
    template<class Head, class Tail>
    friend struct tuple;
};

namespace detail {
typedef void UNUSED;

template<class T>
struct maybe_forwarding
{
    typedef typename traits::forwarding<T>::type type;
    static type default_value();
};
template<>
struct maybe_forwarding<UNUSED>
{
    typedef UNUSED* type;
    static type default_value() {return 0;}
};
}
namespace mp {
// Helper to conveniently define tuple types, and marshall objects into
// tuples.
// Typical usage:
// typedef make_tuple<int, char, slong> maker;
// maker::type my_tuple = maker::make(1, 'a', WORD(2));
// TODO this would be a prime use for variadic templates ...
#define FLINTXX_MAKE_TUPLE_TEMPLATE_ARGS \
    class T1 = detail::UNUSED, class T2 = detail::UNUSED, \
    class T3 = detail::UNUSED, class T4 = detail::UNUSED, \
    class T5 = detail::UNUSED, class T6 = detail::UNUSED, \
    class T7 = detail::UNUSED, class T8 = detail::UNUSED
#define FLINTXX_MAKE_TUPLE_FUNC_ARGS \
    typename detail::maybe_forwarding<T1>::type t1 \
        = detail::maybe_forwarding<T1>::default_value(), \
    typename detail::maybe_forwarding<T2>::type t2 \
        = detail::maybe_forwarding<T2>::default_value(), \
    typename detail::maybe_forwarding<T3>::type t3 \
        = detail::maybe_forwarding<T3>::default_value(), \
    typename detail::maybe_forwarding<T4>::type t4 \
        = detail::maybe_forwarding<T4>::default_value(), \
    typename detail::maybe_forwarding<T5>::type t5 \
        = detail::maybe_forwarding<T5>::default_value(), \
    typename detail::maybe_forwarding<T6>::type t6 \
        = detail::maybe_forwarding<T6>::default_value(), \
    typename detail::maybe_forwarding<T7>::type t7 \
        = detail::maybe_forwarding<T7>::default_value(), \
    typename detail::maybe_forwarding<T8>::type t8 \
        = detail::maybe_forwarding<T8>::default_value()
#define FLINTXX_MAKE_TUPLE_TYPES_APPLYMACRO(func) \
    func(T1), func(T2), func(T3), func(T4), \
    func(T5), func(T6), func(T7), func(T8)
#define FLINTXX_MAKE_TUPLE_FUNC_ARG_NAMES t1, t2, t3, t4, t5, t6, t7, t8
template<FLINTXX_MAKE_TUPLE_TEMPLATE_ARGS>
struct make_tuple
{
    typedef make_tuple<T2, T3, T4, T5, T6, T7, T8> next;
    typedef typename next::type tail_t;
    typedef tuple<T1, tail_t> type;
    static type make(FLINTXX_MAKE_TUPLE_FUNC_ARGS)
    {
        return type(t1, next::make(t2, t3, t4, t5, t6, t7, t8));
    }
};
template<>
struct make_tuple<detail::UNUSED, detail::UNUSED, detail::UNUSED,
    detail::UNUSED, detail::UNUSED, detail::UNUSED, detail::UNUSED,
    detail::UNUSED>
{
    typedef detail::UNUSED* T;
    typedef empty_tuple type;
    // g++-4.4 bolts if we use (...), even though all arguments are PODs
    static empty_tuple make(T=0, T=0, T=0, T=0, T=0, T=0, T=0, T=0)
        {return empty_tuple();}
};


// Indexified access
template<class Tuple, unsigned idx>
struct tuple_get
{
    typedef tuple_get<typename Tuple::tail_t, idx-1> nget;
    typedef typename nget::type type;

    static typename traits::forwarding<type>::type get(const Tuple& t)
    {
        return nget::get(t.tail);
    }
    static typename traits::reference<type>::type get(Tuple& t)
    {
        return nget::get(t.tail);
    }
};

template<class Tuple>
struct tuple_get<Tuple, 0>
{
    typedef typename Tuple::head_t type;
    static typename traits::forwarding<type>::type get(const Tuple& t)
    {
        return t.head;
    }
    static typename traits::reference<type>::type get(Tuple& t)
    {
        return t.head;
    }
};

// Create a tuple backing a tuple of points.
// That is to say, given a tuple like <A*, B*, C*, D*>,
// compute a backing tuple type <A, B, C, D>.
// 
// If one of the types in the tuple is Return*, do not back it
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
    typedef tuple<Head, typename back_tuple<Tail, Return>::type> type;
    static void init(tuple<Head*, Tail>& to, type& from, Return* ret = 0)
    {
        back_tuple<Tail, Return>::init(to.tail, from.tail, ret);
        to.head = &from.head;
    }
};

// Merging case: non-empty tuple <Head, Tail>, and Return is Head*
template<class Head, class Tail>
struct back_tuple<tuple<Head*, Tail>, Head>
{
    typedef typename back_tuple<Tail, void /* no more merging */>::type type;
    static void init(tuple<Head*, Tail>& to, type& from, Head* ret = 0)
    {
        to.head = ret;
        back_tuple<Tail, void>::init(to.tail, from, 0 /* unused */ );
    }
};

// Base case: empty tuple; nothing to do.
template<class Return>
struct back_tuple<empty_tuple, Return>
{
    typedef empty_tuple type;
    static void init(empty_tuple& to, type& from, Return* ret = 0)
    {
    }
};

// Helper to concatenate two tuples.
//
// This has one member type, and three static member functions get_second,
// get_first and doit.
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

    static type doit(const empty_tuple& t1, const Tuple2& t2) {return t2;}
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
    static type doit(const tuple<Head, Tail>& t1, const Tuple2& t2)
    {
        return type(t1.head, concat_tuple<Tail, Tuple2>::doit(t1.tail, t2));
    }
};


// Merging of tuples
//
// This takes two tuples, and computes a tuple type which can store either.
// As usual, the extraction functions require copying which can be amortized
// by the compiler.
//
template<class Tuple1, class Tuple2>
struct merge_tuple;
//{
//    typedef XYZ type;
//    Tuple1 get_first(const type& type);
//    Tuple2 get_second(const type& type);
//};

// General case
// NB: tail is *always* a sub-tuple of the second argument!
template<class Head, class Tail, class Tuple2>
struct merge_tuple<tuple<Head, Tail>, Tuple2>
{
public:
    typedef merge_tuple<Tail, Tuple2> comp1;
    typedef merge_tuple<tuple<Head, empty_tuple>, typename comp1::tail_t> comp2;
    typedef concat_tuple<typename comp1::used_t, typename comp2::type> concater;
    
public:
    typedef typename concater::type type;

    // This is the part of type into which we can still merge.
    typedef typename comp2::tail_t tail_t;

    typedef typename concat_tuple<
        typename comp1::used_t,
        typename comp2::used_t
      >::type used_t;

private:
    static typename comp1::type get_tail(const type& input)
    {
        typedef concat_tuple<
            typename comp1::used_t,
            typename comp1::tail_t
          > myconcat;
        return myconcat::doit(concater::get_first(input),
                    comp2::get_second(concater::get_second(input)));
    }

public:
    static tuple<Head, Tail> get_first(const type& input)
    {
        Head h = comp2::get_first(concater::get_second(input)).first();
        return tuple<Head, Tail>(h, comp1::get_first(get_tail(input)));
    }

    static Tuple2 get_second(const type& input)
    {
        return comp1::get_second(get_tail(input));
    }
};

// First argument a singleton, no merging
template<class T, class U, class Tail>
struct merge_tuple<tuple<T, empty_tuple>, tuple<U, Tail> >
{
private:
    typedef merge_tuple<tuple<T, empty_tuple>, Tail> comp;
    typedef concat_tuple<
        typename comp::used_t,
        tuple<U, typename comp::tail_t>
      > concater;

public:
    typedef typename comp::used_t used_t;
    typedef tuple<U, typename comp::tail_t> tail_t;
    typedef typename concater::type type;

private:
    static typename comp::type get_tail(const type& input)
    {
        typedef concat_tuple<
            typename comp::used_t,
            typename comp::tail_t
          > myconcat;
        return myconcat::doit(concater::get_first(input),
                concater::get_second(input).tail);
    }

public:
    static tuple<T, empty_tuple> get_first(const type& input)
    {
        return comp::get_first(get_tail(input));
    }

    static tuple<U, Tail> get_second(const type& input)
    {
        return tuple<U, Tail>(
                concater::get_second(input).head,
                comp::get_second(get_tail(input)));
    }
};

// First argument a singleton, with merging
template<class T, class Tail>
struct merge_tuple<tuple<T, empty_tuple>, tuple<T, Tail> >
{
    typedef tuple<T, Tail> type;
    typedef tuple<T, empty_tuple> used_t;
    typedef Tail tail_t;

    static tuple<T, empty_tuple> get_first(const type& input)
    {
        return make_tuple<T>::make(input.head);
    }

    static tuple<T, Tail> get_second(const type& input)
    {
        return input;
    }
};

// Termination case 1
template<class Tuple2>
struct merge_tuple<empty_tuple, Tuple2>
{
    typedef Tuple2 type;
    typedef type tail_t;
    typedef empty_tuple used_t;

    static empty_tuple get_first(const type& input)
    {
        return empty_tuple();
    }

    static Tuple2 get_second(const type& input)
    {
        return input;
    }
};

// It seems like this code path is unnecessary and in fact ambiguous.
// I am fairly convinced by now this is correct.
// However, in case issues ever come up, it seemed useful to me to retain this.
// -- Tom Bachmann (15/10/2013)
#if 0
// Termination case 2
template<class Tuple1>
struct merge_tuple<Tuple1, empty_tuple>
{
    typedef Tuple1 type;
    typedef type tail_t;
    typedef empty_tuple used_t;

    static Tuple1 get_first(const type& input)
    {
        return input;
    }

    static empty_tuple get_second(const type& input)
    {
        return empty_tuple();
    }
};
#endif

// Termination case 3
template<class T>
struct merge_tuple<tuple<T, empty_tuple>, empty_tuple>
{
    typedef tuple<T, empty_tuple> type;
    typedef empty_tuple tail_t;
    // NB: we "use" T here - it cannot be merged into it any longer!
    typedef type used_t;

    static empty_tuple get_second(const type& input)
    {
        return empty_tuple();
    }

    static tuple<T, empty_tuple> get_first(const type& input)
    {
        return input;
    }
};

// Termination case 4
template<>
struct merge_tuple<empty_tuple, empty_tuple>
{
    typedef empty_tuple type;
    typedef empty_tuple tail_t;
    typedef empty_tuple used_t;

    static empty_tuple get_first(const type& input) {return empty_tuple();}
    static empty_tuple get_second(const type& input) {return empty_tuple();}
};


// Creation and manipulation of homogeneous tuples

// Build a tuple type of "n" repetitions of "T".
template<class T, unsigned n>
struct make_homogeneous_tuple
{
    typedef tuple<T, typename make_homogeneous_tuple<T, n-1>::type> type;
};
template<class T>
struct make_homogeneous_tuple<T, 0>
{
    typedef empty_tuple type;
};

namespace htuples {
namespace hdetail {
template<unsigned n, class Tuple>
struct extract_helper
{
    typedef typename Tuple::head_t T;
    typedef typename Tuple::tail_t tail_t;
    typedef typename make_homogeneous_tuple<T, n>::type ht;

    static ht get_noskip(const Tuple& tuple)
    {
        return ht(tuple.head,
                extract_helper<n-1, tail_t>::get_noskip(
                    tuple.tail));
    }
};
template<class Tuple>
struct extract_helper<0, Tuple>
{
    static empty_tuple get_noskip(const Tuple&) {return empty_tuple();}
};

template<class Tuple>
struct remove_helper
{
    static const unsigned n = Tuple::len;
    typedef typename Tuple::tail_t tail_t;
    typedef typename Tuple::head_t T;
    static tail_t get(const Tuple& tuple, T res)
    {
        if(tuple.head == res)
            return tuple.tail;
        return tail_t(tuple.head, remove_helper<tail_t>::get(tuple.tail, res));
    }
};
template<class T>
struct remove_helper<tuple<T, empty_tuple> >
{
    static empty_tuple get(const tuple<T, empty_tuple>&, T)
    {
        return empty_tuple();
    }
};
} // hdetail

// Extract the first "n" values from the *homogeneous* tuple "tuple".
template<unsigned n, class Tuple>
inline typename make_homogeneous_tuple<typename Tuple::head_t, n>::type extract(
        const Tuple& tuple)
{
    return hdetail::extract_helper<n, Tuple>::get_noskip(tuple);
}

// Remove one element from the *homogeneous* tuple "tuple", if possible "res".
// Example:
//   t1 = (1, 2, 3, 4)
//   t2 = (1, 1, 1, 1)
//   t3 = (2, 3, 4, 5)
//
//   removeres(t1, 1) -> (2, 3, 4)
//   removeres(t2, 1) -> (1, 1, 1)
//   removeres(t3, 1) -> (2, 3, 4) or any other three element subset
template<class Tuple>
inline typename Tuple::tail_t
removeres(const Tuple& tuple, typename Tuple::head_t res)
{
    return hdetail::remove_helper<Tuple>::get(tuple, res);
}
} // htuples
} // mp

namespace traits {
// Compute if "Tuple" is of the form (U, U, .. U), or empty.
template<class Tuple, class U>
struct is_homogeneous_tuple : mp::and_<
    is_homogeneous_tuple<typename Tuple::tail_t, U>,
    mp::equal_types<typename Tuple::head_t, U> > { };
template<class U>
struct is_homogeneous_tuple<empty_tuple, U> : true_ { };

// Compute if T is a tuple
template<class T>
struct is_tuple : false_ { };
template<class Head, class Tail>
struct is_tuple<tuple<Head, Tail> > : true_ { };
template<>
struct is_tuple<empty_tuple> : true_ { };
}
} // flint

#endif
