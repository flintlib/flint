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

#ifndef CXX_FORWARDING_H
#define CXX_FORWARDING_H

#include <iosfwd>
#include <string>

#include "evaluation_tools.h"
#include "mp.h"
#include "rules.h"
#include "tuple.h"

// XXX This file is UNFINISHED. It is/was meant to use for extending flintxx by
// compisition, in particular for the NTL wrapper.
// Ask Tom (e_mc_h2@web.de) if you want to use it.

namespace flint {
namespace forwarding {
// Specialise this to enable forwarding of your class.
template<class T>
struct enable : mp::false_ { };
} // forwarding

namespace rules {

template<class T>
struct print<T, typename mp::enable_if<traits::is_implemented<
        print<typename forwarding::enable<T>::underlying_t> > >::type>
{
    static void doit(const T& t, std::ostream& o)
    {
        typedef typename forwarding::enable<T> et;
        print<typename et::underlying_t>::doit(et::get_underlying(t), o);
    }
};

template<class T>
struct to_string<T, typename mp::enable_if<traits::is_implemented<
        to_string<typename forwarding::enable<T>::underlying_t> > >::type>
{
    static std::string get(const T& t, int base)
    {
        typedef typename forwarding::enable<T> et;
        return to_string<typename et::underlying_t>::get(
                et::get_underlying(t), base);
    }
};

namespace rdetail {
template<class T, class Enable = void>
struct maybe_forward
{
    typedef T underlying_t;
    static const T& get_underlying(const T& t) {return t;}
    static T& get_underlying(T& t) {return t;}
};
template<class T>
struct maybe_forward<T, typename mp::enable_if<forwarding::enable<T> >::type>
    : forwarding::enable<T> { };

template<template<class, class, class> class Wrapped, class T, class U>
struct fwd_enable2
    : mp::enable_if<mp::and_<
          traits::is_implemented<Wrapped<
              typename maybe_forward<T>::underlying_t,
              typename maybe_forward<U>::underlying_t,
              void> >,
          mp::or_<
              typename forwarding::enable<T>,
              typename forwarding::enable<U> > > >
{
    typedef Wrapped<
              typename maybe_forward<T>::underlying_t,
              typename maybe_forward<U>::underlying_t,
              void> wrapped_t;
};
} // rdetail

template<class T, class U>
struct assignment<T, U, typename rdetail::fwd_enable2<rules::assignment, T, U>::type>
{
    static void doit(T& t, const U& u)
    {
        rdetail::fwd_enable2<rules::assignment, T, U>::wrapped_t::doit(
                rdetail::maybe_forward<T>::get_underlying(t),
                rdetail::maybe_forward<U>::get_underlying(u));
    }
};

template<class T, class U>
struct cmp<T, U, typename rdetail::fwd_enable2<rules::cmp, T, U>::type>
{
    static int get(const T& t, const U& u)
    {
        return rdetail::fwd_enable2<rules::cmp, T, U>::wrapped_t::get(
                rdetail::maybe_forward<T>::get_underlying(t),
                rdetail::maybe_forward<U>::get_underlying(u));
    }
};

template<class T, class U>
struct equals<T, U, typename rdetail::fwd_enable2<rules::equals, T, U>::type>
{
    static bool get(const T& t, const U& u)
    {
        return rdetail::fwd_enable2<rules::equals, T, U>::wrapped_t::get(
                rdetail::maybe_forward<T>::get_underlying(t),
                rdetail::maybe_forward<U>::get_underlying(u));
    }
};

template<class T, class U>
struct conversion<T, U, typename mp::enable_if<traits::is_implemented<
    conversion<T, typename forwarding::enable<U>::underlying_t> > >::type>
{
    static T get(const U& u)
    {
        return conversion<T, typename forwarding::enable<U>::underlying_t>::get(
                forwarding::enable<U>::get_underlying(u));
    }
};

// TODO automatically convert fwd<U> to fwd<T> ?

// TODO instantiate_temporaries ?


namespace rdetail {
template<class Data>
struct needs_forwarding : mp::false_ { };

template<class Data>
struct needs_forwarding<const Data&> : forwarding::enable<Data> { };

template<class Head, class Tail>
struct needs_forwarding<tuple<Head, Tail> >
    : mp::or_<needs_forwarding<Head>, needs_forwarding<Tail> > { };

template<class Data>
struct translate
{
    typedef Data type;
    static type doit(const Data& d) {return d;}
};

template<class Wrapped>
struct translate<const Wrapped&>
{
    typedef maybe_forward<Wrapped> fw;
    typedef typename fw::underlying_t type;
    static type doit(const Wrapped& d)
    {
        return fw::get_underlying(d);
    }
};

template<class Head, class Tail>
struct translate<tuple<Head, Tail> >
{
    typedef translate<Head> htranslator;
    typedef translate<Tail> ttranslator;
    typedef tuple<typename htranslator::type, typename ttranslator::type> type;
    static type doit(const tuple<Head, Tail>& d)
    {
        return type(htranslator::doit(d.head), ttranslator::doit(d.tail));
    }
};

template<class Op, class Data, bool result_is_temporary, class Enable = void>
struct can_forward
{
    typedef rdetail::translate<Data> translator;
    typedef typename translator::type translated_t;
    typedef typename mp::find_evaluation<
        Op, translated_t, result_is_temporary>::type rule_t;
    static const bool val = traits::is_implemented<rule_t>::val;
};

template<class Op, class Data, bool result_is_temporary>
struct can_forward<Op, Data, result_is_temporary,
    typename mp::disable_if<needs_forwarding<Data> >::type>
{
    static const bool val = false;
};
}

template<class Op, class Data, bool result_is_temporary>
struct evaluation<Op, Data, result_is_temporary, 1,
    typename mp::enable_if<
        rdetail::can_forward<Op, Data, result_is_temporary> >::type>
{
    typedef rdetail::translate<Data> translator;
    typedef typename translator::type translated_t;
    typedef typename mp::find_evaluation<
        Op, translated_t, result_is_temporary>::type rule_t;

    typedef typename rule_t::temporaries_t temporaries_t;
    // TODO repackaging
    typedef typename rule_t::return_t return_t;

    static void doit(const Data& input, temporaries_t temps, return_t* output)
    {
        rule_t::doit(translator::doit(input), temps, output);
    }
};

} //rules
} // flint

#endif
