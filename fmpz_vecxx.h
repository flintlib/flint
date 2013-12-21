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

#ifndef CXX_FMPZ_VECXX_H
#define CXX_FMPZ_VECXX_H

#include "fmpzxx.h"
#include "fmpz_vec.h"
#include "flintxx/vector.h"

namespace flint {
namespace detail {
struct fmpz_vector_data
{
    slong size;
    fmpz* array;

    fmpz_vector_data(slong n)
        : size(n), array(_fmpz_vec_init(n)) {}

    ~fmpz_vector_data() {_fmpz_vec_clear(array, size);}

    fmpz_vector_data(const fmpz_vector_data& o)
        : size(o.size), array(_fmpz_vec_init(o.size))
    {
        _fmpz_vec_set(array, o.array, size);
    }

    fmpzxx_ref at(slong i) {return fmpzxx_ref::make(array + i);}
    fmpzxx_srcref at(slong i) const {return fmpzxx_srcref::make(array + i);}
};
} // detail

typedef vector_expression<
    detail::wrapped_vector_traits<fmpzxx, slong, fmpzxx_ref, fmpzxx_srcref, fmpz>,
    operations::immediate,
    detail::fmpz_vector_data> fmpz_vecxx;
// TODO references

template<>
struct enable_vector_rules<fmpz_vecxx> : mp::false_ { };

namespace rules {
// TODO hack to make code look like references are implemented
template<class T> struct FMPZ_VECXX_COND_S : mp::equal_types<T, fmpz_vecxx> { };
#define FMPZ_VECXX_COND_T FMPZ_VECXX_COND_S

template<class T>
struct to_string<T, typename mp::enable_if<FMPZ_VECXX_COND_S<T> >::type>
{
    static std::string get(const fmpz_vecxx& e, int base)
    {
        // TODO use _fmpz_vec_print somehow?
        std::ostringstream o;
        o << e.size();
        if(e.size() == 0)
        {
            return o.str();
        }
        o << "  ";
        for(slong i = 0;i < e.size();++i)
        {
            o << e[i].to_string(base);
            if(i != e.size() - 1)
                o << " ";
        }
        return o.str();
    }
};

// TODO references
FLINT_DEFINE_GET(equals, bool, fmpz_vecxx,
        e1.size() == e2.size()
        && _fmpz_vec_equal(e1._data().array, e2._data().array, e1.size()))

FLINT_DEFINE_BINARY_EXPR_COND2(plus, fmpz_vecxx,
        FMPZ_VECXX_COND_S, FMPZ_VECXX_COND_S,
        _fmpz_vec_add(to._data().array, e1._data().array, e2._data().array,
            e1.size()))

FLINT_DEFINE_DOIT_COND2(assignment, FMPZ_VECXX_COND_T, FMPZ_VECXX_COND_S,
        _fmpz_vec_set(to._array(), from._array(), to.size()))

// TODO more
} // rules
} // flint

#endif
