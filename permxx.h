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

#ifndef PERMXX_H
#define PERMXX_H

#include "perm.h"

#include "flintxx/frandxx.h"
#include "flintxx/mp.h"

namespace flint {
class permxx
{
private:
    slong* data;
    slong sz;

public:
    permxx(slong n) {data = _perm_init(n);sz = n;}
    ~permxx() {_perm_clear(data);}

    permxx(const permxx& o)
    {
        sz = o.size();
        data = _perm_init(sz);
        _perm_set(data, o._data(), sz);
    }
    permxx& operator=(const permxx& o)
    {
        sz = o.size();
        _perm_set(_data(), o._data(), sz);
        return *this;
    }

    bool operator==(const permxx& o)
        {return size() == o.size() && _perm_equal(_data(), o._data(), size());}
    bool operator!=(const permxx& o) {return !(*this == o);}

    static permxx one(slong n) {return permxx(n);}
    static permxx randtest(slong n, frandxx& state)
        {permxx res(n);res.set_randtest(state);return res;}

    void set_one() {_perm_set_one(_data(), size());}
    int set_randtest(frandxx& state)
        {return _perm_randtest(_data(), size(), state._data());}

    slong* _data() {return data;}
    const slong* _data() const {return data;}
    slong size() const {return sz;}

    slong& operator[](slong idx) {return data[idx];}
    slong operator[](slong idx) const {return data[idx];}

    int parity() const {return _perm_parity(_data(), size());}

    permxx operator*(const permxx& o) const
    {
        permxx res(o.size());
        _perm_compose(res._data(), _data(), o._data(), size());
        return res;
    }
    permxx& operator*=(const permxx& o)
    {
        _perm_compose(_data(), _data(), o._data(), size());
        return *this;
    }

    void set_inv(const permxx& o) {_perm_inv(_data(), o._data(), o.size());}
    permxx inv() const {permxx res(size());res.set_inv(*this);return res;}
};

inline permxx compose(const permxx& p1, const permxx& p2) {return p1*p2;}
inline int parity(const permxx& p) {return p.parity();}
inline permxx inv(const permxx& o) {return o.inv();}

inline slong* maybe_perm_data(permxx* p) {return p ? p->_data() : 0;}
inline slong* maybe_perm_data(int zero) {return 0;}

namespace traits {
template<class T> struct is_maybe_perm
    : mp::or_<mp::equal_types<T, permxx*>, mp::equal_types<T, int> > { };
template<class T> struct is_permxx : mp::equal_types<T, permxx> { };
} // traits

inline int print(const permxx& p)
{
    return _perm_print(p._data(), p.size());
}
} // flint

#endif
