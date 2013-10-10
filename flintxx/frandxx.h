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

#ifndef CXX_FRANDXX_H
#define CXX_FRANDXX_H

#include "../flint.h"

// This class contains a first-class wrapper of flint_rand_t.
// Note that frandxx is not copyable.

namespace flint {
class frandxx
{
private:
    flint_rand_t inner;

    // not copyable
    frandxx(const frandxx&);

public:
    frandxx() {flint_randinit(inner);}
    ~frandxx() {flint_randclear(inner);}

    flint_rand_t& _data() {return inner;}
    const flint_rand_t& _data() const {return inner;}
};
} // flint

#endif
