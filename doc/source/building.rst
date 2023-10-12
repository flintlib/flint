.. _building:

**Building, testing and installing**
===============================================================================

Quick start
-------------------------------------------------------------------------------

Building FLINT requires:

* GMP (https://gmplib.org/)
* MPFR (https://mpfr.org/)
* Either of the following build systems:

  * GNU Make together with GNU Autotools
  * CMake

On a typical Linux or Unix-like system where Autotools is available (see below
for instructions using CMake), FLINT can be built and installed as follows:

.. code-block:: bash

    ./bootstrap.sh
    ./configure --disable-static
    make -j
    make install

We also recommend that you run ``make check`` before installing.

For a complete list of build settings, write

.. code-block:: bash

    ./configure --help

An example of a custom configuration command would be

.. code-block:: bash

    ./configure                                         \
        --enable-assert                                 \
        --enable-avx2                                   \
        --disable-static                                \
        --with-gmp-include=/home/user1/builds/includes/ \
        --with-gmp-lib=/home/user1/builds/lib/          \
        --with-mpfr=/usr                                \
        --prefix=/home/user1/installations/             \
        CC=clang                                        \
        CFLAGS="-Wall -O3 -march=alderlake"

Library and install paths
-------------------------------------------------------------------------------

If you intend to install the FLINT library and header files, you can specify
where they should be placed by passing ``--prefix=path`` to ``configure``, where
``path`` is the directory under which the ``lib`` and ``include`` directories
exist into which you wish to place the FLINT files when it is installed.

If GMP and MPFR are not installed in the default search path of your compiler
(e.g. ``/usr/include/`` and ``/usr/lib/``), you must specify where they are by
passing their location to configure ``--with-gmp=ABSOLUTE_PATH`` for GMP and
``--with-mpfr=ABSOLUTE_PATH`` for MPFR.
Note that the FLINT build system can handle GMP and MPFR as installed at some
location and as source builds (built from source but not installed).  Though, to
run the FLINT tests, GMP and MPFR needs to be properly installed.

Testing FLINT
-------------------------------------------------------------------------------

The full FLINT test suite can be run using

.. code-block:: bash

    make check

or in parallel on a multicore system using

.. code-block:: bash

    make -j check

Number of test iterations
...............................................................................

The number of test iterations can be changed with the
``FLINT_TEST_MULTIPLIER`` environment variable. For example, the
following will only run 10% of the default iterations::

    export FLINT_TEST_MULTIPLIER=0.1
    make check

Conversely, ``FLINT_TEST_MULTIPLIER=10`` will stress test FLINT
by performing 10x the default number of iterations.

Testing single modules
...............................................................................

If you wish to simply check a single module of FLINT you can pass the option
``MOD=modname`` to ``make check``. You can also pass a list of module names:

.. code-block:: bash

    make check MOD=ulong_extras
    make check MOD="fft fmpz_mat"

Test coverage
...............................................................................

To obtain coverage statistics for the FLINT test suite, assuming
that ``gcov`` and ``lcov`` are installed, configure
FLINT with ``--enable-coverage``. Then run:

.. code-block:: bash

    make -j check
    make coverage

This will place a coverage report in ``build/coverage``.


Static or dynamic library only
-------------------------------------------------------------------------------

FLINT builds static and shared libraries by default, except on
platforms where this is not supported. If you do not require either a shared
or static library then you may pass ``--disable-static`` or
``--disable-shared`` to ``configure``. This can substantially speed up the
build.

AVX2 instructions
-------------------------------------------------------------------------------

On x86-64 machines with AVX2 support, compiling FLINT with the ``--enable-avx2``
option can improve performance substantially, notably by enabling
the small-prime FFT. Currently this option is not enabled by default.

TLS, reentrancy and single mode
-------------------------------------------------------------------------------

FLINT uses thread local storage by default (``--enable-tls``). However, if
reentrancy is required on systems that do not support this, one can pass
``--disable-tls`` and mutexes will be used instead (requires POSIX). As most
modern systems support thread local storage, it is not recommended to build
FLINT without TLS.

There are two modes in which FLINT may installed: the default ``single`` mode,
which is faster, but makes use of thread local storage for its memory manager
and to handle threading, and a slower but less complicated ``reentrant`` mode.
The later is useful when debugging a program where tracing allocations is
important.

If you wish to select the single mode, pass the ``--single`` option to
configure, though note that this is the default. The reentrant mode is selected
by passing the option ``--reentrant`` to configure.

ABI and architecture support
-------------------------------------------------------------------------------

On some systems, e.g. Sparc and some Macs, more than one ABI is available.
FLINT chooses the ABI based on the CPU type available, however its default
choice can be overridden by passing either ``ABI=64`` or ``ABI=32`` to
configure.

To build on MinGW64 it is necessary to pass ``ABI=64`` to configure, as FLINT
is otherwise unable to distinguish it from MinGW32.

In some cases, it is necessary to override the CPU/OS defaults. This can be
done by passing ``--build=cpu-os`` to configure.

The available choices for CPU include ``x86_64``, ``x86``, ``ia64``, ``sparc``,
``sparc64``, ``ppc``, ``ppc64``. Other CPU types are unrecognised and FLINT
will build with generic code on those machines.

The choices for OS include ``Linux``, ``MINGW32``, ``MINGW64``, ``CYGWIN32``,
``CYGWIN64``, ``Darwin``, ``FreeBSD``, ``SunOS`` and numerous other operating
systems.

It is also possible to override the default CC, AR and CFLAGS used by FLINT by
passing ``CC=full_path_to_compiler``, etc., to FLINT's configure.


CMake build
-------------------------------------------------------------------------------

If you wish to install FLINT with CMake, simply type:

.. code-block:: bash

    mkdir build && cd build
    cmake .. -DBUILD_SHARED_LIBS=ON
    cmake --build . --target install

Uninstalling FLINT
-------------------------------------------------------------------------------

To uninstall FLINT with GNU make, type:

.. code-block:: bash

    make uninstall

Now to use FLINT, simply include the appropriate header files for the FLINT
modules you wish to use in your C program.  Then compile your program,
linking against the FLINT library, GMP/MPIR, MPFR and pthreads with the
options ``-lflint -lmpfr -lgmp -lpthread``.

To clean up the local build files, use:

.. code-block:: bash

    make clean
    make distclean

Assertion checking
-------------------------------------------------------------------------------

FLINT has an assert system. If you want a debug build you can pass
``--enable-assert`` to configure. However, this will slow FLINT considerably,
so asserts should not be enabled (``--disable-assert``, the default) for
deployment.

Linking and running code
-------------------------------------------------------------------------------

Here is an example program to get started using FLINT:

.. code-block:: c

    #include "flint/flint.h"
    #include "flint/arb.h"

    int main()
    {
        arb_t x;
        arb_init(x);
        arb_const_pi(x, 50 * 3.33);
        arb_printn(x, 50, 0); flint_printf("\n");
        flint_printf("Computed with FLINT-%s\n", flint_version);
        arb_clear(x);
    }

Compile it with::

    gcc test.c -lflint

You may also have to pass the flags ``-lmpfr`` and ``-lgmp`` to the compiler.
If the FLINT header and library files are not in a standard location
such as ``/usr/local``, you may also have to provide flags such as::

    -I/path/to/flint -L/path/to/flint

Finally, to run the program, make sure that the linker
can find ``libflint``. If it is installed in a
nonstandard location, you can for example add this path to the
``LD_LIBRARY_PATH`` environment variable.

The output of the example program should be something like the following::

    [3.1415926535897932384626433832795028841971693993751 +/- 4.43e-50]
    Computed with flint-3.0.0

Header file conflicts
-------------------------------------------------------------------------------

If you have any difficulties with conflicts with system headers on your
machine, you can do the following in your code:

.. code-block:: C

    #undef ulong
    #define ulong ulongxx
    #include <stdio.h>
    // other system headers
    #undef ulong
    #define ulong mp_limb_t

This prevents FLINT's definition of ``ulong`` interfering with your system
headers.

