.. _building:

**Configuring and building**
===============================================================================

Configuring Flint
-------------------------------------------------------------------------------

There are three ways to build Flint. One way is through a standard set of CMake
files provided in the distribution. A second way is with a custom build
system based on a configuration script and Makefiles also provided in the
distribution. A third way using MSVC solution files is documented in a section
further down.

The easiest way to use FLINT is to build a shared library. Simply download
the FLINT tarball and untar it on your system.

FLINT requires either MPIR (version 2.6.0 or later) or GMP (version 5.1.1 or
later). If MPIR is used, MPIR must be built with the ``--enable-gmpcompat``
option. FLINT also requires MPFR 3.0.0 or later and a pthread implementation.

Some of the input/output tests require ``fork`` and ``pipe``, however
these are disabled on Windows when a posix implementation is not provided.

If it happens that GMP/MPIR and MPFR are not in a standard location on your
system (e.g. not in ``/usr/include/`` and ``/usr/lib/``), you need to tell the
configure script where they are with the options ``--with-gmp=/path/to/gmp``
or ``--with-mpir=/path/to/mpir`` and ``--with-mpfr=/path/to/mpfr``, e.g.

.. code-block:: bash

   ./configure --with-gmp=/home/user1/local --with-mpfr=/home/user1/local

FLINT can also build against a source build of GMP/MPIR and MPFR. Though
programs using FLINT may require GMP/MPIR and MPFR to be installed (via
``make install`` if built from sources).

Note that FLINT builds static and shared libraries by default, except on
platforms where this is not supported. If you do not require either a shared
or static library then you may pass ``--disable-static`` or
``--disable-shared`` to ``configure``. This can substantially speed up the
build.

If you intend to install the FLINT library and header files, you can specify
where they should be placed by passing ``--prefix=path`` to configure, where
``path`` is the directory under which the ``lib`` and ``include`` directories
exist into which you wish to place the FLINT files when it is installed.

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

C++ wrapper
-------------------------------------------------------------------------------

If you wish to enable the test functions for the FLINT ``C++`` wrapper
``flintxx`` you must pass ``--enable-cxx`` to configure.

The ``C++`` wrapper is always available, but tests will only run if
this option is selected. It is disabled by default (``--disable-cxx``)
because some ``C++`` compilers internally segfault when compiling the
tests, or exhaust memory due to the complexity of the ``C++`` code.

Building, testing, installing and using FLINT
-------------------------------------------------------------------------------

Once FLINT is configured, in the main directory of the FLINT directory
tree simply type:

.. code-block:: bash

    make
    make check

GNU make or CMake is required to build FLINT. For GNU make, this is simply
``make`` on Linux, Darwin, MinGW and Cygwin systems. However, on some unixes
the command is ``gmake``. For CMake, this is ``cmake``, which which may need
to be installed with your package manager.

If you wish to install FLINT with GNU make, simply type:

.. code-block:: bash

    make install

If you wish to install FLINT with CMake, simply type:

.. code-block:: bash

    mkdir build && cd build
    cmake .. -DBUILD_SHARED_LIBS=ON
    cmake --build . --target install

Now to use FLINT, simply include the appropriate header files for the FLINT
modules you wish to use in your C program.  Then compile your program,
linking against the FLINT library, GMP/MPIR, MPFR and pthreads with the
options ``-lflint -lmpfr -lgmp -lpthread``.

To uninstall FLINT with GNU make, type:

.. code-block:: bash

    make uninstall

Note that you may have to set ``LD_LIBRARY_PATH`` or equivalent for your
system to let the linker know where to find these libraries. Please refer to
your system documentation for how to do this.

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

The FLINT custom make system responds to the standard commands

.. code-block:: bash

    make
    make library
    make check
    make clean
    make distclean
    make install

If your system supports parallel builds, FLINT will build in parallel, e.g:

.. code-block:: bash

    make -j4 check

On some systems, parallel builds appear to be available but buggy.

Testing a single module or file
-------------------------------------------------------------------------------

If you wish to simply check a single module of FLINT you can pass the option
``MOD=modname`` to ``make check``. You can also pass a list of module names in
inverted commas, e.g:

.. code-block:: bash

    make check MOD=ulong_extras
    make check MOD="fft fmpz_mat"

To specify an individual test(s) for any module you can add it (or comma
separated test list) after chosen module name followed by the colon, e.g.:

.. code-block:: bash

    make check MOD=ulong_extras:clog,factor,is_prime
    make check MOD="fft fmpz_mat:add_sub,charpoly fq_vec:add"

Assertion checking
-------------------------------------------------------------------------------

FLINT has an assert system. If you want a debug build you can pass
``--enable-assert`` to configure. However, this will slow FLINT considerably,
so asserts should not be enabled (``--disable-assert``, the default) for
deployment.

Exceptions
-------------------------------------------------------------------------------

When FLINT encounters a problem, mostly illegal input, it currently aborts.
There is an experimental interface for generating proper exceptions
``flint_throw``, but this is currently rarely used and experimental - you
should expect this to change.

At the end, all of FLINT's exceptions call ``abort()`` to terminate
the program. Using ``flint_set_abort(void (*abort_func)(void))``, the
user can install a function that will be called instead. Similar
to the exceptions, this should be regarded as experimental.

Building FLINT2 with Microsoft Visual Studio using solution files
-------------------------------------------------------------------------------

Brian Gladman has kindly provided the build scripts for building
Flint with Microsoft Visual Studio.

Building FLINT2 with Microsoft Visual Studio requires 
Visual Studio 2015 Community (or higher version) and:

- an installed version of Python 3

- an installed version of Python Tools for Visual Studio
  <https://github.com/Microsoft/PTVS>

Obtain FLINT2 by cloning it using GIT from Brian Gladman's repository:

  ``git@github.com:BrianGladman/flint.git``

FLINT2 depends on the MPIR, MPFR and PTHREADS libraries that have
to be installed and built using Visual Studio before FLINT2 can
be built.  The application directories are assumed to be in the
same root directory with the names and layouts:

.. code ::

    mpir
       lib
       dll
    mpfr  
       lib
       dll
    pthreads  
       lib
       dll
    flint
       build.vc
       lib
       dll
   
Here the ``lib`` and ``dll`` sub-directories for each application hold
the  static and dynamic link library outputs which will be used when 
Flint is built.  They each contain up to four sub-directories for
the normal configurations for building on Windows:

    ``Win32\Release``

    ``Win32\Debug``

    ``x64\Release``

    ``x642\Debug``
    
To build FLINT2 for a particular configuration requires that each of the 
three libraries on which FLINT2 depends must have been previously built
for the same configuration.

Opening the solution file ``flint\build.vc\flint.sln`` provides the
following build projects:

    ``flint_config``  - a Python program for creating the Visual Studio build files

    ``build_tests``   - a Python program for building the FLINT2 tests (after they have been created)

    ``run_tests``     - a Python program for running the FLINT2 tests (after they have been built)

The first step in building FLINT2 is to generate the Visual Studio build
files for the version of Visual Studio being used. This is done by
running the Python application ``flint_config.py``, either from within
Visual Studio or on the command line. It is run with a single input
parameter which is the last two digits of the Visual Studio version
selected for building FLINT2 (the default is 19 if no input is given).

This creates a build directory in the Flint root directory, for 
example:

   ``flint\build.vs19``
   
that contains the file ``flint.sln`` which can now be loaded into
Visual Studio and used to build the FLINT2 library.

Once the FLINT2 library has been built, the FLINT2 tests can now be 
built and run by returning to the Visual Studio solution:

  ``flint\build.vc\flint.sln``
  
and running the ``build_tests`` and ``run_tests`` Python applications.
  
After building FLINT2, the libraries and the header files that 
you need to use FLINT2 are placed in the directories:

- ``lib\<Win32|x64>\<Debug|Release>``

- ``dll\<Win32|x64>\<Debug|Release>``

depending on the version(s) that have been built.


