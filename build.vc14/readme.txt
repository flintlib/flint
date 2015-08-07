Building FLINT2 with Microsoft Visual Studio 2015
-------------------------------------------------

Building FLINT2 with Microsoft Visual Studio requires Visual
Visual Studio 2015 Community (or higher version) and:

    a) an installed version of Python 3
    b) an installed version of Python Tools for 
       Visual Studio (http://pytools.codeplex.com/)

Obtain FLINT2 either as a released distribution or clone it using
GIT from:

    git@github.com:BrianGladman/flint2.git

FLINT2 depends on the MPIR, MPFR and PTHREADS libraries that have
to be installed and built using Vissual Studio before FLINT2 can
be built.  The application directories are assumed to be in the
same root directory with the names and layouts:
   
    mpir
       build.vc14
       lib
       dll
    mpfr  
       build.vc14
       lib
       dll
    pthreads  
       build.vc14
       lib
       dll
    flint2
       build.vc14
       lib
       dll
   
where the build.vc14 directories hold the Visual Studio build
files and the lib and dll directories hold the static and dynamic
library outputs for each package.  To libraries on which FLINT2
depends have to be built for the same configuration that will be 
used to build FLINT2 before FLINT2 itself can be built:
   
    <Static Library|Dynamic Link Library> 
    <Win32|x64>
    <Release|Debug>
   
where <a|b> shows the choices (a or b) that have to be made.   

Opening the solution file flint.sln in Visual Studio 2015 provides
the following build projects:

    dll_flint     - a Visual Studio build project for
                    FLINT2 as a Dynamic Link Library
    lib_flint     - a Visual Studio build project for
                    FLINT2 as a Static Library
    flint_config  - a Python program for creating the Visual 
                    Studio build files
    build_tests   - a Python program for building the FLINT2
                    tests (after they have been created)
     run_tests     - a Python program for running the FLINT2
                     tests (after they have been built)

The projects lib_flint and dll_flint can be used immediately to
build FLINT2 as a Static and Dynamic Link Library respectively.
Before building one or both of these, you need to select the
architecture (Win32 or x64) and the build type (Release or Debug).

To run the FLINT2 tests, the necessary Visual Studiop build files
have to be created.  If you have Python and Python Tools for
Visual Studio (PTVS) installed, this is done by setting the 
project flint_config (loaded into Visual Studio by the solution
file flint.sln) as the start-up project and then running it.
If you don't have PTVS installed but you do have Python, you
can run flint_config.py directly without Visual Studio. 
   
By default flint_config creates only the FLINT2 tests and profiling.
But it can also recreate the Visual Studio 2015 build files for the
FLINT2 DLL and Static Libraries by changing the defines at the
start of flint_config.py: 

    build_lib = False
    build_dll = False
    build_tests = Tru
     build_profiles = True

Rebuilding the library build files in this way may be necessary
if FLINT2 has been updated since it was first downloaded. 
      
After the FLINT2 tests have been created using flint_config.py,
they can then be built by setting build_tests.py as the start up
project and then running it.

There are also a number of Visual Studio solution files that
provide an *alternative* way of building the FLINT2 tests and
profiling.  However, their use is not recommended because each
of the multiple solution files flint-tests<NN>.sln (where NN
is a number) has to be loaded and built by Visual Studio (this
approach is used because it takes Visual Studio too long to
load the tests from a single solution file).

Once the tests have been built, the Python project run_tests can
be set as the start-up project and started to run all the tests 
(or the file run_tests.py can be run outside Visual Studio).
   
After building FLINT2, the libraries and the header files that 
you need to use FLINT2 are placed in the directories:
   
    lib\<Win32|x64>\<Debug|Release>
    dll\<Win32|x64>\<Debug|Release>

depending on the version(s) that have been built.
   
      Brian Gladman
      7th August 2015
