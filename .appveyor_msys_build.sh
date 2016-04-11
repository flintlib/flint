export PATH=/c/msys64/mingw$ABI/bin:$PATH

cd /c/projects/flint2

./.build_dependencies
./configure ABI=$ABI --with-mpir=mpir-2.7.0 --with-mpfr=mpfr-3.1.4 --disable-shared
make -j4
make -j4 check
