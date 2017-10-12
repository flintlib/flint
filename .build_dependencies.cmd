echo int flint_test_multiplier(){return 1;} > test_helpers.c

cd ..

curl -fsS -o mpir.zip http://mpir.org/mpir-2.7.2.zip
7z x mpir.zip > NUL
rename mpir-2.7.2 mpir
cd mpir\build.vc14
msbuild.exe lib_mpir_gc/lib_mpir_gc.vcxproj /p:Configuration=Release /p:Platform=%PLATFORM% /p:PlatformToolset=v140 /verbosity:minimal
cd ..\..

curl -fsSL -o mpfr.zip https://github.com/BrianGladman/mpfr/archive/master.zip
7z x mpfr.zip > NUL
rename mpfr-master mpfr
cd mpfr\build.vc14
msbuild.exe lib_mpfr.sln /p:Configuration=Release /p:Platform=%PLATFORM% /p:PlatformToolset=v140 /verbosity:minimal
cd ..\..

curl -fsSL -o pthreads.zip https://github.com/BrianGladman/pthreads/archive/master.zip
7z x pthreads.zip > NUL
rename pthreads-master pthreads
cd pthreads\build.vc14
msbuild.exe pthreads.sln /p:Configuration=Release /p:Platform=%PLATFORM% /p:PlatformToolset=v140 /verbosity:minimal
cd ..\..
 
cd flint2
