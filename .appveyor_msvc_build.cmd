cd build.vc14\flint_config
python %cd%\flint_config.py --build-lib True --build-tests False --build-profiles False

cd ..
msbuild.exe lib_flint/lib_flint.vcxproj /p:Configuration=Release /p:Platform=%PLATFORM% /p:PlatformToolset=v140 /verbosity:minimal
copy lib_flint\%PLATFORM%\Release\lib_flint.lib ..\lib\%PLATFORM%\Release\lib_flint.lib

cd build_tests
python %cd%\build_tests.py --interfaces-tests False --platform %PLATFORM%

REM cd ..\run_tests
REM python %cd%\run_tests.py 0
REM cd ..\..
