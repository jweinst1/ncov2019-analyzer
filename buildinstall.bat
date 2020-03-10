echo "Building the cmake targets"
mkdir build_c_bin
cd build_c_bin
cmake -DWITH_testing=ON ..
msbuild ALL_BUILD.vcxproj
msbuild INSTALL.vcxproj
cd ..
rmdir /s /q build_c_bin
echo "Finished build and run"