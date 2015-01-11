compiler=clang++
#compiler=icpc
#compiler=g++
type=Debug
#type=Release

mkdir build
rm -rf ./build/*
cd build
cmake -DCMAKE_BUILD_TYPE=$type -DCMAKE_CXX_COMPILER=$compiler ../src
make -j4
make install
