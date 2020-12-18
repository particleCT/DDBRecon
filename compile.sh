export G4WORKDIR=$PWD
mkdir build
cd build

rm -rf *

cmake -DCMAKE_INSTALL_PREFIX=$G4WORKDIR ../

make -j8

make install

cd ../

