mkdir build
cd build
cmake -DCXXAPI=ON .. -DCMAKE_INSTALL_PREFIX=.
make
make install
export LD_LIBRARY_PATH=`pwd`/lib:$LD_LIBRARY_PATH
