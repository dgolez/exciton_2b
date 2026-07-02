module purge
module load  foss/2022a HDF5/1.12.2-gompi-2022a Eigen/3.4.0-GCCcore-12.2.0 CMake/3.24.3-GCCcore-12.2.0

mkdir cbuild
cd cbuild

CC=mpicc CXX=mpicxx
cmake \
 -DCMAKE_INSTALL_PREFIX="$HOME/opt" \
 -DCMAKE_INSTALL_NAME_DIR="$HOME/opt/lib" \
 -DCMAKE_BUILD_TYPE=Release \
 -Domp=ON \
 -Dmpi=ON \
 -Dhdf5=ON \
 -DCMAKE_INCLUDE_PATH="$HOME/opt/include"  \
 -DCMAKE_LIBRARY_PATH="$HOME/opt/lib" \
 -DCMAKE_CXX_FLAGS="-std=c++11 -O3 -lfftw3_omp -lfftw3 -lm" \
 ..

make -j

make install
