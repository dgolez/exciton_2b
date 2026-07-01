module purge
module load  foss/2022a HDF5/1.12.2-gompi-2022a Eigen/3.4.0-GCCcore-12.2.0 CMake/3.24.3-GCCcore-12.2.0
#Boost/1.81.0-GCC-12.2.0


#module purge
#module load GCCcore/11.2.0 CMake/3.22.1-GCCcore-11.2.0 HDF5 Python/3.9.6-GCCcore-11.2.0
#module load HDF5/1.12.0-gompi-2021a
#module load CMake
#module load Eigen
#module load Boost

# module load GCC
# module load OpenMPI
# module load HDF5/1.10.7-iimpi-2020b
# module load Eigen
# module load Boost
# module load CMake

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
 -DCMAKE_INCLUDE_PATH="/home/golez/opt/include/"  \
 -DCMAKE_LIBRARY_PATH="/home/golez/opt/lib" \
 -DCMAKE_CXX_FLAGS="-std=c++11 -O3" \
 .. && make -j

make install
