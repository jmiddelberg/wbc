#!/bin/sh
sudo apt-get -y install git cmake build-essential libboost-system-dev libboost-program-options-dev libboost-thread-dev libboost-test-dev pkg-config libeigen3-dev libboost-filesystem-dev

# cmake makros
git clone https://github.com/rock-core/base-cmake.git   
mkdir base-cmake/build && cd base-cmake/build
cmake .. && make -j8 && sudo make install && cd ../..

# Logging 
git clone https://github.com/rock-core/base-logging.git
mkdir base-logging/build && cd base-logging/build
cmake .. && make -j8 && sudo make install && cd ../..

# Base Types
git clone https://github.com/rock-core/base-types.git
mkdir base-types/build && cd base-types/build
cmake .. -DUSE_SISL=OFF -DBINDINGS_RUBY=OFF -DROCK_VIZ_ENABLED=OFF
make -j8 && sudo make install && cd ../..

# URDF
sudo apt-get -y install liburdfdom-headers-dev liburdfdom-dev 
# compatability with urdfdom >= 4.0.0 (Ubuntu24.04)´
sudo apt-get -y install libtinyxml2-dev

# Clone WBC repo to have the patches 
git clone https://github.com/ARC-OPT/wbc.git

# RBDL
git clone --branch v3.2.1 --recurse-submodules https://github.com/rbdl/rbdl.git
cd rbdl
git apply ../wbc/patches/rbdl.patch --ignore-whitespace
mkdir build && cd build
cmake .. -DRBDL_BUILD_ADDON_URDFREADER=ON
make -j8 && sudo make install && cd ../..

# Pinocchio
git clone --branch v2.6.8 --recurse-submodules https://github.com/stack-of-tasks/pinocchio.git
cd pinocchio
mkdir build && cd build
cmake .. -DBUILD_PYTHON_INTERFACE=OFF -DBUILD_UNIT_TESTS=OFF 
make -j8 && sudo make install && cd ../..

# If not done yet, setup a ssh key pair using the command `ssh-keygen` and add the 
# key from `~/.ssh/id_rsa.pub `to the keys in your Gitlab account.

# qpOASES
git clone https://github.com/coin-or/qpOASES.git -b releases/3.2.0
cd qpOASES
mkdir patches && cp ../wbc/patches/qpOASES.patch patches
git apply patches/qpOASES.patch
mkdir build && cd build
cmake .. && make -j8 && sudo make install && cd ../..

# eiquadprog
git clone --recurse-submodules https://github.com/stack-of-tasks/eiquadprog.git -b v1.2.5
cd eiquadprog
cp ../wbc/patches/eiquadprog.patch . && git apply eiquadprog.patch
mkdir build && cd build
cmake ..
make -j8 && sudo make install && cd ../..

# qpSWIFT
git clone https://github.com/qpSWIFT/qpSWIFT.git
cd qpSWIFT
cp ../wbc/patches/qpSWIFT.patch . && git apply qpSWIFT.patch
mkdir build && cd build
cmake .. 
make -j8 && sudo make install && cd ../.. 

# proxQP
# check out the same proxsuite version as acados, otherwise the two copies collide:
# proxsuite is header-only, so both wbc-solvers-proxqp and acados' libproxsuite_c
# end up carrying their own compiled copy of the same C++ symbols, and when both are
# loaded into one process the dynamic linker picks a single winner for all callers.
# acados pins FranekStark's fork of 0.7.3, which carries a not-yet-merged fix for a
# preconditioner bug that corrupts the dense backend's box-constraint multipliers
# when a QP object is re-solved.
git clone --branch fix/dense-box-i-scaled-not-reset-on-keep --recurse-submodules https://github.com/FranekStark/proxsuite.git proxqp
cd proxqp
git checkout 0bd12daa6c54f22361744fcaa2ca5228e7f719d8
git submodule update --recursive --init
mkdir build && cd build
cmake .. -DBUILD_TESTING=OFF -DBUILD_PYTHON_INTERFACE=OFF -DBUILD_WITH_VECTORIZATION_SUPPORT=OFF
make -j8 && sudo make install && cd ../..

# OSQP
git clone https://github.com/osqp/osqp.git
cd osqp
mkdir build && cd build
cmake ..
make -j8 && sudo make install && cd ../..
git clone https://github.com/robotology/osqp-eigen.git
cd osqp-eigen
mkdir build && cd build
cmake ..
make -j8 && sudo make install && cd ../..

# Clarabel
# Rust toolchain (required to build Clarabel). Installs the latest stable rustc/cargo.
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
. "$HOME/.cargo/env"
# Clarabel.cpp provides no `make install`, so build the Rust C-library (the C++ headers are
# committed in include/) and copy the headers + shared library into /usr/local manually.
git clone --branch v0.11.1 --recurse-submodules https://github.com/oxfordcontrol/Clarabel.cpp.git
cd Clarabel.cpp/rust_wrapper
cargo build --release
cd ../..
sudo cp -r Clarabel.cpp/include /usr/local/include/clarabel
sudo cp Clarabel.cpp/rust_wrapper/target/release/libclarabel_c.so /usr/local/lib/

# WBC
mkdir wbc/build && cd wbc/build
cmake .. -DROBOT_MODEL_RBDL=ON -DSOLVER_PROXQP=ON -DSOLVER_EIQUADPROG=ON -DSOLVER_QPSWIFT=ON -DSOLVER_OSQP=ON -DSOLVER_CLARABEL=ON -DCMAKE_BUILD_TYPE=RELEASE
# DAQP
git clone https://github.com/darnstrom/daqp.git
cd daqp
git checkout dc13508
mkdir build && cd build
cmake ..
make -j8 && sudo make install && cd ../..

# WBC
mkdir wbc/build && cd wbc/build
cmake .. -DROBOT_MODEL_RBDL=ON -DSOLVER_PROXQP=ON -DSOLVER_EIQUADPROG=ON -DSOLVER_QPSWIFT=ON -DSOLVER_OSQP=ON -DSOLVER_DAQP=ON -DCMAKE_BUILD_TYPE=RELEASE
make -j8 && sudo make install && cd ..

sudo ldconfig
