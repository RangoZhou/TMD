#!/bin/bash

tar -jxf boost_1_71_0.tar.bz2
boost_folder=boost_1_71_0
mkdir ${boost_folder}_installation
cd ${boost_folder}
./bootstrap.sh --show-libraries
./bootstrap.sh --prefix=../${boost_folder}_installation --with-libraries=filesystem,system,thread,serialization,program_options
./b2 install 
