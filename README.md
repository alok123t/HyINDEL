# indel-detect

[![Cpp Standard](https://img.shields.io/badge/C%2B%2B-17-blue.svg)](https://en.wikipedia.org/wiki/C%2B%2B17)

## Requirements
* [CMake](https://cmake.org/download/)
* [g++](https://launchpad.net/~jonathonf/+archive/ubuntu/gcc-7.3) 7.2 or higher (or) [clang++](http://apt.llvm.org/) 5.0 or higher

[Refer here](#help) for installation help

## Installation
```shell
git clone --recurse-submodules https://github.com/alok123t/indel-detect.git
cd indel-detect
rm -rf build
mkdir build
cd build
cmake ..
make
```

## Usage
```shell
./bin/indel --help
```
### Example
```shell
./bin/indel -s 300 -d 3 -i /path/to/test/test.bam
./bin/indel -s 300 -d 3 -i ../test/test.bam
```

#help
### Ubuntu installation of g++-7
```shell
sudo add-apt-repository ppa:jonathonf/gcc-7.3
sudo apt-get update
```
### Mac installation of g++-7
```shell
brew install gcc@7
```
Ensure g++/clang++ is working

cmake can take a different C++ compiler
```shell
# here /path/to/g++-7 is something like /usr/local/bin/g++-7
# can be found using
# which g++
cmake .. -DCMAKE_CXX_COMPILER=/path/to/g++-7
```
