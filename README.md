# indel-detect

[![Cpp Standard](https://img.shields.io/badge/C%2B%2B-11-blue.svg)](https://en.wikipedia.org/wiki/C%2B%2B11)

## Requirements
* [CMake](https://cmake.org/download/)

## Installation
```shell
git clone --recursive https://github.com/alok123t/indel-detect.git
cd indel-detect
rm -rf build
mkdir build
cd build
cmake ..
make
```

## Usage
```shell
bin/indel --help
```
### Example
```shell
bin/indel -s 300 -d 3 -i /path/to/test/test.bam
bin/indel -s 300 -d 3 -i ../test/test.bam
```
