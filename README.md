# HyINDEL

[![Build Status](https://travis-ci.com/alok123t/HyINDEL.svg?token=4hAKK2irggAzvcM7yK4z&branch=master)](https://travis-ci.com/alok123t/HyINDEL)
[![Cpp Standard](https://img.shields.io/badge/C%2B%2B-11-blue.svg)](https://en.wikipedia.org/wiki/C%2B%2B11)

## Installation
### Requirements
Tested with GCC (minimum tested version 4.9.2 check using gcc -v), C++11 support required
* [CMake](https://cmake.org/download/)
* [Minia](https://github.com/GATB/minia#instructions)
* [Minimap2](https://github.com/lh3/minimap2#getting-started)
* [Mosdepth](https://github.com/brentp/mosdepth#installation)
* [Samtools](https://github.com/samtools/samtools#building-samtools)

```shell
git clone --recursive https://github.com/alok123t/HyINDEL.git
cd HyINDEL && mkdir -p build && cd build
# with root access
cmake .. && make -j 4 install

# without root access, install in a local directory
# Executable path will be /path/to/install/dir/bin/HyINDEL
mkdir -p /path/to/install/dir
cmake .. -DCMAKE_INSTALL_PREFIX=/path/to/install/dir && make -j 4 install
```

Note: 
* Make sure paths to minia, minimap2, mosdepth and samtools are on [PATH](https://www.wikihow.com/Change-the-Path-Variable-in-Linux), otherwise write the absolute paths before installation in `scripts/Assemble.sh`, `scripts/Pre.sh` and `scripts/Post.sh`
* Will install/replace bamtools if already present in installation directory

## Usage
```shell
HyINDEL --help
```

### Test Example
```shell
HyINDEL -i ../test/input.bam -o ../test/output -r ../test/ref.fa -s 350 -d 20 -l 100 -c 30 -t 4
```
Compare output with file `test/expected-output.vcf`

### Input/Output
* Input file: Coordinate sorted BAM file with index (.bai)
* Output file: `output.vcf` in output directory

Output directory will be created, if it doesn't exist

### Parameters
| Options Short | Options Long | Description | Attributes | Mandatory |
| --- | --- | --- | --- | --- |
| `-i PATH` | `--inp=PATH` | Input File | Path | <ul><li>[x] yes</li></ul> |
| `-o PATH` | `--out=PATH` | Output Folder | Path | <ul><li>[x] yes</li></ul> |
| `-s VAL` | `--insSz=VAL` | Insert Size | Integer | <ul><li>[x] yes</li></ul> |
| `-d VAL` | `--stdDev=VAL` | Standard Deviation | Integer | <ul><li>[x] yes</li></ul> |
| `-l VAL` | `--readLen=VAL` | Read Length | Integer | <ul><li>[x] yes</li></ul> |
| `-c VAL` | `--cov=VAL` | Coverage | Integer | <ul><li>[x] yes</li></ul> |
| `-t VAL` | `--threads=VAL` | Threads | Integer | <ul><li>[ ] no</li></ul> |
