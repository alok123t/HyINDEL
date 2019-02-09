# indel-detect

[![Build Status](https://travis-ci.com/alok123t/indel-detect.svg?token=4hAKK2irggAzvcM7yK4z&branch=master)](https://travis-ci.com/alok123t/indel-detect)
[![Cpp Standard](https://img.shields.io/badge/C%2B%2B-11-blue.svg)](https://en.wikipedia.org/wiki/C%2B%2B11)

## Requirements
* [CMake](https://cmake.org/download/)
* [Bedtools](https://bedtools.readthedocs.io/en/latest/content/installation.html)
* [Mosdepth](https://github.com/brentp/mosdepth#installation)

## Installation
```shell
git clone --recursive https://github.com/alok123t/indel-detect.git
cd indel-detect && mkdir -p build && cd build
cmake .. && make -j4 && make test
make install
```
For a different install directory, use the following
```shell
cmake .. -DCMAKE_INSTALL_PREFIX=/path/to/dir
```

### Dependencies
```shell
conda install -c bioconda bedtools samtools
conda install mosdepth
```

## Usage
```shell
bin/indel --help
```

### Example
```shell
cd /path/to/indel-detect
# Create output directory
mkdir output
# Pre-process
bash scripts/preProcess.sh -i ../test/test_filter.bam -o output -m 1
# Run program
bin/indel -s 300 -d 3 -i ../test/test_filter.bam -o output
# Post-process
bash scripts/postProcess.sh -i ../test/test_filter.bam -o output -c 30 -l 15 -m 10 -s 15 -q 20
```

### PreProcess parameters
| Options Short | Description | Attributes | Mandatory |
| --- | --- | --- | --- | 
| `-i PATH` | Input File | Absolute path | <ul><li>[x] yes</li></ul> |
| `-o PATH` | Output Folder | Absolute path | <ul><li>[x] yes</li></ul> |
| `-m VAL` | Use Mosdepth or Samtools for depth | `1` for Mosdepth, `0` for Samtools | <ul><li>[x] yes</li></ul> |

### Detect parameters
| Options Short | Options Long | Description | Attributes | Mandatory |
| --- | --- | --- | --- | --- |
| `-s VAL` | `--insSz=VAL` | Insert Size | Integer | <ul><li>[x] yes</li></ul> |
| `-d VAL` | `--stdDev=VAL` | Standard Deviation | Integer | <ul><li>[x] yes</li></ul> |
| `-i PATH` | `--inp=PATH` | Input Files | Absolute path, comma seperated | <ul><li>[x] yes</li></ul> |
| `-o PATH` | `--out=PATH` | Output Folder | Absolute Path, folder should exist | <ul><li>[x] yes</li></ul> |
| `-t VAL` | `--threads=VAL` | Threads | Integer, minimum `1` | <ul><li>[ ] no</li></ul> |
| `-v VAL` | `--verbose=VAL` | Verbose | Bool, `0` or `1` | <ul><li>[ ] no</li></ul> |

### PostProcess parameters
| Options Short | Description | Attributes | Mandatory |
| --- | --- | --- | --- | 
| `-i PATH` | Input File | Absolute path | <ul><li>[x] yes</li></ul> |
| `-o PATH` | Output Folder | Absolute path | <ul><li>[x] yes</li></ul> |
| `-c VAL` | Coverage | Integer | <ul><li>[x] yes</li></ul> |
| `-l VAL` | Support Large | Integer | <ul><li>[x] yes</li></ul> |
| `-m VAL` | Support imprecise large | Integer | <ul><li>[x] yes</li></ul> |
| `-s VAL` | Support small | Integer | <ul><li>[x] yes</li></ul> |
| `-q VAL` | Mapping quality | Integer | <ul><li>[x] yes</li></ul> |
