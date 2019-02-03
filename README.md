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
```

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
bash scripts/preProcess.sh ../test/test_filter.bam output
# Run program
bin/indel -s 300 -d 3 -i ../test/test_filter.bam -o output
# Post-process
bash scripts/postProcess.sh ../test/test_filter.bam output 30 15 10 15 20
```

| Options Short | Options Long | Description | Attributes | Mandatory |
| --- | --- | --- | --- | --- |
| `-s VAL` | `--insSz=VAL` | Insert Size | Integer | <ul><li>[x] yes</li></ul> |
| `-d VAL` | `--stdDev=VAL` | Standard Deviation | Integer | <ul><li>[x] yes</li></ul> |
| `-i PATH` | `--inp=PATH` | Input Files | Absolute path, comma seperated | <ul><li>[x] yes</li></ul> |
| `-o PATH` | `--out=PATH` | Output Folder | Path, folder should exist | <ul><li>[ ] no</li></ul> |
| `-t VAL` | `--threads=VAL` | Threads | Integer, minimum `1` | <ul><li>[ ] no</li></ul> |
| `-v VAL` | `--verbose=VAL` | Verbose | Bool, `0` or `1` | <ul><li>[ ] no</li></ul> |
