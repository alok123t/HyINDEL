# indel-detect

[![Build Status](https://travis-ci.com/alok123t/indel-detect.svg?token=4hAKK2irggAzvcM7yK4z&branch=master)](https://travis-ci.com/alok123t/indel-detect)
[![Cpp Standard](https://img.shields.io/badge/C%2B%2B-11-blue.svg)](https://en.wikipedia.org/wiki/C%2B%2B11)

## Installation
### Requirements
* [CMake](https://cmake.org/download/)
* [Mosdepth](https://github.com/brentp/mosdepth#installation)
* [Samtools](https://github.com/samtools/samtools#building-samtools)

```shell
conda install -c bioconda mosdepth samtools
```

Note: 
* Edit path to mosdepth and samtools before installation in `scripts/preProcess.sh` and `scripts/postProcess.sh` (if installed from source, ignore if installed using conda)
* Will install/replace bamtools if already present in installation directory

```shell
git clone --recursive https://github.com/alok123t/indel-detect.git
cd indel-detect && mkdir -p build && cd build
cmake .. && make -j 4 install
```
For a different install directory, use the following
```shell
cmake .. -DCMAKE_INSTALL_PREFIX=/path/to/install/dir && make -j 4 install
```

## Usage
```shell
plusminus --help
```

### Example
```shell
plusminus -i ../test/input.bam -o ../test/output -s 350 -d 50 -l 100 -c 5 -t 4
```

### Output
`output.vcf` file is created in output directory

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
