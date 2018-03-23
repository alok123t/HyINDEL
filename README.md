# indel-detect

## Requirements
* [Bamtools](https://github.com/pezmaster31/bamtools)
* [Intel TBB](https://github.com/01org/tbb/releases)

## Installation
Run `make`

## Usage
`./Cluster <input BAM> <insert size> <std dev> > <output file>`
### Example
`./Cluster Data/NA12878/Alignment/chr1.bam 300 10 > Data/chr1.txt`
