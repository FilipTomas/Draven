# Draven

Draven is a de novo genome assembler for long uncorrected or corrected reads.

## Usage
To build Draven run the following commands (< 30s):

```bash
git clone https://github.com/lbcb-sci/raven && cd raven && mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release .. && make
```
To use a specific model use:
```bash
cp models/model.hpp /third_party/catboost_model/catboost_model.hpp
```
Then rebuild.

which will create raven executable and unit tests (running `make install` will install the executable to your system). Running the executable will display the following usage:

```bash
usage: raven [options ...] <sequences>

  # default output is to stdout in FASTA format
  <sequences>
    input file in FASTA/FASTQ format (can be compressed with gzip)

  options:
    -p, --polishing-rounds <int>
      default: 2
      number of times racon is invoked
    -m, --match <int>
      default: 3
      score for matching bases
    -n, --mismatch <int>
      default: -5
      score for mismatching bases
    -g, --gap <int>
      default: -4
      gap penalty (must be negative)
    --graphical-fragment-assembly <string>
      prints the assembly graph in GFA format
    --resume
      resume previous run from last checkpoint
    --disable-checkpoints
      disable checkpoint file creation
    -t, --threads <int>
      default: 1
      number of threads
    --version
      prints the version number
    -h, --help
      prints the usage


#### Dependencies
- gcc 4.8+ | clang 4.0+
- cmake 3.11+
- zlib 1.2.8+

###### Hidden
- lbcb-sci/racon/tree/library 3.0.1
- rvaser/bioparser 3.0.13
- (racon_test) google/googletest 1.10.0
