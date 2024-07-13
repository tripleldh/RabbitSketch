![RabbitSketch](sketch.png)

## Getting Started

A Linux system on a recent x86_64 CPU is required.

### Installing (C++ interface) 


```bash
cd RabbitSketch
mkdir build
cd build
cmake -DCXXAPI=ON .. -DCMAKE_INSTALL_PREFIX=.
make
make install
export LD_LIBRARY_PATH=`pwd`/lib:$LD_LIBRARY_PATH
```


### Testing (C++)

```bash
cd ../examples/
#default install dir: ../build/
make 
./exe_main genome.fna
```

We will get the value of jaccard and distance.

or:
```bash
./exe_SKETCH_ALGORITHM FILE_PATH threshold(0.05) thread_num 
```
We will get the distance among large-scale genome sequences.

### PYTHON bind
**pip install:**
```bash
cd RabbitSketch
pip3 install . --user
```
or
```bash
#pypi available (not up to date)
#pip3 install rabbitsketch --user
```
**cmake install**
```bash
cd RabbitSketch
mkdir build
cd build
cmake .. #default with pybind support
make
```
**test using bpython or python**

```bash
cd examples
python pysketch.py #require fastx
```
We will get the Jaccard index among large-scale genome sequences with Python API. To change the algorithm, simply modify 
```bash 
sketch.SKETCH_ALGORITHM.
```
** case study for multi-thread sketch building with Python API
```bash
python multi_minhash.py #require pymp

