![RabbitSketch](sketch.png)
RabbitSketch is a highly optimized sketching library that exploits the power of modern multi-core CPUs. It supports various sketching algorithms including MinHash, OrderMinHash, and HyperLogLog. RabbitSketch achieves significant speedups compared to existing implementations, ranging from 2.30x to 49.55x.In addition, we provide flexible and easy-to-use interfaces for both Python and C++. The similarity analysis of 455GB genomic data can be completed in about 5 minutes using RabbitSketch with Python code.
Detailed API documentation at https://rabbitsketch.readthedocs.io/en/latest
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

If the Kssd algorithm is used, the shuffled file must first be generated. You can generate the shuffled file in the `shuf_file/` directory by running `exe_generate_shuf_file`. Here, L represents the drlevel, and K represents halfk. By default, we use `L3K10.shuf`.

```bash
cd ../examples/
#default install dir: ../build/
make 
#./exe_generate_shuf_file
./exe_SKETCH_ALGORITHM FILE_PATH threshold(0.05) thread_num 
```
We will get the distance among large-scale genome sequences.

```bash
./exe_generate_shuf_file
./exe_main genome1.fna genome2.fna
```
We will get the distance between genome1 and genome2 with different algorithm


### PYTHON bind

## ⚠️ Note on `fastx` Installation
The current version of `fastx` (0.0.3) may fail to install on recent Python versions (e.g., Python 3.10+ or 3.12) due to an invalid `python_requires` specifier in its `setup.py` (`'>=3.5.*'` is not a valid version constraint).
This is due to incompatibility with newer versions of pip and setuptools, which enforce stricter PEP 440 validation.
To work around this issue, you can downgrade `pip` and `setuptools` as recommended in `requirement.txt` before installing.
Python < 3.12 is required.
**pip install:**
```bash
cd RabbitSketch
pip install -r requirement.txt
pip install . --user
```
or
```bash
#python 3.9 is require in this version
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
#in your current conda_env
conda install rabbitsketch
```

**test using bpython or python**

```bash
#pip install -r requirement.txt 
cd examples
python rabbitsketch_pymp.py #require fastx
```
We will get the Jaccard index among large-scale genome sequences with Python API. 


## Tested Platforms and Compilation Fix

We have conducted extensive deployment tests to ensure cross-platform compatibility. The following operating system versions have been successfully tested:

- **Debian**: 9.11, 9.9, 10.2, 11.1, 11.3, 12.0, 12.9  
- **Ubuntu**: 14.04, 16.04, 18.04, 22.04, 24.04  
- **AlmaLinux**: 8.10, 9.5  
- **Rocky Linux**: 8.6, 9.5  
- **CentOS Stream**: 8, 9  
- **CentOS**: 7.6, 7.9  
- **Fedora**: 39, 40  

If you encounter any compatibility issues on these or other platforms, please report them in our GitHub issues section.

### **Fixing `CMake uv_spawn` Failure on Fedora 39**

During our tests, we identified an issue on Fedora 39 where running `cmake` may fail.
This issue is caused by an **incompatible or outdated `libuv`** version provided by the system. Manually compiling and installing the latest `libuv` resolves this problem.

#### **Solution: Manually Compile `libuv`**
1. Install build dependencies:
   ```bash
   sudo dnf install -y autoconf automake libtool gcc gcc-c++

2. Clone and compile the latest `libuv`
   To manually compile and install the latest `libuv`, follow these steps:

   ```bash
   git clone https://github.com/libuv/libuv.git
   cd libuv
   sh autogen.sh
   ./configure --prefix=/usr
   make -j$(nproc)
   sudo make install
   
