MinHash API
===========

The MinHash class provides a powerful sketching algorithm for estimating set similarity and distances. Below is the detailed documentation for both C++ and Python interfaces.

Contents:
---------

- **C++ Interface**: Methods and attributes provided for C++ development.
- **Python Interface**: Methods and attributes exposed to Python using pybind11.

C++ Interface
-------------

The C++ interface provides direct access to the MinHash class with full control over the sketching process.

### Methods

.. cpp:function:: void update(char * seq)

   Updates the MinHash with a new sequence. This method supports streaming updates.

   **Parameters**:
     - `seq`: A pointer to the sequence to be added.

.. cpp:function:: void merge(MinHash& msh)

   Merges another MinHash into this MinHash, combining their hash values.

   **Parameters**:
     - `msh`: The MinHash object to be merged.

.. cpp:function:: double jaccard(MinHash * msh)

   Computes the Jaccard index between this MinHash and another.

   **Parameters**:
     - `msh`: The MinHash object to compare against.

   **Returns**:
     - `double`: The Jaccard index.

.. cpp:function:: double distance(MinHash * msh)

   Computes the mutation distance between this MinHash and another, as defined in Mash.

   **Parameters**:
     - `msh`: The MinHash object to compare against.

   **Returns**:
     - `double`: The mutation distance.

.. cpp:function:: MashLite toLite(vector<uint64_t> hashL) const

   Converts the current MinHash object into a `MashLite` representation, for index_dict method.


.. cpp:function:: void printMinHashes()

   Prints the MinHash values for debugging purposes.

.. cpp:function:: uint64_t getTotalLength()

   Returns the total sequence length, including multiple updates.

   **Returns**:
     - `uint64_t`: The total sequence length.

### Attributes

.. cpp:member:: int getKmerSize()

   Returns the k-mer size used by this MinHash.

   **Returns**:
     - `int`: The k-mer size.

.. cpp:member:: uint32_t getSeed()

   Returns the hash seed used for generating the MinHash.

   **Returns**:
     - `uint32_t`: The hash seed.

.. cpp:member:: uint32_t getMaxSketchSize()

   Returns the maximum sketch size.

   **Returns**:
     - `uint32_t`: The maximum sketch size.

.. cpp:member:: bool isEmpty()

   Checks if the MinHash is empty.

   **Returns**:
     - `bool`: `true` if the MinHash is empty, `false` otherwise.

---

Python Interface
----------------

The Python interface exposes MinHash functionality via pybind11, enabling easy use in Python projects.

### Constructor

.. py:class:: MinHash(kmer=21, size=1000, seed=42)

   Creates a MinHash object with the specified parameters.

   **Parameters**:
     - `kmer` (int): Size of the k-mers (default: 21).
     - `size` (int): Maximum number of hashes to store (default: 1000).
     - `seed` (int): Random seed for reproducibility (default: 42).

### Methods

.. py:method:: update(seq: str)

   Updates the MinHash with a new sequence.

   **Parameters**:
     - `seq` (str): The sequence to add.

.. py:method:: merge(other: MinHash)

   Merges another MinHash into this MinHash, combining their hash values.

   **Parameters**:
     - `other` (MinHash): The MinHash object to merge.

.. py:method:: jaccard(other: MinHash) -> float

   Computes the Jaccard index between this MinHash and another.

   **Parameters**:
     - `other` (MinHash): The MinHash object to compare against.

   **Returns**:
     - `float`: The Jaccard index.

.. py:method:: distance(other: MinHash) -> float

   Computes the mutation distance between this MinHash and another.

   **Parameters**:
     - `other` (MinHash): The MinHash object to compare against.

   **Returns**:
     - `float`: The mutation distance.

.. py:method:: get_total_length() -> int

   Returns the total sequence length, including multiple updates.

   **Returns**:
     - `int`: The total sequence length.

.. py:method:: print_min_hashes()

   Prints the MinHash values for debugging purposes.

.. py:method:: count() -> int

   Estimates the cardinality count of the set represented by the MinHash.

   **Returns**:
     - `int`: The estimated cardinality count.

### Attributes

.. py:attribute:: kmer_size

   Returns the k-mer size used by this MinHash.

   **Returns**:
     - `int`: The k-mer size.

.. py:attribute:: seed

   Returns the hash seed used for generating the MinHash.

   **Returns**:
     - `int`: The hash seed.

.. py:attribute:: max_sketch_size

   Returns the maximum sketch size.

   **Returns**:
     - `int`: The maximum sketch size.

.. py:attribute:: is_empty

   Checks if the MinHash is empty.

   **Returns**:
     - `bool`: `True` if the MinHash is empty, `False` otherwise.

---

Examples
--------

**C++ Example**:

.. code-block:: cpp

   #include "MinHash.h"

   MinHash mh;
   mh.update("ACGTACGT");
   std::vector<uint64_t> hashes = mh.storeMinHashes();

   MinHash other;
   other.update("TGCA");

   double jaccardIndex = mh.jaccard(&other);
   std::cout << "Jaccard Index: " << jaccardIndex << std::endl;

**Python Example**:

.. code-block:: python

   from rabbitsketch import MinHash

   mh = MinHash(kmer=21, size=1000, seed=42)
   mh.update("ACGTACGT")

   other = MinHash(kmer=21)
   other.update("TGCA")

   jaccard_index = mh.jaccard(other)
   print(f"Jaccard Index: {jaccard_index}")

