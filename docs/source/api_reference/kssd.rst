Kssd API
========

The `Kssd` class provides functionalities for processing and analyzing k-mer based sketches, with methods for updating, storing, comparing, and managing sketches.


C++ Interface
-------------

.. cpp:function:: kssd_parameter_t(int half_k_, int half_subk_, int drlevel_, string shuffle_file)

   Constructs a `kssd_parameter_t` object with the specified parameters.

   **Parameters**:
     - `half_k_` (*int*): Half of the k-mer size used for hashing.
     - `half_subk_` (*int*): Half of the sub k-mer size for dimensionality reduction.
     - `drlevel_` (*int*): Dimensionality reduction level.
     - `shuffle_file` (*string*): Path to the file containing the shuffled dimension map.


.. cpp:function:: Kssd(kssd_parameter_t params)

   Constructs a `Kssd` object with the specified parameters.

   **Parameters**:
     - `params` (*kssd_parameter_t*): A parameter object containing configuration for `Kssd`.

.. cpp:function:: void update(const char* seq)

   Updates the `Kssd` object with the provided sequence.

   **Parameters**:
     - `seq` (*const char**): A pointer to the sequence to be processed.

.. cpp:function:: vector<uint32_t> storeHashes()

   Retrieves the stored 32-bit hashes from the `Kssd` object.

   **Returns**:
     - `vector<uint32_t>`: A vector of 32-bit hashes.

.. cpp:function:: vector<uint64_t> storeHashes64()

   Retrieves the stored 64-bit hashes from the `Kssd` object.

   **Returns**:
     - `vector<uint64_t>`: A vector of 64-bit hashes.

.. cpp:function:: void loadHashes(std::vector<uint32_t> hashArr)

   Loads a list of 32-bit hashes into the `Kssd` object.

   **Parameters**:
     - `hashArr` (*std::vector<uint32_t>*): A vector of 32-bit hashes to load.

.. cpp:function:: void loadHashes64(std::vector<uint64_t> hashArr)

   Loads a list of 64-bit hashes into the `Kssd` object.

   **Parameters**:
     - `hashArr` (*std::vector<uint64_t>*): A vector of 64-bit hashes to load.

.. cpp:function:: double jaccard(Kssd* kssd)

   Computes the Jaccard similarity between this `Kssd` object and another.

   **Parameters**:
     - `kssd` (*Kssd***): A pointer to another `Kssd` object to compare against.

   **Returns**:
     - `double`: The Jaccard similarity value.

.. cpp:function:: double distance(Kssd* kssd)

   Computes the Mash distance between this `Kssd` object and another.

   **Parameters**:
     - `kssd` (*Kssd***): A pointer to another `Kssd` object to compare against.

   **Returns**:
     - `double`: The Mash distance value.

.. cpp:function:: void transSketches(vector<KssdLite>& sketches, sketchInfo_t& info, string dictFile, string indexFile, int numThreads)

   Transforms `KssdLite` sketches for indexing purposes.

   **Parameters**:
     - `sketches` (*vector<KssdLite>&*): A reference to the vector of `KssdLite` sketches.
     - `info` (*sketchInfo_t&*): Metadata associated with the sketches.
     - `dictFile` (*string*): Path to the dictionary file.
     - `indexFile` (*string*): Path to the index file.
     - `numThreads` (*int*): The number of threads to use.

.. cpp:function:: void index_tridist(vector<KssdLite>& sketches, sketchInfo_t& info, string refSketchOut, string outputFile, int kmer_size, double maxDist, int isContainment, int numThreads)

   Computes the sketch index using a distance-based method.

   **Parameters**:
     - `sketches` (*vector<KssdLite>&*): A reference to the vector of `KssdLite` sketches.
     - `info` (*sketchInfo_t&*): Metadata associated with the sketches.
     - `refSketchOut` (*string*): Path to save the reference sketches.
     - `outputFile` (*string*): Path to the output file.
     - `kmer_size` (*int*): Size of the k-mers.
     - `maxDist` (*double*): Maximum allowed distance.
     - `isContainment` (*int*): Whether to use containment comparison.
     - `numThreads` (*int*): The number of threads to use.

.. cpp:function:: void saveSketches(vector<KssdLite>& sketches, sketchInfo_t& info, string outputFile)

   Saves the `KssdLite` sketches to a specified file.

   **Parameters**:
     - `sketches` (*vector<KssdLite>&*): A reference to the vector of `KssdLite` sketches.
     - `info` (*sketchInfo_t&*): Metadata associated with the sketches.
     - `outputFile` (*string*): Path to the output file.

.. cpp:function:: bool isSketchFile(string inputFile)

   Checks if the given file has a `.sketch` extension.

   **Parameters**:
     - `inputFile` (*string*): The file name to check.

   **Returns**:
     - `bool`: `true` if the file is a `.sketch` file, otherwise `false`.

.. cpp:function:: void printInfos(vector<Sketch::Kssd*>& sketches, string outputFile)

   Prints basic information about the sketches to the specified file.

   **Parameters**:
     - `sketches` (*vector<Sketch::Kssd*>&*): A reference to the vector of `Sketch::Kssd` objects.
     - `outputFile` (*string*): Path to the output file.

.. cpp:function:: void printSketches(vector<Sketch::Kssd*>& sketches, string outputFile)

   Prints detailed information about the sketches, including all stored hashes.

   **Parameters**:
     - `sketches` (*vector<Sketch::Kssd*>&*): A reference to the vector of `Sketch::Kssd` objects.
     - `outputFile` (*string*): Path to the output file.

.. cpp:function:: KssdLite toLite() const

   Converts the current `Kssd` object into a lightweight `KssdLite` representation.

   **Returns**:
     - `KssdLite`: A lightweight representation of the `Kssd` object.
 ---

Python Interface
-------------

.. cpp:function:: Kssd(kssd_parameter_t params)

   Constructs a `Kssd` object with the specified parameters.

   **Parameters**:
     - `params` (*kssd_parameter_t*): A parameter object containing configuration for `Kssd`.

.. py:class:: Kssd(params: kssd_parameter_t)

   Constructs a `Kssd` object in Python.

   **Parameters**:
     - `params` (*kssd_parameter_t*): Configuration parameters for the `Kssd` object.

.. cpp:function:: void update(const char* seq)

   Updates the `Kssd` object with the provided sequence.

   **Parameters**:
     - `seq` (*const char**): A pointer to the sequence to be processed.

.. py:method:: update(seq: str)

   Updates the sketch with the given sequence.

   **Parameters**:
     - `seq` (*str*): The sequence to be processed.

.. cpp:function:: double jaccard(Kssd* kssd)

   Computes the Jaccard similarity between this `Kssd` object and another.

   **Parameters**:
     - `kssd` (*Kssd***): A pointer to another `Kssd` object to compare against.

   **Returns**:
     - `double`: The Jaccard similarity value.

.. py:method:: jaccard(other: Kssd) -> float

   Computes the Jaccard similarity between this sketch and another.

   **Parameters**:
     - `other` (*Kssd*): The other `Kssd` object to compare.

   **Returns**:
     - `float`: The Jaccard similarity value.

.. cpp:function:: double distance(Kssd* kssd)

   Computes the Mash distance between this `Kssd` object and another.

   **Parameters**:
     - `kssd` (*Kssd***): A pointer to another `Kssd` object to compare against.

   **Returns**:
     - `double`: The Mash distance value.

.. py:method:: distance(other: Kssd) -> float

   Computes the Mash distance between this sketch and another.

   **Parameters**:
     - `other` (*Kssd*): The other `Kssd` object to compare.

   **Returns**:
     - `float`: The Mash distance value.

.. cpp:function:: vector<uint32_t> storeHashes()

   Retrieves the stored 32-bit hashes from the `Kssd` object.

   **Returns**:
     - `vector<uint32_t>`: A vector of 32-bit hashes.

.. py:method:: store_hashes() -> list[int]

   Retrieves the stored 32-bit hashes.

   **Returns**:
     - `list[int]`: A list of 32-bit hashes.

.. cpp:function:: KssdLite toLite() const

   Converts the current `Kssd` object into a lightweight `KssdLite` representation.

   **Returns**:
     - `KssdLite`: A lightweight representation of the `Kssd` object.

.. py:method:: toLite() -> KssdLite

   Converts the current `Kssd` object to a lightweight `KssdLite` representation.

   **Returns**:
     - `KssdLite`: A lightweight representation of the `Kssd` object.

.. cpp:function:: kssd_parameter_t(int half_k_, int half_subk_, int drlevel_, string shuffle_file)

   Constructs a `kssd_parameter_t` object with the specified parameters.

   **Parameters**:
     - `half_k_` (*int*): Half of the k-mer size used for hashing.
     - `half_subk_` (*int*): Half of the sub k-mer size for dimensionality reduction.
     - `drlevel_` (*int*): Dimensionality reduction level.
     - `shuffle_file` (*string*): Path to the shuffle file containing dimension map.

.. py:class:: kssd_parameter_t(half_k: int, half_subk: int, drlevel: int, shuffle_file: str)

   Represents the parameter set used for configuring the `Kssd` object.

   **Parameters**:
     - `half_k` (*int*): Half of the k-mer size used for hashing.
     - `half_subk` (*int*): Half of the sub k-mer size for dimensionality reduction.
     - `drlevel` (*int*): Dimensionality reduction level.
     - `shuffle_file` (*str*): Path to the shuffle file containing dimension map.


.. py:function:: save_sketches(sketches: list[KssdLite], info: SketchInfo, filename: str)

   Saves a list of `KssdLite` sketches to a specified file.

   **Parameters**:
     - `sketches` (*list[KssdLite]*): A list of sketches to save.
     - `info` (*SketchInfo*): Metadata associated with the sketches.
     - `filename` (*str*): Path to the file where the sketches will be saved.

.. py:function:: trans_sketches(sketches: list[KssdLite], info: SketchInfo, dict_file: str, index_file: str, num_threads: int)

   Transforms `KssdLite` sketches into a format suitable for indexing.

   **Parameters**:
     - `sketches` (*list[KssdLite]*): A list of sketches to transform.
     - `info` (*SketchInfo*): Metadata associated with the sketches.
     - `dict_file` (*str*): Path to the dictionary file used for transformation.
     - `index_file` (*str*): Path to the index file used for transformation.
     - `num_threads` (*int*): Number of threads for parallel processing.

.. py:function:: index_dict(sketches: list[KssdLite], info: SketchInfo, ref_sketch_out: str, output_file: str, kmer_size: int, max_dist: float, is_containment: bool, num_threads: int)

   Computes the sketch index using a distance-based dictionary approach.

   **Parameters**:
     - `sketches` (*list[KssdLite]*): A list of sketches to process.
     - `info` (*SketchInfo*): Metadata associated with the sketches.
     - `ref_sketch_out` (*str*): Path to save the reference sketches.
     - `output_file` (*str*): Path to the file where results will be saved.
     - `kmer_size` (*int*): Size of the k-mers used for sketching.
     - `max_dist` (*float*): Maximum allowed distance for comparisons.
     - `is_containment` (*bool*): Whether to use containment comparisons.
     - `num_threads` (*int*): Number of threads for parallel processing.
