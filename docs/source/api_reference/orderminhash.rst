OrderMinHash API
================

The `OrderMinHash` class provides functionality for sketching sequences and computing similarities or distances between them. This class is particularly suited for biological sequence analysis.

C++ Interface
-------------

.. cpp:function:: OrderMinHash()

   Constructs an `OrderMinHash` object with default parameters.

.. cpp:function:: void buildSketch(char* seqNew)

   Builds the `OrderMinHash` sketch. If `seqNew` is `NULL`, the sketch is rebuilt using the existing data. This is useful when parameters are changed, and a new sketch needs to be created.

   **Parameters**:
     - `seqNew` (*char***): A pointer to the sequence for sketching. If `NULL`, rebuilds the sketch with existing data.

.. cpp:function:: double similarity(OrderMinHash& omh2)

   Computes the similarity between two `OrderMinHash` sketches. This is a proxy for edit distance and does not calculate Jaccard similarity.

   **Parameters**:
     - `omh2` (*OrderMinHash&*): The other `OrderMinHash` object to compare.

   **Returns**:
     - `double`: The similarity score.

.. cpp:function:: double distance(OrderMinHash& omh2)

   Computes the distance between two `OrderMinHash` sketches. The distance is defined as \(1.0 - \text{similarity}\).

   **Parameters**:
     - `omh2` (*OrderMinHash&*): The other `OrderMinHash` object to compare.

   **Returns**:
     - `double`: The distance score.

.. cpp:function:: void setK(int k)

   Sets the k-mer size parameter. Default is 21.

   **Parameters**:
     - `k` (*int*): The k-mer size.

.. cpp:function:: void setL(int l)

   Sets the `l` parameter, typically between 2 and 5. Default is 2.

   **Parameters**:
     - `l` (*int*): The `l` parameter value.

.. cpp:function:: void setM(int m)

   Sets the `m` parameter. Default is 500.

   **Parameters**:
     - `m` (*int*): The `m` parameter value.

.. cpp:function:: void setSeed(uint64_t seedNew)

   Sets the seed value for the random generator. Default is 32.

   **Parameters**:
     - `seedNew` (*uint64_t*): The new seed value.

.. cpp:function:: void setReverseComplement(bool isRC)

   Specifies whether to deal with reverse complement sequences. Default is `false`.

   **Parameters**:
     - `isRC` (*bool*): `true` to consider reverse complements, `false` otherwise.

.. cpp:function:: int getK()

   Returns the k-mer size parameter.

   **Returns**:
     - `int`: The k-mer size.

.. cpp:function:: int getL()

   Returns the `l` parameter value.

   **Returns**:
     - `int`: The `l` parameter value.

.. cpp:function:: int getM()

   Returns the `m` parameter value.

   **Returns**:
     - `int`: The `m` parameter value.

.. cpp:function:: uint64_t getSeed()

   Returns the random generator seed value.

   **Returns**:
     - `uint64_t`: The seed value.

.. cpp:function:: bool isReverseComplement()

   Indicates whether reverse complement sequences are considered.

   **Returns**:
     - `bool`: `true` if reverse complements are considered, `false` otherwise.

---

Python Interface
----------------

.. py:class:: OrderMinHash

   Represents the Python interface for `OrderMinHash`.

.. py:method:: build_sketch(seqNew: Optional[str])

   Builds the `OrderMinHash` sketch. If `seqNew` is `None`, rebuilds the sketch using existing data.

   **Parameters**:
     - `seqNew` (*Optional[str]*): The sequence to sketch. If `None`, rebuilds the sketch with existing data.

.. py:method:: similarity(omh2: OrderMinHash) -> float

   Computes the similarity between two `OrderMinHash` sketches.

   **Parameters**:
     - `omh2` (*OrderMinHash*): The other `OrderMinHash` object to compare.

   **Returns**:
     - `float`: The similarity score.

.. py:method:: distance(omh2: OrderMinHash) -> float

   Computes the distance between two `OrderMinHash` sketches.

   **Parameters**:
     - `omh2` (*OrderMinHash*): The other `OrderMinHash` object to compare.

   **Returns**:
     - `float`: The distance score.

.. py:method:: set_k(k: int)

   Sets the k-mer size parameter. Default is 21.

   **Parameters**:
     - `k` (*int*): The k-mer size.

.. py:method:: set_l(l: int)

   Sets the `l` parameter, typically between 2 and 5. Default is 2.

   **Parameters**:
     - `l` (*int*): The `l` parameter value.

.. py:method:: set_m(m: int)

   Sets the `m` parameter. Default is 500.

   **Parameters**:
     - `m` (*int*): The `m` parameter value.

.. py:method:: set_seed(seed: int)

   Sets the seed value for the random generator. Default is 32.

   **Parameters**:
     - `seed` (*int*): The new seed value.

.. py:method:: set_reverse_complement(isRC: bool)

   Specifies whether to deal with reverse complement sequences. Default is `False`.

   **Parameters**:
     - `isRC` (*bool*): `True` to consider reverse complements, `False` otherwise.


