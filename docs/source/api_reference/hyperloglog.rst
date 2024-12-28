HyperLogLog API
===============

The `HyperLogLog` class provides functionality for cardinality estimation and similarity measurements using HyperLogLog sketches. This class supports operations such as updating, merging, and comparing sketches.

C++ Interface
-------------

.. cpp:function:: HyperLogLog(int np)

   Constructs a `HyperLogLog` object with the specified precision.

   **Parameters**:
     - `np` (*int*): The number of bits used for precision. \(2^{np}\) is the size of the sketch.

.. cpp:function:: void update(char* seq)

   Updates the HyperLogLog sketch with a given sequence.

   **Parameters**:
     - `seq` (*char***): A pointer to the sequence for updating the sketch.

.. cpp:function:: HyperLogLog merge(const HyperLogLog& other) const

   Merges the current sketch with another `HyperLogLog` sketch.

   **Parameters**:
     - `other` (*const HyperLogLog&*): The other `HyperLogLog` object to merge with.

   **Returns**:
     - `HyperLogLog`: A new `HyperLogLog` object resulting from the merge.

.. cpp:function:: void printSketch()

   Prints the content of the HyperLogLog sketch for debugging purposes.

.. cpp:function:: double distance(const HyperLogLog& h2) const

   Computes the distance between two `HyperLogLog` sketches. The distance is defined as \(1.0 - \text{Jaccard index}\).

   **Parameters**:
     - `h2` (*const HyperLogLog&*): The other `HyperLogLog` object to compare.

   **Returns**:
     - `double`: The distance score.

.. cpp:function:: double jaccard_index(const HyperLogLog& h2) const

   Computes the Jaccard index between two `HyperLogLog` sketches.

   **Parameters**:
     - `h2` (*const HyperLogLog&*): The other `HyperLogLog` object to compare.

   **Returns**:
     - `double`: The Jaccard index.

---

Python Interface
----------------

.. py:class:: HyperLogLog(np: int = 10)

   Represents the Python interface for the `HyperLogLog` sketch.

   **Parameters**:
     - `np` (*int*): The number of bits used for precision. Default is 10, resulting in a sketch size of \(2^{10}\).

.. py:method:: update(seq: str)

   Updates the HyperLogLog sketch with a given sequence.

   **Parameters**:
     - `seq` (*str*): The sequence to update the sketch.

.. py:method:: merge(other: HyperLogLog) -> HyperLogLog

   Merges the current sketch with another `HyperLogLog` sketch.

   **Parameters**:
     - `other` (*HyperLogLog*): The other `HyperLogLog` object to merge with.

   **Returns**:
     - `HyperLogLog`: A new `HyperLogLog` object resulting from the merge.

.. py:method:: jaccard(other: HyperLogLog) -> float

   Computes the Jaccard index between two `HyperLogLog` sketches.

   **Parameters**:
     - `other` (*HyperLogLog*): The other `HyperLogLog` object to compare.

   **Returns**:
     - `float`: The Jaccard index.

.. py:method:: distance(other: HyperLogLog) -> float

   Computes the distance between two `HyperLogLog` sketches. The distance is defined as \(1.0 - \text{Jaccard index}\).

   **Parameters**:
     - `other` (*HyperLogLog*): The other `HyperLogLog` object to compare.

   **Returns**:
     - `float`: The distance score.

