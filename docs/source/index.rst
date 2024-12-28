.. rabbitsketch documentation master file, created by
   sphinx-quickstart on 2024-12-27.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to rabbitsketch's Documentation!
=========================================

**rabbitsketch** is a high-performance Python extension library based on pybind11, designed for complex Sketch algorithms and data structures.

Features:
---------
- **save_sketches**: Save Sketches to a file.
- **trans_sketches**: Transform Sketches with multithreading support.
- **index_dict**: Generate an index dictionary for Sketches.
- Provides support for multiple data structures, including `MinHash`, `OSketch`, `OrderMinHash`, `HyperLogLog`, `Kssd`, and more.

.. image:: _static/logo.png
   :align: center
   :alt: rabbitsketch logo

Contents:
---------

.. toctree::
   :maxdepth: 2
   :caption: Documentation Contents:

   installation
   usage
   api_reference
   contributing
   changelog

Installation
------------

To install rabbitsketch, use the following methods:

.. code-block:: bash

   pip install rabbitsketch

Or install from source:
.. code-block:: bash

   git clone https://github.com/yourusername/rabbitsketch.git
   cd rabbitsketch
   pip install .

Quick Start
-----------

Here is a quick example to get started:

.. code-block:: python

   import rabbitsketch

   # Use the save_sketches function
   rabbitsketch.save_sketches(sketches, info, "output_file.sketch")

   # Create a MinHash object and update it with data
   minhash = rabbitsketch.MinHash(kmer=21, size=1000, seed=42)
   minhash.update("ACGTACGT")

   # Compute Jaccard similarity
   similarity = minhash.jaccard(other_minhash)
   print(f"Jaccard similarity: {similarity}")

API Reference
-------------

For detailed documentation on all modules and classes, refer to :ref:`api_reference`.

Contributing
------------

We welcome contributions! Please refer to :ref:`contributing` for contribution guidelines and related information.

Changelog
---------

Check out :ref:`changelog` for the latest updates to the project.

---

Indices and Tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

