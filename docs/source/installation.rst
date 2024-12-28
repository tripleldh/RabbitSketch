Installation
============

How to install rabbitsketch.


Install the C++ interface, follow these steps:

.. code-block:: bash

   cd RabbitSketch
   mkdir build
   cd build
   cmake -DCXXAPI=ON .. -DCMAKE_INSTALL_PREFIX=.
   make
   make install
   export LD_LIBRARY_PATH=`pwd`/lib:$LD_LIBRARY_PATH


To install the Python interface using `pip`, follow these steps:

.. code-block:: bash

   cd RabbitSketch
   pip3 install . 
