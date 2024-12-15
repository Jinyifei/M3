Getting Started
================

The updated M3 code can be downloaded on Github at https://github.com/Jinyifei/M3/ .
After downloading and unziping the code, changing to the ``/M3`` directory.
There are two subdirectories, ``/src`` and ``/lab``, under the parent directory ``/M3``.


.. rubric:: Before you start
   :name: before-start

M3 is written in Fortran. A morden Fortran compiler is required.
M3 is a MPI code, a MPI standard, MPICH, is required.

.. rubric:: Compile the code
   :name: compile-code

The ``Makefile`` is provided under ``M3/src``, compiled by ``mpif77``.

.. code :: bash

    $ cd ~/M3/src/
    
    $ make build
