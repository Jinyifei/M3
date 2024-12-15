Example Models
==============

On Github, two example models are stored under ``m3/lab/exampe/`` 

1: *a black body ionized pure gas nebula* **(example1)**,

2: *an input stellar ionized dusty nebula* **(example2)**.

Pure gas HII region model
-------------------------

.. rubric:: Set up the model example1
   :name: model-setup1
   
Go to ``/example``, review the input files and make the subdirectory ``result1`` to store the output files.

.. code :: bash

    $ cd ~/M3/lab/example
    
    $ ls example1*
    
    $ mkdir results1
    
.. rubric:: Run the model example1
   :name: run-model1
   
Go back to ``m3/lab``, and run the model example1.

.. code :: bash

    $ cd ../
    
    $ ./m3 < example/example1.in
   

Dusty gas HII region model
--------------------------

.. rubric:: Set up the model example2
   :name: model-setup2
   
Go to ``/example``, review the input files and make the subdirectory ``result2`` to store the output files.

.. code :: bash

    $ cd ~/M3/lab/example
    
    $ ls example2*
    
    $ mkdir results2
    
.. rubric:: Run the model example2
   :name: run-model2
   
Go back to ``m3/lab``, and run the model example1.

.. code :: bash

    $ cd ../
    
    $ ./m3 < example/example2.in
    
    
Glance at the output data
--------------------------

Check the model in ``m3/lab/example/result1``. The 3D distribution of some physical variables are stored in ``structure.out``. Assuming we use ``python`` to analyse the data, the commands are below:

.. code :: bash

    $ cd m3/lab/example/result1
    
    $ python3
    
    $ >>>import numpy as np
        
    $ >>>data = np.loadtxt('structure.out')
    
    $ >>>x = data[:,0] 
    
    $ >>>r = np.sqrt(data[:,0]**2+data[:,1]**2+data[:,2]**2)
    

Plot the electron temperature profile:
    
.. code :: bash

    $ >>>import matplotlib.pyplot as plt
        
    $ >>>plt.plot(r, data[:,3], '.')

Plot the electron density profile:            

.. code :: bash
        
    $ >>>plt.plot(r, data[:,], '.')    
    

Go through the same process for the dusty HII region model.


*The output profiles are below:*

For the **Pure HII region model**

.. figure:: example_1.png
   :width: 800
   
For the **Dusty HII region model**  
   
.. figure:: example_2.png
   :width: 800
   
       
