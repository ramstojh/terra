*terra* package
---------------
*terra* is a tool for reproducing the abundance patterns of stars. *terra* can be used to modelate planet engulfment events. For more details, see our papers  `Yana Galarza et al. 2016 <https://ui.adsabs.harvard.edu/abs/2016A%26A...589A..65G/abstract>`_ and  `Yana Galarza et al. 2021 <https://ui.adsabs.harvard.edu/abs/2021arXiv210900679G>`_

Installation
------------
The only way to install *terra* is through pip::

    pip install terra-2.0

Dependencies
------------
The main dependencies of *terra* are `pandas <https://pandas.pydata.org/>`_, `NumPy <https://numpy.org/>`_, `Astropy <https://www.astropy.org/>`_, `matplotlib <https://matplotlib.org/>`_, `tqdm <https://tqdm.github.io/>`_, and `os <https://docs.python.org/3/library/os.html>`_. 
These are installed using pip::

    pip install pandas numpy astropy matplotlib tqdm numpy 
    
Example usage
-------------

.. code-block:: python

    from terra import pacha
    
    # Computing the convective mass of a star with [Fe/H] = 0.164 dex 
    # and mass = 1.18 solar masses
    # The mass can take values from 0.5 <= M <= 1.3 (solar mass)
    # The [Fe/H] can take values from -1.0 <= [Fe/H] <= 0.3 (dex)
    # By default the code computes the convective mass using the Yale isocrhones of stellar evolution
    terra.cvmass(feh=0.1, mass=1)
    
    # Computing the abundance pattern of a star with [Fe/H] = 0.164 dex and mass = 1.14 M_sun
    # obs_abd.csv is a table containing the observed abundance.
    terra.pacha(feh=0.164, mass=1 data_input='obs_abd.csv')
    
    # If you want to save the outputs (figures and tables) with the star name (e.g., HIP 71726).
    terra.pacha(feh=0.164, mass=1, Mcon=0.01, data_input='obs_abd.csv', data_output='HIP71726')
    
    #For more details please see the terra's tutorial within the file terra_example file
    

Contributing
------------
*terra* is a tool that needs input to improve. Please contact me (ramstojh@alumni.usp.br) if you have questions. Users are welcome to propose new features or report bugs by opening an issue on GitHub. A special thanks to my friend Marilia Carlos, who created the 'yale.txt' table to estimate convective masses.


Author
------
- `Jhon Yana Galarza <https://github.com/ramstojh>`_


Preferred citation
------------------
`Yana Galarza et al. 2016, A&A, 589, A65 <https://ui.adsabs.harvard.edu/abs/2016A%26A...589A..65G/abstract>`_


License & attribution
---------------------

Copyright 2021, Jhon Yana Galarza.

The source code is made available under the terms of the MIT license.

If you make use of this code, please cite this package and its dependencies.
