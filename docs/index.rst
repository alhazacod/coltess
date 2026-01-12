Coltess Documentation
=====================

Curves of Light from TESS (COLTESS): photometry tool for TESS FFIs.

Features
--------

* Aperture photometry with local background subtraction
* Gaia DR3 catalog integration
* Parallel processing for thousands of images
* Lomb-Scargle periodogram analysis

Installation
------------

From PyPI (Not released yet)::

    pip install coltess

From source::

    git clone https://github.com/alhazacod/coltess.git
    cd coltess
    pip install -e .

Quick Example
-------------

.. code-block:: python

   from coltess import create_catalog, process_images_parallel
   
   # Create catalog
   star = create_catalog("lambda tau", output_file="catalog.csv")
   
   # Process images in parallel
   process_images_parallel(
       script_file="sector_5.sh",
       catalog_file="catalog.csv",
       output_dir="results",
       star=star
   )

For detailed examples, see the README on GitHub: https://github.com/alhazacod/coltess

API Reference
=============

Core Module
-----------

.. automodule:: coltess.core
   :members:
   :show-inheritance:

Photometry Module
-----------------

.. automodule:: coltess.photometry
   :members:
   :show-inheritance:

Catalog Module
--------------

.. automodule:: coltess.catalog
   :members:
   :show-inheritance:

Download Module
---------------

.. automodule:: coltess.download
   :members:
   :show-inheritance:

Analysis Module
---------------

.. automodule:: coltess.analysis
   :members:
   :show-inheritance:

Parallel Module
---------------

.. automodule:: coltess.parallel
   :members:
   :show-inheritance:

Indices
=======

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
