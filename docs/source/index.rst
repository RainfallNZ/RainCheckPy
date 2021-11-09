Welcome to RainCheckPy documentation
===================================

**RainCheckPy** is a Python package for checking rain data.

Overview
--------
RainCheckPy provides functions to check the quality of rain data time series. Each function returns timeseries of quality indices.
The quality checks are divided into two main types:

#. Single series checks. These checks use only data associated with a single site.
#. Multi-series checks. These checks compare data with a second site's data.

Single Series Checks
--------------------

*  Outliers
*  Invalid Number
*  Duplicate Date
*  Dry Spells
*  Repeated Values
*  Homogeneity
*  Rain during sub-zero temperatures
*  Flow events associated with rain events

Mult-series Checks
------------------

*  affinity
*  spearman rank correlation
*  outliers
*  dry spells

Check out the :doc:`usage` section for further information, including
how to :ref:`installation` the project.

.. note::

   This project is under active development.

Contents
--------

.. toctree::
   :maxdepth: 2
   
   usage
   Single Series Checks <singleserieschecks>
   Multi-series Checks <multiserieschecks>
