Welcome to RainCheckPy documentation
===================================

**RainCheckPy** is a Python package for checking rain data.

Overview
--------
RainCheckPy provides functions to check the quality of rain data time series. Each function returns timeseries of quality indices.
The quality checks are divided into two main types:

#. Single series checks. These checks use only data associated with a single site.
#. Multi-series checks. These checks compare data with a second site's data.

Contents
--------

.. toctree::

   installation
   Single Series Checks <singleserieschecks>
   Multi-series Checks <multiserieschecks>
