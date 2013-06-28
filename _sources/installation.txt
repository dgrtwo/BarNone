.. Instructions for installation and requirements

Installation
===================================

Download
-----------------

BarNone can be obtained at the `GitHub repository <http://github.com/dgrtwo/BarNone>`_.

Requirements
------------

BarNone requires the following:

- `Python <http://www.python.org>`_, version >=2.6

  - If using a Python version before 2.7, install the `argparse <http://pypi.python.org/pypi/argparse/>`_ package.

- `CMake <http://www.cmake.org/>`_, a cross-platform open-source installer, used to compile the FLAMINGO C++ library. CMake must have access to a C++ compiler on your platform.

Installation
------------

In the root BarNone directory, perform the following commands::

    python setup.py build
    sudo python setup.py install

If the installation appears successful, you can test your installation of BarNone with::

    python setup.py test