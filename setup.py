import os
import subprocess
from unittest import TextTestRunner, TestLoader
from distutils.core import setup, Extension, Command
from distutils.command.build import build as DistutilsBuild
from distutils.command.install import INSTALL_SCHEMES

# need to make sure data_files are installed alongside site-packages
# see http://www.arthurkoziel.com/2008/12/31/
#            including-additional-files-distutils-python-23/
for scheme in INSTALL_SCHEMES.values():
    scheme['data'] = scheme['purelib']

flamingo_dir = os.path.join(".", "src", "flamingo-4.1", "src")


lib_files = [os.path.join(flamingo_dir, "filtertree", "build",
                            "libfiltertree-lib." + ext)
                    for ext in ["so", "dylib"]]
lib_files = [l for l in lib_files if os.path.exists(l)]
data_files = [('.', [lib_files[0]])] if len(lib_files) > 0 else []

module1 = Extension('flamingo',
                    sources=['src/flamingo.cpp'],
                    include_dirs=[os.path.join(flamingo_dir, d) for d in
                                        ["", "filtertree",
                                         os.path.join("filtertree", "src")]],
                    library_dirs=[os.path.join(flamingo_dir, d, "build")
                                      for d in ["common", "util", "listmerger",
                                                "filtertree"]],
                    #runtime_library_dirs=[os.path.join(flamingo_dir,
                    #                                  "filtertree", "build")],
                    runtime_library_dirs=["$ORIGIN"],
                    libraries=['filtertree-lib'])


class BarNoneBuild(DistutilsBuild):
    def run(self):
        """compile appropriate FLAMINGO library before installing"""
        os.chdir(os.path.join("src", "flamingo-4.1", "src", "filtertree"))
        subprocess.call(["cmake", '"-DCMAKE_OSX_ARCHITECTURES=i386;x86_64"',
                            "."])
        subprocess.call(["make"])
        os.chdir(os.path.join("..", "..", "..", ".."))
        DistutilsBuild.run(self)


class TestCommand(Command):
    """Run test suite using 'python setup.py test'"""
    user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        """Run test suite in BarNone.tests"""
        try:
            from BarNone import tests
        except ImportError:
            raise Exception("BarNone is not installed, cannot run tests")
        tests = TestLoader().loadTestsFromNames(["BarNone.tests"])
        t = TextTestRunner(verbosity=1)
        t.run(tests)

setup(name="BarNone",
      author="David Robinson",
      author_email="dgrtwo@princeton.edu",
      description="Match and count barcodes from a barcoded sequencing run, " +
                  "allowing inexact matching.",
      version="0.1",
      packages=["BarNone"],
      package_dir={"BarNone": os.path.join("src", "BarNone")},
      scripts=[os.path.join("scripts", "BarNone")],
      ext_modules=[module1],
      cmdclass={"build": BarNoneBuild, "test": TestCommand},
      requires=["Levenshtein"],
      data_files=data_files
      )
