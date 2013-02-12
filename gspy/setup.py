# setup.py

import os
import numpy as np

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

module = 'compeconcy'

setup(cmdclass={'build_ext': build_ext},
      name=module,
      version='1.0',
      ext_modules = [
              Extension("compeconcy",
                        ["compeconcy.pyx"],
                        include_dirs=[np.get_include()])
                    ]
      )
