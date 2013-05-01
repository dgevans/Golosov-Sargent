# setup.py

import os
import numpy

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

module = 'compeconcy'
module2 = 'bellmancy'
module3 = 'bellmancy2'

setup(cmdclass={'build_ext': build_ext},
      name=module,
      version='1.0',
      ext_modules=[Extension(module, [module + ".pyx"])
                  # , Extension(module2, [module2 + ".pyx"])
                   ],
      include_dirs=[numpy.get_include(),
                    os.path.join(numpy.get_include(), 'numpy')]
      )
