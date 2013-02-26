# setup.py

import os
import numpy

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

module = 'compeconcy'
module2 = 'bellmancy'

# TODO: See pandas for how to specify paths inside Extension call
os.chdir('./cyed')  # move into cython directory so setup can find .pyx files

setup(cmdclass={'build_ext': build_ext},
      name=module,
      version='1.0',
      ext_modules=[Extension(module, [module + ".pyx"])
                  , Extension(module2, [module2 + ".pyx"])
                   ],
      include_dirs=[numpy.get_include(),
                    os.path.join(numpy.get_include(), 'numpy')]
      )
