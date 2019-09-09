from setuptools import setup, find_packages


setup(name = 'particleinterpolator',
      version = '0.1',
      description = 'Interpolate (in parallel) particle data onto a grid.',
      url = 'https://github.com/rennehan/particle-interpolator',
      author = 'Doug Rennehan',
      author_email = 'douglas.rennehan@gmail.com',
      license = 'GPLv3',
      packages = find_packages(),
      zip_safe = False, install_requires=['numpy', 'mpi4py'])
