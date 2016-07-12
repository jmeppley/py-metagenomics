try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

DESCRIPTION = "Python scripts for metagenomics"
LONG_DESCRIPTION = open('README.md').read()
NAME = "py-metagenomics"
AUTHOR = "John Eppley"
AUTHOR_EMAIL = "jmeppley@mit.edu"
MAINTAINER = "John Eppley"
MAINTAINER_EMAIL = "jmeppley@mit.edu"
URL = 'http://eddelong.mit.edu/'
DOWNLOAD_URL = 'http://github.com/jmeppley/py-metagenomics'
LICENSE = 'GPL'
VERSION = '0.2.03'

setup(name=NAME,
      version=VERSION,
      description=DESCRIPTION,
      long_description=LONG_DESCRIPTION,
      author=AUTHOR,
      author_email=AUTHOR_EMAIL,
      maintainer=MAINTAINER,
      maintainer_email=MAINTAINER_EMAIL,
      url=URL,
      download_url=DOWNLOAD_URL,
      license=LICENSE,
      packages=['edl',],
      classifiers=[
          'Development Status :: 4 - Beta',
          'Environment :: Console',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: GPL License',
          'Natural Language :: English',
          'Programming Language :: Python :: 2.7'],
      )
