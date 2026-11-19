import glob
import edl
from setuptools import setup

DESCRIPTION = "Python scripts for metagenomics"
with open('README.md') as readme_handle:
    LONG_DESCRIPTION = readme_handle.read()
NAME = "py-metagenomics"
AUTHOR = "John Eppley"
AUTHOR_EMAIL = "jmeppley@mit.edu"
MAINTAINER = "John Eppley"
MAINTAINER_EMAIL = "jmeppley@mit.edu"
URL = 'http://eddelong.mit.edu/'
DOWNLOAD_URL = 'http://github.com/jmeppley/py-metagenomics'
LICENSE = 'GPL'
VERSION = edl.__version__

setup(name=NAME,
      version=VERSION,
      description=DESCRIPTION,
      long_description=LONG_DESCRIPTION,
      long_description_content_type='text/markdown',
      author=AUTHOR,
      author_email=AUTHOR_EMAIL,
      maintainer=MAINTAINER,
      maintainer_email=MAINTAINER_EMAIL,
      url=URL,
      download_url=DOWNLOAD_URL,
      license=LICENSE,
      packages=['edl', ],
      python_requires='>=3.10',
      install_requires=['numpy', 'biopython>=1.80', 'pandas'],
      extras_require={'plots': ['matplotlib']},
      scripts=[s for s in glob.glob('*.py')
               if s != 'setup.py'],
      classifiers=[
          'Development Status :: 4 - Beta',
          'Environment :: Console',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: GNU General Public License v2 (GPLv2)',
          'Natural Language :: English',
          'Programming Language :: Python :: 3',
          'Programming Language :: Python :: 3 :: Only'],
      )
