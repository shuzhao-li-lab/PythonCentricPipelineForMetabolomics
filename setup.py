from setuptools import setup, find_packages

with open("pcpfm/__init__.py") as f:
    exec([x for x in f.readlines() if '__version__' in x][0])

with open("README.md", "r") as f:
    long_description = f.read()

with open("requirements.txt", "r") as f:
    requirements = f.read()

setup(
  name='pcpfm',
  version=__version__,
  author='Joshua Mitchell',
  author_email='joshua.mitchell@jax.org',
  description='LC-MS metabolomics data processing pipeline',
  long_description=long_description,
  long_description_content_type="text/markdown",
  url='https://github.com/jmmitc06/PythonCentricPipelineForMetabolomics',
  license='MIT',
  keywords='metabolomics bioinformatics mass spectrometry',

  classifiers=[
    'Intended Audience :: Developers',
    'Intended Audience :: Science/Research',
    'License :: OSI Approved :: BSD License',
    'Natural Language :: English',
    'Operating System :: OS Independent',
    'Programming Language :: Python :: 3',
    'Topic :: Scientific/Engineering :: Bio-Informatics',
    'Topic :: Scientific/Engineering :: Chemistry',
    'Topic :: Software Development :: Libraries :: Python Modules',
  ],

  packages=find_packages(
    include=['*', '']
  ),
  include_package_data=True,
  #zip_safe=True,
  entry_points = {
        'console_scripts': ['pcpfm=pcpfm.main:CLI'],
    },

  python_requires='>=3.7',
  install_requires=requirements,
)