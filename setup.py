long_description = "Mandos provides tools for measuring the fit of phylogenies to the stratographic record and hypothesis-testing of direct ancestorship."

from setuptools import setup, find_packages
from setuptools.extension import Extension
from Cython.Build import cythonize
import numpy

#def readme():
#    with open('README.rst') as f:
#        return f.read()

extensions = [
    Extension(
        "mandos.node2",
        ["mandos/node2.pyx"],
        #include_dirs=['/some/path/to/include/'], # not needed for fftw unless it is installed in an unusual place
        #libraries=['fftw3', 'fftw3f', 'fftw3l', 'fftw3_threads', 'fftw3f_threads', 'fftw3l_threads'],
        #library_dirs=['/some/path/to/include/'], # not needed for fftw unless it is installed in an unusual place
    ),
    Extension(
        "mandos.stratoML",
        ["mandos/stratoML.pyx"],
        include_dirs=[numpy.get_include()],
        #include_dirs=['/some/path/to/include/'], # not needed for fftw unless it is installed in an unusual place
        #libraries=['fftw3', 'fftw3f', 'fftw3l', 'fftw3_threads', 'fftw3f_threads', 'fftw3l_threads'],
        #library_dirs=['/some/path/to/include/'], # not needed for fftw unless it is installed in an unusual place
    ),
    Extension(
        "mandos.calc_bm_likelihood",
        ["mandos/calc_bm_likelihood.pyx"],
        include_dirs=[numpy.get_include()],
    ),
    Extension(
        "mandos.tree_likelihood_calculator",
        ["mandos/tree_likelihood_calculator.pyx"],
        include_dirs=[numpy.get_include()],
    ),
]


setup(name='mandos',
      version='0.1',
      description='Phylogenetic analysis of stratigraphic occurrences and morphology',
      long_description=long_description,
      classifiers=[
        'Development Status :: 3 - Alpha',
        'License :: OSI Approved :: GPL License',
        'Programming Language :: Python :: 2.7'
      ],
      keywords='',
      url='http://github.com/carolinetomo/mandos',
      author='Caroline Parins-Fukuchi',
      author_email='cfukuchi@umich.edu',
      license='GPL',
      packages=find_packages(exclude = ['scripts','tests']),#['mandos'],
      ext_modules=cythonize(extensions), 
      #entry_points = {'console_scripts':['search-mandos-trees = mandos.command_line:main'],},
      #package_data = `
      #include_package_data= True,
      install_requires=[
          'scipy',
      ],
      zip_safe=False)


