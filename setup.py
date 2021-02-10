from setuptools import setup

setup(
    name='pyClassRich',
    version='0.0.1',
    description='Chemical class enrichment for mass spectrometry',
    license='BSD-3-Clause',
    packages=['pyclassrich'],
    package_data={'pyclassrich': ['data/ChemOnt_2_1.obo']},
    author='Ricardo R. da Silva',
    scripts=['bin/pyclassrich'],
    author_email='ridasilva@usp.br'
)
