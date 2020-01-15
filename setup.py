from setuptools import setup, find_packages
from os import path

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
#with open(path.join(here, 'README.md'), encoding='utf-8') as f:
#    long_description = f.read()

setup(
    name='PhotoManip',
    version='0.0.5',
    #description=''
    #long_description=long_description,
    #long_description_content_type='text/markdown',    
    author='Maayane T. Soumagnac',
    author_email='maayane.soumagnac@weizmann.ac.il',  # Optional
    keywords='astronomy',  # Optional
    packages=["Photomanip"],
    install_requires=['numpy','matplotlib','pandas','pylab','csv','os'],  # Optional
    python_requires='>=2.7.10',
)
