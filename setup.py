from setuptools import setup, find_packages
from codecs import open
from os import path
from version import __version__

here = path.abspath(path.dirname(__file__))

# Get the long description from the README and CHANGELOG files
with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
    readme = f.read()
with open(path.join(here, 'CHANGELOG.rst'), encoding='utf-8') as f:
    changelog = f.read()
long_description = readme + '\n\n' + changelog

setup(
    name='FEV_KEGG',  # Required

    version=__version__,  # Required

    description='FEV@KEGG allows for easy analysis of metabolic networks of organisms in KEGG.',  # Required

    long_description=long_description,  # Optional
    long_description_content_type="text/x-rst",

    url='https://github.com/ryhaberecht/FEV-KEGG',  # Optional

    author='Robin-Yves Haberecht',  # Optional

    # https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[  # Optional
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 4 - Beta',

        'Intended Audience :: Science/Research',
        'Intended Audience :: Education',
        'Topic :: Scientific/Engineering :: Bio-Informatics',

        'License :: OSI Approved :: MIT License',

        'Programming Language :: Python :: 3 :: Only',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        
        'Environment :: Console',
    ],

    keywords='bioinformatics metabolism graphs KEGG',  # Optional

    packages=find_packages(exclude=['contrib', 'docs', 'tests']),  # Required

    install_requires=['networkx', # Optional
                      'anytree',
                      'jsonpickle',
                      'tqdm',
                      'beautifulsoup4',
                      'retrying',
                      'appdirs',
                      'yattag'],
    
    extras_require={  # Optional
        'python34': ['typing'], # only required in Python == 3.4
        'draw_image': ['pygraphviz'],
        'draw_window': ['matplotlib'],
        },
    
    python_requires='~=3.4', # Python >=3.4, but not 4.x

    project_urls={  # Optional
        'Bug Reports': 'https://github.com/ryhaberecht/FEV-KEGG/issues',
        'Source': 'https://github.com/ryhaberecht/FEV-KEGG/',
    },
)