from setuptools import setup, find_packages
from codecs import open

all_requirements = []
with open('setup-req.txt', encoding='utf-8') as f:
    all_requirements = f.read().splitlines()

setup(
    name='pyprot',

    version='0.1-alpha',

    description='Python Proteins',
    long_description="A package designed to represent and maniupate amino acid sequences.",

    url='https://github.com/StanIsAdmin/PyProt',
	download_url='https://github.com/StanIsAdmin/PyProt/archive/0.1-alpha.tar.gz',

    author='Stanislas Gueniffey',
    author_email='gueniffeystanislas@gmail.com',

    license='GPL-3.0',

    classifiers=[
        'Development Status :: 3 - Alpha',

        'Intended Audience :: Developers',
        'Intended Audience :: Education',
        'Intended Audience :: Science/Research',

        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        
        'Natural Language :: English',

        'Programming Language :: Python :: 3 :: Only',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],

    keywords='biology proteins sequence-alignment',

    packages=find_packages(exclude=[]),

    install_requires=all_requirements,

    extras_require={
        'test': ['nose']
    },
)
