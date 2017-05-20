from setuptools import setup, find_packages
from codecs import open

all_requirements = []
with open('setup-req.txt', encoding='utf-8') as f:
    all_requirements = f.read().splitlines()

setup(
    name='pyprot',

    version='0.1',

    description='Python Proteins',
    long_description="A package designed to represent and maniupate amino acid sequences.",

    url='https://github.com/StanIsAdmin/PyProt',

    author='Stanislas Gueniffey',
    author_email='',

    license='GPL-3.0',

    classifiers=[
        'Development Status :: 3 - Alpha',

        'Intended Audience :: Developers',
        'Intended Audience :: Biologists',
        'Topic :: Biology :: Proteins',

        'License :: OSI Approved :: MIT License',

        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
    ],

    keywords='biology proteins sequence-alignment',

    packages=find_packages(exclude=[]),

    install_requires=all_requirements,

    extras_require={
        'test': ['nose']
    },
)
