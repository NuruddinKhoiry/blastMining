#!/usr/bin/env python3

from setuptools import setup, find_packages

with open('README.md') as f:
    long_description = f.read()

requirements = ["numpy >= 1.22.3",
	"pandas>=1.4.2",
	"fastnumbers>=3.2.1"]

setup(
    name='blastMining',
    version = '0.1.0',
    author = 'Ahmad Nuruddin Khoiri',
    author_email = 'nuruddinkhoiri34@gmail.com',
    long_description = long_description,
    license = 'GPLv3',
    packages=find_packages(),
    entry_points = {'console_scripts': ['blastMining = blastMining.blastMining:main']},
    scripts=['blastMining/script/blastMining_lca.sh', 'blastMining/script/blastMining_lca2.sh'],
    url='https://github.com/NuruddinKhoiry/blastMining.git',
    python_requires='>=3',
    install_requires=requirements,
    #install_requires = [# No dependencies listed here since we need to rely on conda anyway],
    include_package_data=True, 
    classifiers = [
        "Development Status :: 1 - Beta",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GPLv3 License",
        "Natural Language :: English",
        "Programming Language :: Python :: 3.6",
        "Topic :: Scientific/Engineering :: Bio-Informatics"
    ]
)


