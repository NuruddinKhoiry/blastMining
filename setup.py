#!/usr/bin/env python3

from setuptools import setup, find_packages

with open('README.md') as f:
    long_description = f.read()

requirements = ["numpy >= 1.22.3",
	"pandas>=1.4.2",
	"fastnumbers>=3.2.1"]

setup(
    name='blastMining',
    version = '0.1.1',
    author = 'Ahmad Nuruddin Khoiri',
    author_email = 'nuruddinkhoiri34@gmail.com',
    long_description = 'Mining NCBI BLAST output with vote and lca methods',
    license = 'GPLv3',
    packages=find_packages(),
    entry_points = {'console_scripts': ['blastMining = blastMining.blastMining:main']},
    scripts=['blastMining/script/blastMining_lca.sh', 'blastMining/script/blastMining_lca2.sh'],
    url='https://github.com/NuruddinKhoiry/blastMining.git',
    python_requires='>=3',
    install_requires=requirements,
)


