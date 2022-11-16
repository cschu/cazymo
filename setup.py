# coding: utf-8
from setuptools import setup, find_packages
from codecs import open
from os import path
import sys

from cazymo import __version__ as cazymo_version
here = path.abspath(path.dirname("__file__"))

with open(path.join(here, "DESCRIPTION.md"), encoding="utf-8") as description:
	description = long_description = description.read()

	name="cazymo"
	version = cazymo_version

	if sys.version_info.major != 3:
		raise EnvironmentError("""{toolname} requires python3, and is not compatible with python2.""".format(toolname=name))

	setup(
		name=name,
		version=version,
		description=description,
		long_description=long_description,
		url="https://github.com/cschu/cazymo",
		author="Christian Schudoma",
		author_email="christian.schudoma@embl.de",
		license="MIT",
		classifiers=[
			"Development Status :: 4 - Beta",
			"Topic :: Scientific Engineering :: Bio/Informatics",
			"License :: OSI Approved :: MIT License",
			"Operating System :: POSIX :: Linux",
			"Programming Language :: Python :: 3.7",
			"Programming Language :: Python :: 3.8",
			"Programming Language :: Python :: 3.9",
			"Programming Language :: Python :: 3.10",
		],
		zip_safe=False,
		keywords="functional metagenomics profiling tool",
		packages=find_packages(exclude=["test"]),
		install_requires=[
			"intervaltree",
			"numpy",
			"pandas",
			"sqlalchemy",
			"pysam",
		],
		entry_points={
			"console_scripts": [
				"cazymo=cazymo.__main__:main",
				"collate_counts=gq.bin.collate_counts:main",
				"build_domain_database=gq.bin.build_domain_database:main",
				"build_bed_database=gq.bin.build_bed_database:main",
			],
		},
		package_data={},
		include_package_data=True,
		data_files=[],
	)
