from setuptools import setup, find_packages
import glob

setup(
	name="taxa-atlas",
	version="0.0.0.dev0",
	packages=find_packages(),
	##TODO add author and other details
	zip_safe=False,
	scripts=glob.glob("scripts/*.py"),
	install_requires=[
	'scipy',
	'networkx',
	'python-louvain'
	]
	)