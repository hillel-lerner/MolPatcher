from setuptools import setup, find_packages

setup(
    name="molpatcher",
    version="1.0",
    packages=find_packages(),
    entry_points={
        "console_scripts": [
            "molpatcher = mol_patcher.main:main" 
        ]
    }
)