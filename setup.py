import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="iedb_immrep",
    version="0.1",
    author="eve richardson",
    author_email="erichardson@lji.org",
    description="code to parse iedb data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/erichardson97/iedb_immrep25",
    packages=setuptools.find_packages(),
    package_data={
        'iedb_immrep': [
            'dat/*'
        ]
    },
    include_package_data=True,
    install_requires=[
        'pandas',
        'biopython',
        'numpy',
        'tidytcells',
        'requests',
        'pyarrow',
        'openpyxl',
        'stitchr',
        'IMGTgeneDL',
        'airr',
        'wget'
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: Unix",
    ])