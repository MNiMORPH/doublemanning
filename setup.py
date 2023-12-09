import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="doublemanning",
    version="1.0.0-alpha",
    author="Andrew D. Wickert",
    author_email="awickert@umn.edu",
    description="Stage-discharge relationships: double-Manning approach",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/MNiMORPH/rating-curve-2x-manning/",
    install_requires=[
        'pandas',
        'numpy',
        'scipy',
        'matplotlib',
        'scikit-learn',
        ],
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Hydrology",
        "Intended Audience :: Science/Research",
    ],
    entry_points = {
        'console_scripts': [
            'doublemanning-fit = doublemanning.cli:doubleMmanningFit',
            'doublemanning-calc = doublemanning.cli:doubleManningCalc'
        ],
    },
    keywords='fluvial geomorphology sediment transport landscape evolution',
    project_urls={
        'Model page': 'https://csdms.colorado.edu/wiki/Model:GRLP',
    },
)
