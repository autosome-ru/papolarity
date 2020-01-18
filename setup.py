import setuptools
from io import open
import os

version_module = {}
dirname = os.path.dirname(__file__)
with open(os.path.join(dirname, "src/papolarity/version.py")) as fp:
    exec(fp.read(), version_module)
    __version__ = version_module['__version__']

with open(os.path.join(dirname, "README.md"), encoding='utf8') as fh:
    long_description = fh.read()

setuptools.setup(
    name="papolarity",
    version=__version__,
    author="Ilya Vorontsov",
    author_email="vorontsov.i.e@gmail.com",
    description="Papolarity is a tool to analyze polarity of transcriptomic alignments such as Ribo-seq and RNA-seq.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/autosome-ru/papolarity",
    license="Papolarity is licensed under WTFPL, but if you prefer more standard licenses, feel free to treat it as MIT license",
    packages=setuptools.find_packages('src'),
    package_dir={'': 'src'},
    include_package_data=True,
    keywords='bioinformatics NGS coverage alignment polarity Ribo-seq RNA-seq',
    classifiers=[
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Development Status :: 5 - Production/Stable",
        "Operating System :: OS Independent",
        'Programming Language :: Python :: 3.7',
    ],
    python_requires='>=3.7',
    install_requires=['pybedtools >= 0.8.0', 'numpy >= 1.8.0', 'sklearn', 'six', 'matplotlib', 'seaborn'],
    extras_require={
        'dev': ['pytest', 'pytest-benchmark', 'flake8', 'tox', 'wheel', 'twine', 'setuptools_scm'],
    },
    entry_points={
        'console_scripts': [
            'papolarity=papolarity.cli:main',
        ],
    },
    use_scm_version=False,
    setup_requires=['setuptools_scm'],
)
