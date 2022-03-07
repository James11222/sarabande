"""A setuptools based setup module.
See:
https://packaging.python.org/guides/distributing-packages-using-setuptools/
https://github.com/pypa/sampleproject
"""

# Always prefer setuptools over distutils
from setuptools import setup
import pathlib

here = pathlib.Path(__file__).parent.resolve()

setup(
    name='sarabande',  # Required
    version='0.0.10',  # Required
    description='Tool for measuring 3/4 PCFs on discrete periodic data.',  # Optional
    long_description=(here / 'README.md').read_text(encoding='utf-8'),  # Optional
    long_description_content_type='text/markdown',  # Optional (see note above)
    url='https://github.com/James11222/sarabande/',  # Optional

    author='James Sunseri',  # Optional

    author_email='jamessunseri@berkeley.edu',  # Optional

    classifiers=[  # Optional
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 3 - Alpha',

        # Indicate who your project is intended for
        'Intended Audience :: Developers',
        'Topic :: Software Development :: Build Tools',
        # Pick your license as you wish
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        "Programming Language :: Python :: 3.10",
        'Programming Language :: Python :: 3 :: Only'],
    keywords='npcf, 3pcf, 4pcf, cosmology, hydrodynamics',  # Optional
    packages=['sarabande'],
    platforms=['any'],
    license="MIT",
    python_requires='>=3.6, <4',
    setup_requires=['pytest-runner'],
    # install_requires=['requests'],  # Optional
    tests_require=['pytest'],
    include_package_data=True,
    package_data={'': ['*.npy']},
    project_urls={  # Optional
        'Bug Reports': 'https://github.com/James11222/sarabande/issues',
        'Funding': 'https://donate.pypi.org',
        'Check out my Website!': 'http://www.jamessunseri.com',
        'Source': 'https://github.com/James11222/sarabande/',
    }
)
