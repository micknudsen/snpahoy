from setuptools import setup, find_packages

setup(

    name='snpahoy',
    version='0.0',

    packages=find_packages('src'),
    package_dir={'': 'src'},

    test_suite='tests',

    entry_points={
        'console_scripts': ['snpahoy = snpahoy.client:main']
    },

    python_requires='>=3.6',

    install_requires=[
        'pysam',
    ],

    author='Michael Knudsen',
    author_email='micknudsen@gmail.com',
    license='MIT'

)
