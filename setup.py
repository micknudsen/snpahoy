from setuptools import setup, find_packages

setup(

    name='snpahoy',
    version='0.1.1',

    packages=find_packages('src'),
    package_dir={'': 'src'},

    test_suite='tests',

    entry_points={
        'console_scripts': ['snpahoy = snpahoy.client:run']
    },

    python_requires='>=3.6',

    install_requires=[
        'click',
        'pysam',
    ],

    author='Michael Knudsen',
    author_email='micknudsen@gmail.com',
    license='MIT'

)
