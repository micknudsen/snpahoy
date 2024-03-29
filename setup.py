from setuptools import setup, find_packages

setup(

    name='snpahoy',
    version='0.5.3',

    packages=find_packages('src'),
    package_dir={'': 'src'},

    test_suite='tests',

    entry_points={
        'console_scripts': ['snpahoy = snpahoy.client:run']
    },

    python_requires='>=3.7',

    install_requires=[
        'click',
        'pysam >=0.15',
    ],

    author='Michael Knudsen',
    author_email='micknudsen@gmail.com',
    license='MIT'

)
