from setuptools import setup, find_packages

###TODO: fill this out
setup(
    name='vdm-convert',
    version='0.1',
    packages=find_packages(),
    install_requires=[
    ],
    entry_points={
        'console_scripts': [
            'example_script=example_package.example_script:main',
        ],
    },
    author='Nicholas Freitas',
    author_email='nicholas.freitas@ucsf.edu',
    description='For converting VDM parquet files to other formats',
    long_description='',
    url='',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
)