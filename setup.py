from setuptools import setup, find_packages

setup(
    name='vdm-convert',
    version='0.1',
    packages=find_packages(),
    install_requires=[
        'prody',
        'numpy',
        'pandas',
        'pyarrow',
        'fastparquet'
    ],
    entry_points={
        'console_scripts': [
            'vdm_convert=vdm_convert.parquet_convert:main',
            'vdm_display=vdm_convert.launch_pymol:main',
        ],
    },
    author='Nicholas Freitas',
    author_email='nicholas.freitas@ucsf.edu',
    description='For converting VDM parquet files to other formats',
    long_description='https://github.com/njf042/vdm-convert',
    url='',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
)