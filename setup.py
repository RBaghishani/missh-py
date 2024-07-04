from setuptools import setup, find_packages

setup(
    name='missh-py',
    version='0.1.0',
    description='Python port of the MISSH application',
    author='Your Name',
    author_email='your.email@example.com',
    packages=find_packages(where='src'),
    package_dir={'': 'src'},
    install_requires=[
        # Add your dependencies here
    ],
    entry_points={
        'console_scripts': [
            'missh-py=main:main',
        ],
    },
)
