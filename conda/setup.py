from setuptools import setup

setup(
    name='libsmoother',
    version='0.1.0',
    description='A example Python package',
    url='https://github.com/MarkusRainerSchmidt/libSmoother',
    author='Markus Schmidt',
    author_email='markus.rainer.schmidt@gmail.com',
    license='MIT',
    packages=['libsmoother'],
    install_requires=[
        # @todo currently done via conda environment
    ],

    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'License :: MIT License',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
    ],
)