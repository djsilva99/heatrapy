# -*- coding: utf-8 -*-
"""heatrapy: Library for simulating heat transfer processes

Computes heat transfer processes in solids for 1D and 2D models.

"""

from setuptools import setup

description = (
    'This package is a module for simulating dynamic 1D and 2D heat transfer '
    'processes by using the finite difference method. The packages is based '
    'on thermal objects. Two different approaches can be used: single thermal '
    'object and a system of thermal objects that can contact with each other. '
    'The heatrapy package includes the incorporation of phase transitions.'
)


setup(
    name='heatrapy',
    version='2.1.1',
    description='Library for simulating 1D and 2D heat transfer processes',
    long_description=description,
    classifiers=[
        'Development Status :: 3 - Alpha',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.14',
        'Topic :: Scientific/Engineering :: Physics',
    ],
    url='https://github.com/djsilva99/heatrapy',
    author='Daniel Silva',
    author_email='djsilva@gmx.com',
    license='MIT',
    packages=['heatrapy'],
    install_requires=[
        'numpy==2.4.2',
        'matplotlib==3.10.8'
    ],
    include_package_data=True,
    zip_safe=False
)
