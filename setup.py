# -*- coding: utf-8 -*-
"""heatrapy: Library for simulating heat transfer processes

Computes heat transfer processes in solids for 1-dimensional models.
It includes two models for simulating caloric systems.

"""

from setuptools import setup

description = 'This package is a module for simulating dynamic 1D and 2D heat '
description += 'transfer processes by using the finite difference method. The '
description += 'packages is based on thermal objects. Two different '
description += 'approaches can be used: single thermal object and a system of '
description += 'thermal objects that can contact with each other. The heatrapy'
description += 'package includes the modeling of caloric effects and the '
description += 'incorporation of phase transitions.'


setup(name='heatrapy',
      version='2.0.3',
      description='Library for simulating 1D and 2D heat transfer processes',
      long_description=description,
      classifiers=[
          'Development Status :: 3 - Alpha',
          'License :: OSI Approved :: MIT License',
          'Programming Language :: Python :: 3.7',
          'Topic :: Scientific/Engineering :: Physics',
      ],
      url='https://github.com/djsilva99/heatrapy',
      author='Daniel Silva',
      author_email='djsilva@gmx.com',
      license='MIT',
      packages=['heatrapy'],
      install_requires=[
          'numpy>=1.20.1', 'matplotlib==3.3.4'
      ],
      include_package_data=True,
      zip_safe=False)
