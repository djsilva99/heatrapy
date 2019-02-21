# -*- coding: utf-8 -*-
"""heatrapy: Library for simulating heat transfer processes

Computes heat transfer processes in solids for 1-dimensional models.
It includes two models for simulating caloric systems.

"""

from setuptools import setup

description = 'Computes heat transfer processes in solids for 1-dimensional '
description += 'models. It includes two models for simulating caloric '
description += 'systems.'


setup(name='heatrapy',
      version='1.0.3',
      description='Library for simulating heat transfer processes',
      long_description=description,
      classifiers=[
          'Development Status :: 3 - Alpha',
          'License :: OSI Approved :: MIT License',
          'Programming Language :: Python :: 3.5',
          'Topic :: Scientific/Engineering :: Physics',
      ],
      url='https://github.com/djsilva99/heatrapy',
      author='Daniel Silva',
      author_email='djsilva99@gmail.com',
      license='MIT',
      packages=['heatrapy'],
      install_requires=[
          'numpy',
      ],
      include_package_data=True,
      zip_safe=False)
