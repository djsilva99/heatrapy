# heatrapy

<img src="https://github.com/danieljosesilva/heatrapy/blob/master/img/heatrapy.png" alt="Drawing" height="30"/> heatrapy v2.0.0

This package is a module for simulating dynamic 1D and 2D heat transfer processes by using the finite difference method. The packages is based on thermal objects. Two different approaches can be used: single thermal object and a system of thermal objects that can contact with each other. Unlike the majority of software packages, heatrapy includes both the modeling of caloric effects and the incorporation of phase transitions.

For full documentation visit the <a href='https://djsilva99.github.io/heatrapy'>heatrapy website</a>.

author: Daniel Silva (djsilva@gmx.com) <br> current version: v2.0.0

![heatrapy-example](https://github.com/djsilva99/heatrapy/blob/master/img/example.gif)


### Installation

To install heatrapy use the pip package manager:

```bash
$ pip install heatrapy
```

To import the heatrapy module type in the python shell:

```python
>>> import heatrapy as htp
```


### Contributing

Before making a pull request run the text script by going to the parent folder of the cloned repository and executing:
```bash
$ python -m heatrapy.test.unit_tests
```

Besides coding, feel free to report bugs, suggest new features and new systems. A roadmap for future releases can be consulted in the wiki pages.
