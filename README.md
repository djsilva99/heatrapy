# heatrapy

<img src="https://github.com/danieljosesilva/heatrapy/blob/master/img/heatrapy.png" alt="Drawing" height="30"/> heatrapy v1.1.0

This package is a module for simulating dynamic heat transfer processes involving caloric effects in 1.5D systems by using the finite difference method. It is focused on heat conduction, and includes two subpackages for computing caloric systems. The python library <a href='https://github.com/djsilva99/physplotlib'>physplotlib</a> can be used for the visualization of the output data.

For full documentation visit the <a href='https://djsilva99.github.io/heatrapy'>heatrapy website</a>.

author: Daniel Silva (djsilva@gmx.com) <br> current version: v1.1.0

![heatrapy-example](https://github.com/djsilva99/heatrapy/blob/master/img/example.gif)


### Installation

To install heatrapy use the pip package manager:

```bash
$ pip install heatrapy
```

To import the heatrapy module type in the python shell:

```python
>>> import heatrapy as ht
```


### Contributing

Before making a pull request run the text scripts by going to the parent folder of the cloned repository and executing:
```bash
$ python -m heatrapy.test.unit_tests
$ python -m heatrapy.test.integration_tests
```
Note that these test will produce .txt files.

Besides coding, feel free to report bugs, suggest new features and new systems. A roadmap for future releases can be consulted in the wiki pages.
