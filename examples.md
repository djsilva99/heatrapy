---
layout: home
permalink: /examples/
author_profile: true
---
# Example

The following example computes a simple 1-dimensional model with 0.5 m of gadolinium (Gd). The system is initial at 293 K. One end of the system is at a fixed temperature of 300 K, while the other end is insulated. The used time step is 1 second, and the used space step is 0.05 m, so that the overall number of space points is 10. The system is initially deactivated. To create the model we initialize the object `example`:

```python
>>> example = ht.single_object(
...        293, materials=['Gd'], borders=[1,11],
...	   materials_order=[0], 293, materials=['Gd'],
...	   borders=[1,11], materials_order=[0],
...	   dx=0.05, dt=1., file_name='example.txt',
...	   boundaries=[300,0], Q=[], Q0=[],
...	   initial_state=False
...)
```

Then we compute the system for 30000 s, write the output values every 300 s, using the 'implicit_k(x)' solver:

```python
>>> example.compute(30000, 300, solver='implicit_k(x)')
```

Afterwards we activate the whole system:

```python
>>> example.activate(1, 10)
```

and we compute the system for 30000 s one more time:

```python
>>> example.compute(30000, 300, solver='implicit_k(x)')
```

The output data is stored in file example.txt. To visualize the temperature as a function of time for the point index 3 we use the <a href='https://github.com/danieljosesilva/physplotlib'>physplotlib</a> library and write:

```python
>>> import physplotlib as pp
>>> example_visualization = pp.statplot()
>>> example_visualization.loadFile('example.txt')
>>> example_visualization.verticalPlot(
...        [0], [0], [[3]], y_title='temperature (K)'
...)
```

This code will output the following time-dependent temperature plot:
![heatrapy-example](/assets/example.png)
As expected the temperature will tend to the 300 K fixed temperature. The activation of the material at 30000 s makes the temperature to jump to ~ 302 K.

