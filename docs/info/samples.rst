Samples
=======

As an example, the code below computes 100 s of a Cu 2D solid with a
water circle and a aluminum square.

.. code-block:: python

   >>> import heatrapy as htp
   >>> example = htp.SingleObject2D(
   ...     293,
   ...     material='Cu',
   ...     boundaries=(300, 0, 0, 0),
   ...     size=(20, 20)
   ... )
   >>> example.change_material(
   ...     material='Al',
   ...     shape='square',
   ...     initial_point=(5, 2),
   ...     length=(4, 4)
   ... )
   >>> example.change_material(
   ...     material='water',
   ...     shape='circle',
   ...     initial_point=(5, 13),
   ...     length=4
   ... )
   >>> example.compute(100, 10)

The output at the end of the computation is the following:

.. image:: https://raw.githubusercontent.com/djsilva99/heatrapy/refs/heads/master/img/example.png

