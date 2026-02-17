Material database
=================

Builtin database
---------

A built-in database of materials is provided in the database folder, although
it is currently limited. Custom material data can be used by specifying the
absolute path through the ``materials_path`` variable. For example,
``materials_path='<full_path>'`` and ``material='material_a'``. In this case, a
folder named ``material_a`` must exist in the directory ``<full_path>``,
containing all the required material data.

The built-in database includes several standard materials and can be consulted
`here
<https://github.com/djsilva99/heatrapy/tree/master/heatrapy/database>`_. Some
folder names end with _mag, indicating that the state of each point defines
whether a magnetic field is applied.

Contributions are welcome — feel free to add new materials to this database by
opening a Pull Request in the `github repository
<https://github.com/djsilva99/heatrapy>`_.

Material Data
---------

Each material folder should have 10 files:

- ``cp0.txt``: Specific heat capacity under no field (material inactive)
- ``cpa.txt``: Specific heat capacity under field (material active)
- ``k0.txt``: Thermal conductivity under no field (material inactive)
- ``ka.txt``: Thermal conductivity under field (material active)
- ``rho0.txt``: Density under no field (material inactive)
- ``rhoa.txt``: Density under field (material active)
- ``tadi.txt``: Adiabatic temperature change when applying a field (material
  activation)
- ``tadd.txt``: Adiabatic temperature change when removing a field (material
  inactivation)
- ``lheat0.txt``: Latent heat under no field (material inactive)
- ``lheata.txt``: Latent heat under field (material active)

Currently, all files have a .txt extension but are actually TSV files. Each
file contains two columns (without a header): the first column corresponds to
the temperature, and the second column corresponds to the physical property
values. Parameters will be interpolated between datapoints and will stay
constant if above or below datapoint range.

The used units are the following:

- cp (specific heat capacity): :math:`\mathrm{J\,kg^{-1}\,K^{-1}}`
- k (thermal conductivity): :math:`\mathrm{W\,m^{-1}\,K^{-1}}`
- lheat (latent heat): :math:`\mathrm{J\,kg^{-1}}`
- rho (density): :math:`\mathrm{kg\,m^{-3}}`
- tad(d/i) (adiabatic temperature change): :math:`\mathrm{K}`
