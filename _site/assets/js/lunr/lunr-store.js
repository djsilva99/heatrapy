var store = [{
        "title": "Caloric hydraulic system",
        "excerpt":"To compute a hydraulic active caloric regenerative system the function fluid_active_regenerator must be called. The active regenerative processes can be used with the several allowed modes for application and removal of fields. Cascades of materials can also be computed. To run one simulation type &gt;&gt;&gt; ht.fluid_active_regenerator(file_name, amb_temperature=298, fluid_length=160, ... MCM_length=50,...","categories": [],
        "tags": [],
        "url": "http://localhost:4000/docs/caloric_hydraulic_system/",
        "teaser":null},{
        "title": "Caloric solid state system",
        "excerpt":"To compute a fully solid state caloric system the function solid_active_regenerator must be called. The active regenerative processes can be used with the several allowed modes for application and removal of fields. Cascades of materials can also be computed. To run one simulation type &gt;&gt;&gt; ht.solid_active_regenerator(file_name, amb_temperature=293, ... left_thermalswitch_length=2, ......","categories": [],
        "tags": [],
        "url": "http://localhost:4000/docs/caloric_solid_state_system/",
        "teaser":null},{
        "title": "Documentation",
        "excerpt":"This module allows to create thermal objects, establish thermal contact between them, activate or deactivate part of the materials, and compute the respective heat transfer processes, in 1D. The computation uses the finite difference method. At the momento it includes two different systems for caloric devices: a fully solid state...","categories": [],
        "tags": [],
        "url": "http://localhost:4000/docs/documentation/",
        "teaser":null},{
        "title": "Material data",
        "excerpt":" ","categories": [],
        "tags": [],
        "url": "http://localhost:4000/docs/material_data/",
        "teaser":null},{
        "title": "Material state",
        "excerpt":"Following the last subsection, to activate a piece of foo type:   foo.activate(initial_point, final_point)   This command will activate the material from the initial_point to the final_point.   To deactivate a piece of material type:   foo.deactivate(initial_point, final_point)   This command will deactivate the material from the initial_point to the final_point.   ","categories": [],
        "tags": [],
        "url": "http://localhost:4000/docs/material_state/",
        "teaser":null},{
        "title": "Single thermal object",
        "excerpt":"The single_object class solves numerically the heat conduction equation for a single thermal object, in which domains can have different materials. This class inherits the methods of the object class. To create a single thermal object type &gt;&gt;&gt; foo = ht.single_object(amb_temperature, materials=('Cu',), borders=(1, 11), ... materials_order=(0,), dx=0.01, dt=0.1, file_name='data.txt', ......","categories": [],
        "tags": [],
        "url": "http://localhost:4000/docs/single_thermal_object/",
        "teaser":null},{
        "title": "System of thermal objects",
        "excerpt":"The system_objects class can be used to compute heat transfer processes between solids. It creates a system of thermal objects, establishes contact between them and computes the respective thermal processes. To create a system of thermal objects type: &gt;&gt;&gt; foo = ht.system_objects(number_objects=2, materials=('Cu', 'Cu'), ... objects_length=(10, 10), amb_temperature=293, dx=0.01, ......","categories": [],
        "tags": [],
        "url": "http://localhost:4000/docs/system_of_thermal_objects/",
        "teaser":null},{
        "title": "Test",
        "excerpt":"Content   Text here  ","categories": [],
        "tags": [],
        "url": "http://localhost:4000/docs/test/",
        "teaser":null},{
        "title": "Thermal_object_class",
        "excerpt":"1D thermal objects incorporates all thermal physical properties (temperature, specific heat, thermal conductivity and density) and boundary conditions, for a discretized length. Thermal objects can include heat sources and can be in contact with other thermal objects. There are three classes for thermal objects. The object class is responsible for...","categories": [],
        "tags": [],
        "url": "http://localhost:4000/docs/thermal_object_class/",
        "teaser":null},]
