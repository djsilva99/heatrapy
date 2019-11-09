---
layout: home
permalink: /docs/documentation/
author_profile: true
sidebar:
  nav: "docs"
---

This module allows to create thermal objects, establish thermal contact between them, activate or deactivate part of the materials, and compute the respective heat transfer processes, in 1D. The computation uses the finite difference method. At the momento it includes two different systems for caloric devices: a fully solid state device and an hydraulic active regenerative device. The computation is performed by executing the respective function, which is described in the System section.

The package is based on three classes that create generic models:

**object**

This class only creates a single thermal object. It includes two methods: material activation and material deactivation, of part of the object.

**system_objects**

This class creates a system of objects that can be in contact to each other and computes the respective heat transfer processes. It uses the class object for the creation of each thermal object.

**single_object**

This class computes the heat transfer processes involved in only one thermal object. It uses the class object for activating and deactivating the material.

At this moment there are two types of caloric systems (functions) that can be computed:

**solid_active_regenerator**

This function creates and computes an active regenerative system used for
refrigeration and heat pumps. The heat exchanger is the **solid material** itself. It can be used to compute caloric systems, e.g. magnetocaloric, electrocaloric, elastocaloric, and barocaloric.

**fluid_active_regenerator**

This function creates and computes an active regenerative system used for
magnetocaloric refrigeration and heat pumps. The heat exchanger is a **fluid**. It can be used to compute caloric systems, e.g. magnetocaloric, electrocaloric, elastocaloric, and barocaloric.

