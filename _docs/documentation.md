---
layout: home
permalink: /docs/documentation/
author_profile: true
sidebar:
  nav: "docs"
---

The heatrapy module is based on thermal objects that can be created, activated, deactivated, can establish contact with other thermal objects and can be computed (its heat transfer processes). Navigate with the left sidebar for the full documentation. As a brief generic information, the heatrapy module is based on three classes that create generic models:

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

