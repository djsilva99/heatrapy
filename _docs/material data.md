---
layout: home
permalink: /docs/material_data/
author_profile: true
sidebar:
  nav: "docs"
---

When creating thermal objects the materials data must be given as the ```materials``` argument. For each specific material a folder must be created. This folder must contain 6 files: k0.txt, ka.txt, cp0.txt, cpa.txt, rho0.txt, rhoa.txt, lheat0.txt and lheata.txt. With the excepton of lheat0.txt and lheata.txt, all the files includes (temperature, property value) pairs. lheat0.txt and lheata.txt lists all the (temperature, latent heat values) pairs. The 0 and a at the end of the filenames correspond to field 0 and field applied respectively. Finally the ```materials``` argument can be given as a built-in material, present in the database folder of the heatrapy root folder, or by specifying the path of the materials folder with argument ```materials_path```.

