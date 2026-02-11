# Julia MarkovJunior

![A screenshot of the Terrain.jl scene](docs/Screenshot.png)

[![A demo video, sped up 2x](https://img.youtube.com/vi/mjBu7Omch1s/0.jpg)](https://youtu.be/mjBu7Omch1s)

A YouTube video demo

## Description

A work-in-progress recreation of [this awesome procedural generation algorithm](https://github.com/mxgmn/MarkovJunior/).

It can work as both a Julia package (though not on the registry yet),
  and a standalone tool on top of the [B+ game framework](https://github.com/heyx3/B-plus).
The tool is currently written for 2D scenes, with 3D on the horizon,
  but the algorithm (and DSL) can operate in any number of dimensions.
Note that currently it uses the master branch of B+ (and its sub-packages),
  so they need to be cloned and added to this project with `dev`

The package has its own Domain-Specific Language ("DSL"), separate from the original MarkovJunior,
  which will be documented in detail asap.
