# MarkovJunior.jl

A Julia reimagining of [this awesome procedural generation algorithm](https://github.com/mxgmn/MarkovJunior/).

````julia
] add MarkovJunior
> MarkovJunior.main()
````

![A screenshot of the Terrain.jl scene](docs/Screenshot.png)

[![A demo video, sped up 2x](https://img.youtube.com/vi/mjBu7Omch1s/0.jpg)](https://youtu.be/mjBu7Omch1s)

*A YouTube video demo*

## Description

The released 0.1 version is a proof-of-concept, totally functional and fun to play with.
It is missing many of the features of the original library,
  but boasts support for any number of dimensions and comes with the interactive tool.

Math and rendering is all on top of my [B+ game framework](https://github.com/heyx3/B-plus).
The current (v0.1) tool is only written to show 2D scenes.

Currently on main I am building a comprehensive v2.0, including:
* New concise syntax
* Much more flexible rewrite rules, including randomization "Modifiers" and multidimensional patterns
* Continued support for any-dimensional grids (e.g. 4D for animated 3D)
* More operations from the original software
* The ability to run as both a standalone tool and a C-style DLL
* The ability to *extend* the syntax with your own:
  * Operations (like the built-in rewriting, upscaling/downscaling, WFC, etc)
  * Biases (like the original `field` and `observe`)
  * Priorities (deciding how rewrite operations balance different kinds of rules, like "earliest-first")
* A renderer for 3D grids (and 4D, as animated 3D) in the tool
* A way to save images and animations through the tool

## Development

During development this tool uses the master branch of B+ (and its sub-packages).
Clone those four packages (BplusCore, BplusApp, BplusTools, Bplus)
  and run the following from the Julia REPL:

````julia
# Add local B+ sub-packages to B+
] activate ../Bplus.jl
] dev ../BplusCore.jl ../BplusApp.jl ../BplusTools.jl
# Add local B+ to MarkovJunior
] activate .
] dev ../Bplus.jl
````

## Syntax

I've developed my own Domain-Specific Language ("DSL"), separate from the original MarkovJunior.

The final planned version [is documented here](docs/dsl.md).

The current version is undocumented but pretty simple to work out from the sample scenes.

## Scenes

Different algorithm setups can be found in the *scenes/* folder.
They are commented with explanations of how they work.