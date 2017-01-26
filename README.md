# NGSIM.jl

A Julia package for working with the Next Generation Simulation dataset (NGSIM).
Was tested on the [Highway 101](http://www.fhwa.dot.gov/publications/research/operations/07030/) and [I-80](http://www.fhwa.dot.gov/publications/research/operations/06137/) datasets.

[![Build Status](https://travis-ci.org/tawheeler/NGSIM.jl.svg?branch=master)](https://travis-ci.org/tawheeler/NGSIM.jl)
[![Coverage Status](https://coveralls.io/repos/tawheeler/NGSIM.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/tawheeler/NGSIM.jl?branch=master)

This package is fully compatible with [AutomotiveDrivingModels.jl](https://github.com/tawheeler/AutomotiveDrivingModels.jl), providing the `Roadway` and `Trajdata` types from the NGSIM data. Roadway geometry was extracted from the NGSIM CAD files. The vehicle trajectories were filtered to provide better global positions and orientation.

You will currently have to supply the NGSIM trajectory data yourself, until I can find a way to host it on GitHub (and get permission from NGSIM to do so).

## Git It

You just clone it! Note that you also have to clone my [Vec](https://github.com/tawheeler/Vec.jl) package.

```julia
Pkg.clone("https://github.com/tawheeler/Vec.jl.git")
Pkg.clone("https://github.com/tawheeler/AutomotiveDrivingModels.jl.git")
Pkg.clone("https://github.com/tawheeler/NGSIM.jl.git")
```
