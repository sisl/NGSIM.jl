# NGSIM.jl
A Julia package for working with the Next Generation Simulation dataset (NGSIM).
I have tested this on the [Highway 101](http://www.fhwa.dot.gov/publications/research/operations/07030/) and [I-80](http://www.fhwa.dot.gov/publications/research/operations/06137/) datasets.

I provide two primary types: `Roadway` and `Trajdata`. The first represents the roadway at which the data was collected. Information was extracted from the NGSIM CAD files. The `Trajdata` type contains the vehicle positions over time. The raw data is availabe in `TrajdataRaw`, or smoothed and conveniently extracted into types in `Trajdata`. A multitude of accessor and utility functions are provided.

You will currently have to supply the NGSIM trajectory data yourself, until I can find a way to host it on GitHub (and get permission from NGSIM to do so).

## Git It

You just clone it! Note that you also have to clone my [Vec](https://github.com/tawheeler/Vec.jl) package.

```julia
Pkg.clone("https://github.com/tawheeler/Vec.jl.git")
Pkg.clone("https://github.com/tawheeler/NGSIM.jl.git")
```
