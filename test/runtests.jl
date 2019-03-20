using Test
using AutoViz
using NGSIM
using NBInclude

@nbinclude(joinpath(dirname(pathof(NGSIM)), "..", "jnotebooks", "Demo.ipynb"))

# the following test requires downloading the dataset
# td1 = load_trajdata(1)
# scene = get!(Scene(500), td1, 1000)
# render(scene, ROADWAY_101)