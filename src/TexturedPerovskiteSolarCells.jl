module TexturedPerovskiteSolarCells

using ChargeTransport # for the units in the parameter file
using ExtendableGrids # for generating the grids
using NPZ             # for reading in npy files (optical photogeneration rate)
using Interpolations  # for interpolating the photogeneration input onto the FVM mesh
using Triangulate  # for visualizing postprocess data on non-uniform mesh
using SimplexGridFactory # for visualizing postprocess data on non-uniform mesh

# for internal data handling (naming "inspired" by DrWatson)
# we do not export these functions
datadir(args...) = joinpath(pkgdir(TexturedPerovskiteSolarCells), "data", args...)
scriptsdir(args...) = joinpath(pkgdir(TexturedPerovskiteSolarCells), "scripts", args...)


include("photogeneration_reader.jl")
export MaxwellPhotogeneration


include("grid/generate_grid.jl")
export generate_grid


end # module TexturedPerovskiteSolarCells
