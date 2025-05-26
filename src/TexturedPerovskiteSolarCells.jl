module TexturedPerovskiteSolarCells

using Triangulate
using SimplexGridFactory
using Interpolations
using PyPlot
using ExtendableGrids
using ChargeTransport
using NPZ

# for internal data handling (naming "inspired" by DrWatson)
# we do not export these functions
datadir(args...) = joinpath(pkgdir(TexturedPerovskiteSolarCells), "data", args...)
scriptsdir(args...) = joinpath(pkgdir(TexturedPerovskiteSolarCells), "scripts", args...)

include("photogeneration_reader.jl")
export MaxwellPhotogeneration

include("grid/generate_grid.jl")
export generate_grid

end # module TexturedPerovskiteSolarCells
