
include("grid_1D.jl")
include("grid_2D_planar.jl")
include("grid_2D_nanotextured.jl")
include("grid_2D_nanotextured_8p0eM7.jl")

 function generate_grid(;gridDim = 1, type = "nanotextured", #"planar"
                        amplitude = 0.8e-7, parameter_file = nothing, demo_run::Bool)

    if gridDim == 1
        grid = generate_grid1D(;parameter_file, demo_run = demo_run)
    elseif gridDim == 2
        if type == "nanotextured"
            if amplitude < 8.0e-7
                grid = generate_grid2D_nanotextured(amplitude = amplitude, parameter_file = parameter_file, demo_run = demo_run)
            else
                grid = generate_grid2D_nanotextured_8p0eM7(amplitude = amplitude, parameter_file = parameter_file, demo_run = demo_run)
            end
        elseif type == "planar"
            grid = generate_grid2D_planar(parameter_file = parameter_file, demo_run = demo_run)
        end

    end

    return grid

 end
