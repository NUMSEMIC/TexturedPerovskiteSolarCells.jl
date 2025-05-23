

function MaxwellPhotogeneration(;gridDim = 1,
                                ########################
                                typeGrid  = "nanotextured", # "planar"
                                ########################
                                amplitude = 4.0e-7,
                                parameter_file, demo_run
                                )

    G_interpol = calc_interpol_photogen(typeGrid = typeGrid, amplitude = amplitude)

    grid       = generate_grid(gridDim = gridDim, type = typeGrid, amplitude = amplitude, parameter_file = parameter_file, demo_run = demo_run)

    G          = zeros(length(grid[Coordinates][1,:]))

    coord      = grid[Coordinates]
    subg       = subgrid(grid, [regionPero])
    iNode      = subg[NodeParents]

    if gridDim == 2

        for inode in iNode
            x = coord[1, inode]
            y = coord[2, inode]

            G[inode] = G_interpol(x, y)
        end

    elseif gridDim == 1

        for inode in iNode

            y = coord[inode]

            G[inode] = G_interpol(0.1e-7, y)
        end

    end

    return G

end


function calc_interpol_photogen(;typeGrid = "nanotextured", # "planar", #
                                amplitude = 4.0e-7
                                )
    ####################################################
	### read in JCMSuite generation profile
	####################################################

    helpampl = collect(string(amplitude));  helpampl[ findall(x -> x == '.', helpampl)[1] ] = 'p'
    textampl = join(helpampl)

    if typeGrid == "nanotextured"
        typeGridText = "nanotextured-ampl-$textampl"
    elseif typeGrid == "planar"
        typeGridText = "planar"
    end

    X = npzread("data/optics/$(typeGrid)/X_Coords-$(typeGridText).npy")
    Y = npzread("data/optics/$(typeGrid)/Y_Coords-$(typeGridText).npy") .- 1.0e-7
    G = npzread("data/optics/$(typeGrid)/Photogeneration-$(typeGridText).npy")
    G = real.(G[:, :, 1])

    G_interpol = Interpolations.interpolate((X, Y), G, Gridded(Linear()))

    # XXPlot = X' .* ones(length(Y))
    # YYPlot = Y' .* ones(length(X))

    # G_interpol_eval = G_interpol.(XXPlot', YYPlot)

    # figure()
    # contourf(XXPlot'./nm, YYPlot./nm, G_interpol_eval, levels = 40)#,
    # colorbar()

    return G_interpol

end


function test(;gridDim = 2, typeGrid = "planar", # "nanotextured", #
              #######################
              amplitude = 4.0e-7,
              parameter_file, demo_run
            )

    PyPlot.rc("font", family="serif", size=14)
    PyPlot.rc("mathtext", fontset="dejavuserif")
    PyPlot.close("all")

    grid = generate_grid(gridDim = gridDim, type = typeGrid, amplitude = amplitude, demo_run)
    G    = MaxwellPhotogeneration(gridDim = gridDim, typeGrid = typeGrid, amplitude = amplitude, parameter_file, demo_run)

    if gridDim == 1
        coord = grid[Coordinates]'

        plot(coord, G)
        PyPlot.grid()
        xlabel("position [nm]")
        ylabel("photogeneration [(m\$^3\$s)\$^{-1}\$]")
        PyPlot.tick_params(which ="both", labelsize=18)
        PyPlot.tight_layout()

    elseif gridDim == 2

        XX = grid[Coordinates][1, :]
        YY = grid[Coordinates][2, :]

        tricontourf(XX./nm, YY./nm, G, levels = 40)#, norm=matplotlib[:colors][:LogNorm](vmin=minimum(G_interpol_eval), vmax=maximum(G_interpol_eval)))

        colorbar(orientation = "horizontal", label = "photogeneration [(m\$^3\$s)\$^{-1}\$]")
        PyPlot.xlabel("width [nm]", fontsize=17)
        PyPlot.ylabel("length [nm]", fontsize=17)
        PyPlot.tick_params(which ="both", labelsize=18)
        PyPlot.tight_layout()

    end


end