#=

Code for visualizing density ratio of holes and electrons

=#
module Fig5RatioElectronsAndHoles

using PyPlot
using DelimitedFiles
using ChargeTransport
using ExtendableGrids
using TexturedPerovskiteSolarCells
using PyCall

# Import Python's colors module
colors = pyimport("matplotlib.colors")

# for convenience
datadir = TexturedPerovskiteSolarCells.datadir
scriptsdir = TexturedPerovskiteSolarCells.scriptsdir

color_map = "seismic"

include(scriptsdir("SingleJunction.jl"))

## for generating the data in manuscript, the options are:
# vETL = 10,   vHTL = 10,  and then V = "inival" or V = "V-1p2"
# vETL = 2000, vHTL = 500, and then V = "inival" or V = "end" (V = 1.2 V)
function main(;
        vETL = 10,  # 2000  # in cm/s
        vHTL = 10,  # 500   # in cm/s
        V = "V-1p2", # "inival",
        printText = true,
        saveFig = false,
        parameter_set = ParamsSingleJunction,
    )

    # use the destructuring operator to extract all the necessary parameters
    (; regionPero, ipsi, Ca, iphin, iphip, Fcc, regionPero, T, zn, zp, Nn, Np, En, Ep, Ca) = parameter_set()

    @local_unitfactors nm

    (; q, k_B) = ChargeTransport.constants

    PyPlot.rc("font", family = "sans-serif", size = 14)
    PyPlot.rc("mathtext", fontset = "dejavusans")
    PyPlot.close("all")

    pathSol = "vETL=$(vETL)vHTL=$(vHTL)/"

    figsize = (7.2, 5.6)
    #################################################################
    ## grid
    #################################################################

    amplitude = 3.0e-7
    amplitude2 = 6.0e-7

    helpampl = collect(string(amplitude));  helpampl[findall(x -> x == '.', helpampl)[1]] = 'p'
    textampl = join(helpampl)

    helpampl2 = collect(string(amplitude2));  helpampl2[findall(x -> x == '.', helpampl2)[1]] = 'p'
    textampl2 = join(helpampl2)

    grid1, ctsys1 = SingleJunction.main(gridDim = 2, typeGrid = "planar", generation = true, plotPostProcess = true, printText = false, plotting = false)
    subg1 = subgrid(grid1, [regionPero])
    data1 = ctsys1.fvmsys.physics.data
    ########
    grid2, ctsys2 = SingleJunction.main(gridDim = 2, typeGrid = "nanotextured", amplitude = amplitude, generation = true, plotPostProcess = true, printText = false, plotting = false)
    subg2 = subgrid(grid2, [regionPero])
    data2 = ctsys2.fvmsys.physics.data
    ########
    grid3, ctsys3 = SingleJunction.main(gridDim = 2, typeGrid = "nanotextured", amplitude = amplitude2, generation = true, plotPostProcess = true, printText = false, plotting = false)
    subg3 = subgrid(grid3, [regionPero])
    data3 = ctsys3.fvmsys.physics.data

    coord1 = subg1[Coordinates]
    coord2 = subg2[Coordinates]
    coord3 = subg3[Coordinates]

    coord1 = coord1 ./ nm
    coord2 = coord2 ./ nm
    coord3 = coord3 ./ nm

    #########################################################################################################
    #########################################################################################################

    sol1 = readdlm(datadir("sol", "$pathSol/Sol-2D-forw-planar-generation-Maxwell-reco-all-$V.dat"))'
    sol2 = readdlm(datadir("sol", "$pathSol/Sol-2D-forw-nanotextured-ampl-$textampl-generation-Maxwell-reco-all-$V.dat"))'
    sol3 = readdlm(datadir("sol", "$pathSol/Sol-2D-forw-nanotextured-ampl-$textampl2-generation-Maxwell-reco-all-$V.dat"))'

    nn1 = Nn[regionPero] .* Fcc[iphin].(zn * (q * (view(sol1[iphin, :], subg1) .- view(sol1[ipsi, :], subg1)) .+ En[regionPero]) ./ (k_B * T))
    nn2 = Nn[regionPero] .* Fcc[iphin].(zn * (q * (view(sol2[iphin, :], subg2) .- view(sol2[ipsi, :], subg2)) .+ En[regionPero]) ./ (k_B * T))
    nn3 = Nn[regionPero] .* Fcc[iphin].(zn * (q * (view(sol3[iphin, :], subg3) .- view(sol3[ipsi, :], subg3)) .+ En[regionPero]) ./ (k_B * T))

    np1 = Np[regionPero] .* Fcc[iphip].(zp * (q * (view(sol1[iphip, :], subg1) .- view(sol1[ipsi, :], subg1)) .+ Ep[regionPero]) ./ (k_B * T))
    np2 = Np[regionPero] .* Fcc[iphip].(zp * (q * (view(sol2[iphip, :], subg2) .- view(sol2[ipsi, :], subg2)) .+ Ep[regionPero]) ./ (k_B * T))
    np3 = Np[regionPero] .* Fcc[iphip].(zp * (q * (view(sol3[iphip, :], subg3) .- view(sol3[ipsi, :], subg3)) .+ Ep[regionPero]) ./ (k_B * T))

    # println(" ")
    if printText
        @show minimum(np1 ./ nn1), maximum(np1 ./ nn1)
        @show minimum(np2 ./ nn2), maximum(np2 ./ nn2)
        @show minimum(np3 ./ nn3), maximum(np3 ./ nn3)
    end

    #################################################################
    ## Plotting
    #################################################################

    subg1[Coordinates] = subg1[Coordinates] ./ nm
    subg2[Coordinates] = subg2[Coordinates] ./ nm
    subg3[Coordinates] = subg3[Coordinates] ./ nm

    # THX at JF!!!
    # https://github.com/j-fu/GridVisualize.jl/blob/1f2b299a436b7750702ccca282fa14152d80ebf9/src/pyplot.jl#L86
    function tridata(grid::ExtendableGrid)
        coord = grid[Coordinates]
        cellnodes = Matrix(grid[CellNodes])
        return coord[1, :], coord[2, :], transpose(cellnodes .- 1)
    end

    vmin = 1 / 50.0e1; vmax = 50.0e1
    norm = matplotlib.colors.LogNorm(vmin = vmin, vmax = vmax)

    figure(figsize = figsize)
    tripcolor(tridata(subg1)..., vcat(np1 ./ nn1...), norm = norm, shading = "gouraud", cmap = color_map, rasterized = true)
    xlabel(" \$x\$ [nm]", fontsize = 17)
    ylabel(" \$y\$ [nm]", fontsize = 17)
    axis([-20, 770, 20, 800])
    title("Ratio np/nn -- V = $V")
    cbar = colorbar(orientation = "vertical", label = "\$ n_{\\mathrm{n}}/ n_{\\mathrm{p}} \$", extend = "both", spacing = "proportional")
    cbar.ax.set_yscale("log")
    tight_layout()

    if saveFig
        savefig(datadir("2D-planar-ratio-nn-np-generation-Maxwell-forw-$V.pdf"))
    end

    #####################
    figure(figsize = figsize)
    tripcolor(tridata(subg2)..., vcat(np2 ./ nn2...), norm = norm, shading = "gouraud", cmap = color_map, rasterized = true)
    xlabel("\$ x\$ [nm]", fontsize = 17)
    ylabel("\$ y\$ [nm]", fontsize = 17)
    axis([-20, 770, 20, 800])
    title("Ratio np/nn -- V = $V")
    cbar = colorbar(orientation = "vertical", label = "\$ n_{\\mathrm{n}}/ n_{\\mathrm{p}} \$", extend = "both", spacing = "proportional")
    cbar.ax.set_yscale("log")
    tight_layout()

    if saveFig
        savefig(datadir("2D-nanotextured-ampl-$textampl-ratio-nn-np-generation-Maxwell-forw-$V.pdf"))
    end

    #####################
    figure(figsize = figsize)
    tripcolor(tridata(subg3)..., vcat(np3 ./ nn3...), norm = norm, shading = "gouraud", cmap = color_map, rasterized = true)
    xlabel(" \$x\$ [nm]", fontsize = 17)
    ylabel(" \$y\$ [nm]", fontsize = 17)
    axis([-20, 770, 20, 800])
    title("Ratio np/nn -- V = $V")
    cbar = colorbar(orientation = "vertical", label = "\$ n_{\\mathrm{n}}/ n_{\\mathrm{p}} \$", extend = "both", spacing = "proportional")
    cbar.ax.set_yscale("log")
    tight_layout()

    if saveFig
        savefig(datadir("2D-nanotextured-ampl-$textampl2-ratio-nn-np-generation-Maxwell-forw-$V.pdf"))
    end

    return nothing

end # main

end # module
