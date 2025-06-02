#=

Code for visualizing density ratio

=#
module FigS8RatioElectronsAndHoles

using PyPlot
using DelimitedFiles
using ChargeTransport
using ExtendableGrids
using TexturedPerovskiteSolarCells

# for convenience
datadir = TexturedPerovskiteSolarCells.datadir
scriptsdir = TexturedPerovskiteSolarCells.scriptsdir

color_map = "seismic"

include(scriptsdir("SingleJunction.jl"))

function main(;scanrate  = 1000.0,   # "10p0" # "0p001"
              typeReco  = "all",    # "radiative"
              generationUniform = false,
              IVDirection = "forw", # "rev"
              V = "end", # "inival"
              saveFig = false,
              parameter_file = scriptsdir("params_single_junction.jl")
              )

    include(parameter_file)

    if generationUniform
        generation = "uniform"
    else
        generation = "Maxwell"
    end

    PyPlot.rc("font", family="sans-serif", size=14)
    PyPlot.rc("mathtext", fontset="dejavusans")
    PyPlot.close("all")

    helpSR  = collect(string(scanrate));  helpSR[ findall(x -> x == '.', helpSR)[1] ] = 'p'
    textSR  = join(helpSR)

    pathSol = "parameter-$paramsname/$textSR"

    #################################################################
    ## grid
    #################################################################

    amplitude  = 3.0e-7
    amplitude2 = 5.0e-7
    amplitude3 = 7.0e-7

    helpampl = collect(string(amplitude));  helpampl[ findall(x -> x == '.', helpampl)[1] ] = 'p'
    textampl = join(helpampl)

    helpampl2 = collect(string(amplitude2));  helpampl2[ findall(x -> x == '.', helpampl2)[1] ] = 'p'
    textampl2 = join(helpampl2)

    helpampl3 = collect(string(amplitude3));  helpampl3[ findall(x -> x == '.', helpampl3)[1] ] = 'p'
    textampl3 = join(helpampl3)

    grid1, ctsys1 = SingleJunction.main(gridDim = 2, typeGrid = "planar", generation = true, generationUniform = generationUniform, plotPostProcess = true, printText = false, plotting = false)
    subg1         = subgrid(grid1, [regionPero])
    data1         = ctsys1.fvmsys.physics.data
    ########
    grid2, ctsys2 = SingleJunction.main(gridDim = 2, typeGrid = "nanotextured", amplitude = amplitude, generation = true, generationUniform = generationUniform, plotPostProcess = true, printText = false, plotting = false)
    subg2         = subgrid(grid2, [regionPero])
    data2         = ctsys2.fvmsys.physics.data
    ########
    grid3, ctsys3 = SingleJunction.main(gridDim = 2, typeGrid = "nanotextured", amplitude = amplitude2, generation = true, generationUniform = generationUniform, plotPostProcess = true, printText = false, plotting = false)
    subg3         = subgrid(grid3, [regionPero])
    data3         = ctsys3.fvmsys.physics.data
    ########
    grid4, ctsys4 = SingleJunction.main(gridDim = 2, typeGrid = "nanotextured", amplitude = amplitude3, generation = true, generationUniform = generationUniform, plotPostProcess = true, printText = false, plotting = false)
    subg4         = subgrid(grid4, [regionPero])
    data4         = ctsys4.fvmsys.physics.data

    coord1  = subg1[Coordinates]; coord2  = subg2[Coordinates]
    coord3  = subg3[Coordinates]; coord4  = subg4[Coordinates]

    coord1  = coord1./nm; coord2  = coord2./nm
    coord3  = coord3./nm; coord4  = coord4./nm

    #########################################################################################################
    #########################################################################################################

    sol1 = readdlm(datadir("sol", "$pathSol/Sol-2D-$IVDirection-planar-generation-$generation-reco-$typeReco-$V.dat"))'
    sol2 = readdlm(datadir("sol", "$pathSol/Sol-2D-$IVDirection-nanotextured-ampl-$textampl-generation-$generation-reco-$typeReco-$V.dat"))'
    sol3 = readdlm(datadir("sol", "$pathSol/Sol-2D-$IVDirection-nanotextured-ampl-$textampl2-generation-$generation-reco-$typeReco-$V.dat"))'
    sol4 = readdlm(datadir("sol", "$pathSol/Sol-2D-$IVDirection-nanotextured-ampl-$textampl3-generation-$generation-reco-$typeReco-$V.dat"))'

    nn1  = Nn[regionPero] .* Fcc[iphin].(zn*( q*(view(sol1[iphin, :], subg1).-view(sol1[ipsi, :], subg1)) .+ En[regionPero])./(kB*T))
    nn2  = Nn[regionPero] .* Fcc[iphin].(zn*( q*(view(sol2[iphin, :], subg2).-view(sol2[ipsi, :], subg2)) .+ En[regionPero])./(kB*T))
    nn3  = Nn[regionPero] .* Fcc[iphin].(zn*( q*(view(sol3[iphin, :], subg3).-view(sol3[ipsi, :], subg3)) .+ En[regionPero])./(kB*T))
    nn4  = Nn[regionPero] .* Fcc[iphin].(zn*( q*(view(sol4[iphin, :], subg4).-view(sol4[ipsi, :], subg4)) .+ En[regionPero])./(kB*T))

    np1  = Np[regionPero] .* Fcc[iphip].(zp*( q*(view(sol1[iphip, :], subg1).-view(sol1[ipsi, :], subg1)) .+ Ep[regionPero])./(kB*T))
    np2  = Np[regionPero] .* Fcc[iphip].(zp*( q*(view(sol2[iphip, :], subg2).-view(sol2[ipsi, :], subg2)) .+ Ep[regionPero])./(kB*T))
    np3  = Np[regionPero] .* Fcc[iphip].(zp*( q*(view(sol3[iphip, :], subg3).-view(sol3[ipsi, :], subg3)) .+ Ep[regionPero])./(kB*T))
    np4  = Np[regionPero] .* Fcc[iphip].(zp*( q*(view(sol4[iphip, :], subg4).-view(sol4[ipsi, :], subg4)) .+ Ep[regionPero])./(kB*T))

    # println(" ")
    # @show minimum(nn1./np1), maximum(nn1./np1)
    # @show minimum(nn2./np2), maximum(nn2./np2)
    # @show minimum(nn3./np3), maximum(nn3./np3)
    # @show minimum(nn4./np4), maximum(nn4./np4)
    # @show minimum(nn5./np5), maximum(nn5./np5)

    #################################################################
    ## Plotting
    #################################################################

    subg1[Coordinates] = subg1[Coordinates]./nm; subg2[Coordinates] = subg2[Coordinates]./nm
    subg3[Coordinates] = subg3[Coordinates]./nm; subg4[Coordinates] = subg4[Coordinates]./nm

    # THX at JF!!!
    # https://github.com/j-fu/GridVisualize.jl/blob/1f2b299a436b7750702ccca282fa14152d80ebf9/src/pyplot.jl#L86
    function tridata(grid::ExtendableGrid)
        coord = grid[Coordinates]
        cellnodes = Matrix(grid[CellNodes])
        coord[1, :], coord[2, :], transpose(cellnodes .- 1)
    end

    if IVDirection == "forw"
        if V == "end"
            # uniform: min: 0.004904686640594403  max: 868.4271359428168
            # Maxwell: min: 0.0062429228146097865  max: 877.2227218127446
            vmin = 1.0e-2; vmax = 1.0e2
        elseif V == "inival"
            # uniform: min: 8.600091971551417e-6  max: 608986.5350550369
            # Maxwell: min: 1.9062557051494214e-5  max: 2.9789428140320955e6
            vmin = 1.0e-3; vmax = 1.0e3
        end
    else
        vmin = 0.005; vmax = 868
    end

    tripcolor(tridata(subg1)..., vcat(nn1./np1...), norm=matplotlib[:colors][:LogNorm](vmin=vmin, vmax=vmax), shading="gouraud", cmap = color_map, rasterized=true)
    xlabel(" x [nm]", fontsize=17)
    ylabel(" y [nm]", fontsize=17)
    axis([-20, 770, 20, 800])
    title("Ratio nn/np -- V = $V")
    cbar = colorbar(orientation = "vertical", label = "\$ n_{\\mathrm{n}}/ n_{\\mathrm{p}} \$", extend="both", spacing = "proportional")
    cbar.ax.set_yscale("log")
    tight_layout()

    if saveFig
        savefig(datadir("2D-planar-ratio-nn-np-scanrate-$textSR-generation-$generation-$IVDirection-$V.pdf"))
    end

    #####################
    figure()
    tripcolor(tridata(subg2)..., vcat(nn2./np2...), norm=matplotlib[:colors][:LogNorm](vmin=vmin, vmax=vmax), shading="gouraud", cmap = color_map, rasterized=true)
    xlabel(" x [nm]", fontsize=17)
    ylabel(" y [nm]", fontsize=17)
    axis([-20, 770, 20, 800])
    title("Ratio nn/np -- V = $V")
    cbar = colorbar(orientation = "vertical", label = "\$ n_{\\mathrm{n}}/ n_{\\mathrm{p}} \$", extend="both", spacing = "proportional")
    cbar.ax.set_yscale("log")
    tight_layout()

    if saveFig
        savefig(datadir("2D-nanotextured-ampl-$textampl-ratio-nn-np-scanrate-$scanrate-generation-$generation-$IVDirection-$V.pdf"))
    end

    #####################
    figure()
    tripcolor(tridata(subg3)..., vcat(nn3./np3...), norm=matplotlib[:colors][:LogNorm](vmin=vmin, vmax=vmax), shading="gouraud", cmap = color_map, rasterized=true)
    xlabel(" x [nm]", fontsize=17)
    ylabel(" y [nm]", fontsize=17)
    axis([-20, 770, 20, 800])
    title("Ratio nn/np -- V = $V")
    cbar = colorbar(orientation = "vertical", label = "\$ n_{\\mathrm{n}}/ n_{\\mathrm{p}} \$", extend="both", spacing = "proportional")
    cbar.ax.set_yscale("log")
    tight_layout()

    if saveFig
        savefig(datadir("2D-nanotextured-ampl-$textampl2-ratio-nn-np-scanrate-$textSR-generation-$generation-$IVDirection-$V.pdf"))
    end

    #####################
    figure()
    tripcolor(tridata(subg4)..., vcat(nn4./np4...), norm=matplotlib[:colors][:LogNorm](vmin=vmin, vmax=vmax), shading="gouraud", cmap = color_map, rasterized=true)
    xlabel(" x [nm]", fontsize=17)
    ylabel(" y [nm]", fontsize=17)
    axis([-20, 770, 20, 800])
    title("Ratio nn/np -- V = $V")
    cbar = colorbar(orientation = "vertical", label = "\$ n_{\\mathrm{n}}/ n_{\\mathrm{p}} \$", extend="both", spacing = "proportional")
    cbar.ax.set_yscale("log")
    tight_layout()

    if saveFig
        savefig(datadir("2D-nanotextured-ampl-$textampl3-ratio-nn-np-scanrate-$textSR-generation-$generation-$IVDirection-$V.pdf"))
    end

    return nothing
end # main

end # module