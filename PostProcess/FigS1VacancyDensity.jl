#=

Code for visualizing spatial vacancy density

=#

module FigS1VacancyDensity

using PyPlot
using DelimitedFiles
using PyCall
using ChargeTransport
using ExtendableGrids
using TexturedPerovskiteSolarCells

# for convenience
datadir = TexturedPerovskiteSolarCells.datadir
scriptsdir = TexturedPerovskiteSolarCells.scriptsdir

# https://stackoverflow.com/questions/29443369/how-to-make-a-custom-colormap-using-pyplot-not-matplotlib-proper
@pyimport matplotlib.colors as matcolors

# https://github.com/BIDS/colormap/blob/master/parula.py
cm_data = [[0.2081, 0.1663, 0.5292], [0.2116238095, 0.1897809524, 0.5776761905],
 [0.212252381, 0.2137714286, 0.6269714286], [0.2081, 0.2386, 0.6770857143],
 [0.1959047619, 0.2644571429, 0.7279], [0.1707285714, 0.2919380952,
  0.779247619], [0.1252714286, 0.3242428571, 0.8302714286],
 [0.0591333333, 0.3598333333, 0.8683333333], [0.0116952381, 0.3875095238,
  0.8819571429], [0.0059571429, 0.4086142857, 0.8828428571],
 [0.0165142857, 0.4266, 0.8786333333], [0.032852381, 0.4430428571,
  0.8719571429], [0.0498142857, 0.4585714286, 0.8640571429],
 [0.0629333333, 0.4736904762, 0.8554380952], [0.0722666667, 0.4886666667,
  0.8467], [0.0779428571, 0.5039857143, 0.8383714286],
 [0.079347619, 0.5200238095, 0.8311809524], [0.0749428571, 0.5375428571,
  0.8262714286], [0.0640571429, 0.5569857143, 0.8239571429],
 [0.0487714286, 0.5772238095, 0.8228285714], [0.0343428571, 0.5965809524,
  0.819852381], [0.0265, 0.6137, 0.8135], [0.0238904762, 0.6286619048,
  0.8037619048], [0.0230904762, 0.6417857143, 0.7912666667],
 [0.0227714286, 0.6534857143, 0.7767571429], [0.0266619048, 0.6641952381,
  0.7607190476], [0.0383714286, 0.6742714286, 0.743552381],
 [0.0589714286, 0.6837571429, 0.7253857143],
 [0.0843, 0.6928333333, 0.7061666667], [0.1132952381, 0.7015, 0.6858571429],
 [0.1452714286, 0.7097571429, 0.6646285714], [0.1801333333, 0.7176571429,
  0.6424333333], [0.2178285714, 0.7250428571, 0.6192619048],
 [0.2586428571, 0.7317142857, 0.5954285714], [0.3021714286, 0.7376047619,
  0.5711857143], [0.3481666667, 0.7424333333, 0.5472666667],
 [0.3952571429, 0.7459, 0.5244428571], [0.4420095238, 0.7480809524,
  0.5033142857], [0.4871238095, 0.7490619048, 0.4839761905],
 [0.5300285714, 0.7491142857, 0.4661142857], [0.5708571429, 0.7485190476,
  0.4493904762], [0.609852381, 0.7473142857, 0.4336857143],
 [0.6473, 0.7456, 0.4188], [0.6834190476, 0.7434761905, 0.4044333333],
 [0.7184095238, 0.7411333333, 0.3904761905],
 [0.7524857143, 0.7384, 0.3768142857], [0.7858428571, 0.7355666667,
  0.3632714286], [0.8185047619, 0.7327333333, 0.3497904762],
 [0.8506571429, 0.7299, 0.3360285714], [0.8824333333, 0.7274333333, 0.3217],
 [0.9139333333, 0.7257857143, 0.3062761905], [0.9449571429, 0.7261142857,
  0.2886428571], [0.9738952381, 0.7313952381, 0.266647619],
 [0.9937714286, 0.7454571429, 0.240347619], [0.9990428571, 0.7653142857,
  0.2164142857], [0.9955333333, 0.7860571429, 0.196652381],
 [0.988, 0.8066, 0.1793666667], [0.9788571429, 0.8271428571, 0.1633142857],
 [0.9697, 0.8481380952, 0.147452381], [0.9625857143, 0.8705142857, 0.1309],
 [0.9588714286, 0.8949, 0.1132428571], [0.9598238095, 0.9218333333,
  0.0948380952], [0.9661, 0.9514428571, 0.0755333333],
 [0.9763, 0.9831, 0.0538]]

parula_map = matcolors.LinearSegmentedColormap.from_list("parula", cm_data)

include(scriptsdir("SingleJunction.jl"))

function main(;scanrate  = 1000.0,   # "10p0" # "0p001"
              typeReco  = "all",    # "radiative"
              generationUniform = false,
              IVDirection = "forw", # "rev" #
              V = "inival", # "end"
              printText = true,
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

    Ea1           = data1.params.bandEdgeEnergy[iphia, :]; Ea2           = data2.params.bandEdgeEnergy[iphia, :]
    Ea3           = data3.params.bandEdgeEnergy[iphia, :]; Ea4           = data4.params.bandEdgeEnergy[iphia, :]
    ################################################################################
    ################################################################################

    coord1  = subg1[Coordinates]; coord2  = subg2[Coordinates]
    coord3  = subg3[Coordinates]; coord4  = subg4[Coordinates]

    coord1  = coord1./nm; coord2  = coord2./nm
    coord3  = coord3./nm; coord4  = coord4./nm

    #################################################################
    ## Plotting
    #################################################################

    # THX at JF!!!
    # https://github.com/j-fu/GridVisualize.jl/blob/1f2b299a436b7750702ccca282fa14152d80ebf9/src/pyplot.jl#L86
    function tridata(grid::ExtendableGrid)
        coord = grid[Coordinates]
        cellnodes = Matrix(grid[CellNodes])
        coord[1, :], coord[2, :], transpose(cellnodes .- 1)
    end

    #########################################################################################################
    #########################################################################################################
    #########################################################################################################

    sol1 = readdlm(datadir("sol", "$pathSol/Sol-2D-$IVDirection-planar-generation-$generation-reco-$typeReco-$V.dat"))'
    sol2 = readdlm(datadir("sol", "$pathSol/Sol-2D-$IVDirection-nanotextured-ampl-$textampl-generation-$generation-reco-$typeReco-$V.dat"))'
    sol3 = readdlm(datadir("sol", "$pathSol/Sol-2D-$IVDirection-nanotextured-ampl-$textampl2-generation-$generation-reco-$typeReco-$V.dat"))'
    sol4 = readdlm(datadir("sol", "$pathSol/Sol-2D-$IVDirection-nanotextured-ampl-$textampl3-generation-$generation-reco-$typeReco-$V.dat"))'

    na1  = Na[regionPero] .* Fcc[iphia].(za*( q*(view(sol1[iphia, :], subg1).-view(sol1[ipsi, :], subg1)) .+ Ea1[regionPero])./(kB*T))
    na2  = Na[regionPero] .* Fcc[iphia].(za*( q*(view(sol2[iphia, :], subg2).-view(sol2[ipsi, :], subg2)) .+ Ea2[regionPero])./(kB*T))
    na3  = Na[regionPero] .* Fcc[iphia].(za*( q*(view(sol3[iphia, :], subg3).-view(sol3[ipsi, :], subg3)) .+ Ea3[regionPero])./(kB*T))
    na4  = Na[regionPero] .* Fcc[iphia].(za*( q*(view(sol4[iphia, :], subg4).-view(sol4[ipsi, :], subg4)) .+ Ea4[regionPero])./(kB*T))

    # println(" ")
    # @show minimum(na1), maximum(na1)
    # @show minimum(na2), maximum(na2)
    # @show minimum(na3), maximum(na3)
    # @show minimum(na4), maximum(na4)

    #########################################################################
    #########################################################################

    if printText


        ## calculate the measure of the perovskite region
        mOmega = 0.0
        for icellVol in subg1[CellVolumes]
            mOmega = mOmega + icellVol
        end

        Int = ChargeTransport.integrate(ctsys1, storage!, sol1)./q
        println("Avg na for planar is:                     $(Int[iphia, regionPero]/mOmega)")

        ####################################################

        ## calculate the measure of the perovskite region
        mOmega = 0.0
        for icellVol in subg2[CellVolumes]
            mOmega = mOmega + icellVol
        end

        Int = ChargeTransport.integrate(ctsys2, storage!, sol2)./q
        println("Avg na for textured ($textampl nm) is:       $(Int[iphia, regionPero]/mOmega)")

        #####################################################

        ## calculate the measure of the perovskite region
        mOmega = 0.0
        for icellVol in subg3[CellVolumes]
            mOmega = mOmega + icellVol
        end

        Int = ChargeTransport.integrate(ctsys3, storage!, sol3)./q
        println("Avg na for textured ($textampl2 nm) is:       $(Int[iphia, regionPero]/mOmega)")

        #####################################################

        ## calculate the measure of the perovskite region
        mOmega = 0.0
        for icellVol in subg4[CellVolumes]
            mOmega = mOmega + icellVol
        end

        Int = ChargeTransport.integrate(ctsys4, storage!, sol4)./q
        println("Avg na for textured ($textampl3 nm) is:       $(Int[iphia, regionPero]/mOmega)")

    end

    #########################################################################
    #########################################################################

    subg1[Coordinates] = subg1[Coordinates]./nm; subg2[Coordinates] = subg2[Coordinates]./nm
    subg3[Coordinates] = subg3[Coordinates]./nm; subg4[Coordinates] = subg4[Coordinates]./nm

    #################################################################
    ## Plotting
    #################################################################

    if V == "inival"
        if IVDirection == "forw"
            # min: 1.6197422779145101e19   max: 7.185201519359629e23
            vmin = 3.0e20; vmax = 5.5e23
        elseif IVDirection == "rev"
            # min: 5.0316508005117356e20   max: 3.368722968477364e23
            vmin = 3.0e20; vmax = 5.5e23
        end
    end

    if V == "end" && IVDirection == "forw"
         # min: 7.73268508735298e18   max: 4.6453505984481594e23
    end

    vmin = 7.5e18; vmax = 7.0e23

    tripcolor(tridata(subg1)..., vcat(na1...), norm=matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax), shading="gouraud", cmap = parula_map, rasterized=true)
    xlabel("\$ x\$ [nm]", fontsize=17)
    ylabel("\$ y\$ [nm]", fontsize=17)
    axis([-20, 770, 20, 800])
    title("na -- V = $V")
    cbar = colorbar(orientation = "vertical", label = " \$ n_\\mathrm{a} \$  [\$\\mathrm{m}^{-3}\$]", extend="both")
    tight_layout()

    if saveFig
        savefig(datadir("2D-planar-na-scanrate-$textSR-generation-$generation-$IVDirection-$V.pdf"))
    end

    #####################
    figure()
    tripcolor(tridata(subg2)..., vcat(na2...), norm=matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax), shading="gouraud",  cmap = parula_map, rasterized=true)
    xlabel(" \$x\$ [nm]", fontsize=17)
    ylabel(" \$y \$ [nm]", fontsize=17)
    axis([-20, 770, 20, 800])
    title("na -- V = $V")
    cbar = colorbar(orientation = "vertical", label = " \$ n_\\mathrm{a} \$  [\$\\mathrm{m}^{-3}\$]", extend="both")
    tight_layout()

    if saveFig
        savefig(datadir("2D-nanotextured-ampl-$textampl-na-scanrate-$textSR-generation-$generation-$IVDirection-$V.pdf"))
    end

    #####################
    figure()
    tripcolor(tridata(subg3)..., vcat(na3...), norm=matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax), shading="gouraud", cmap = parula_map, rasterized=true)
    xlabel("\$ x\$ [nm]", fontsize=17)
    ylabel("\$ y\$ [nm]", fontsize=17)
    axis([-20, 770, 20, 800])
    title("na -- V = V")
    cbar = colorbar(orientation = "vertical", label = " \$ n_\\mathrm{a} \$  [\$\\mathrm{m}^{-3}\$]", extend="both")
    tight_layout()

    if saveFig
        savefig(datadir("2D-nanotextured-ampl-$textampl2-na-scanrate-$textSR-generation-$generation-$IVDirection-$V.pdf"))
    end

    #####################
    figure()
    tripcolor(tridata(subg4)..., vcat(na4...), norm=matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax), shading="gouraud", cmap = parula_map, rasterized=true)
    xlabel(" \$x\$ [nm]", fontsize=17)
    ylabel("\$ y\$ [nm]", fontsize=17)
    axis([-20, 770, 20, 800])
    title("na -- V = $V")
    cbar = colorbar(orientation = "vertical", label = " \$ n_\\mathrm{a} \$  [\$\\mathrm{m}^{-3}\$]", extend="both")
    tight_layout()

    if saveFig
        savefig(datadir("2D-nanotextured-ampl-$textampl3-na-scanrate-$textSR-generation-$generation-$IVDirection-$V.pdf"))
    end

    return nothing
end

end