#=

Code for visualizing spatial band-edges

=#
module FigS5BandEdges

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
matcolors = pyimport("matplotlib.colors")

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

function main(;scanrate         = 1000.0,   # "10p0" # "0p001"
              typeReco          = "all",    # "radiative"
              generationUniform = false,
              IVDirection       = "forw", # "rev" #
              V                 = "end",  # "inival",
              printText         = true,
              saveFig           = false,
              parameter_file = scriptsdir("params_single_junction.jl"))

    include(parameter_file)

    PyPlot.rc("font", family="sans-serif", size=14)
    PyPlot.rc("mathtext", fontset="dejavusans")
    PyPlot.close("all")

    if generationUniform
        generation = "uniform"
    else
        generation = "Maxwell"
    end

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
    subg1        = subgrid(grid1, [regionPero]); subgn1 = subgrid(grid1, [regionETL1]); subgp1 = subgrid(grid1, [regionHTL])
    data1        = ctsys1.fvmsys.physics.data
    ########
    grid2, ctsys2 = SingleJunction.main(gridDim = 2, typeGrid = "nanotextured", amplitude = amplitude, generation = true, generationUniform = generationUniform, plotPostProcess = true, printText = false, plotting = false)
    subg2        = subgrid(grid2, [regionPero]); subgn2 = subgrid(grid2, [regionETL1]); subgp2 = subgrid(grid2, [regionHTL])
    data2        = ctsys2.fvmsys.physics.data
    ########
    grid3, ctsys3 = SingleJunction.main(gridDim = 2, typeGrid = "nanotextured", amplitude = amplitude2, generation = true, generationUniform = generationUniform, plotPostProcess = true, printText = false, plotting = false)
    subg3        = subgrid(grid3, [regionPero]); subgn3 = subgrid(grid3, [regionETL1]); subgp3 = subgrid(grid3, [regionHTL])
    data3        = ctsys3.fvmsys.physics.data
    ########
    grid4, ctsys4 = SingleJunction.main(gridDim = 2, typeGrid = "nanotextured", amplitude = amplitude3, generation = true, generationUniform = generationUniform, plotPostProcess = true, printText = false, plotting = false)
    subg4         = subgrid(grid4, [regionPero]); subgn4 = subgrid(grid4, [regionETL1]); subgp4 = subgrid(grid4, [regionHTL])
    data4         = ctsys4.fvmsys.physics.data
    #################################################################
    #################################################################

    coord1  = subg1[Coordinates]; coord2 = subg2[Coordinates]
    coord3  = subg3[Coordinates]; coord4 = subg4[Coordinates]

    coord1  = coord1./nm; coord2 = coord2./nm
    coord3  = coord3./nm; coord4 = coord4./nm

    #################################################################
    ## Plotting
    #################################################################

    subg1[Coordinates] = subg1[Coordinates]./nm
    subg2[Coordinates] = subg2[Coordinates]./nm
    subg3[Coordinates] = subg3[Coordinates]./nm
    subg4[Coordinates] = subg4[Coordinates]./nm

    # THX at JF!!!
    # https://github.com/j-fu/GridVisualize.jl/blob/1f2b299a436b7750702ccca282fa14152d80ebf9/src/pyplot.jl#L86
    function tridata(grid::ExtendableGrid)
        coordx = grid[Coordinates]
        cellnodes = Matrix(grid[CellNodes])
        coordx[1, :], coordx[2, :], transpose(cellnodes .- 1)
    end

    #########################################################################################################
    #########################################################################################################
    #########################################################################################################

    #################################################################
    ## read in solution
    #################################################################

    sol1    = readdlm(datadir("sol", "$pathSol/Sol-2D-$IVDirection-planar-generation-$generation-reco-$typeReco-$V.dat"))'
    sol2    = readdlm(datadir("sol", "$pathSol/Sol-2D-$IVDirection-nanotextured-ampl-$textampl-generation-$generation-reco-$typeReco-$V.dat"))'
    sol3    = readdlm(datadir("sol", "$pathSol/Sol-2D-$IVDirection-nanotextured-ampl-$textampl2-generation-$generation-reco-$typeReco-$V.dat"))'
    sol4    = readdlm(datadir("sol", "$pathSol/Sol-2D-$IVDirection-nanotextured-ampl-$textampl3-generation-$generation-reco-$typeReco-$V.dat"))'

    psi1    = view(sol1[ipsi, :], subg1); psi2 = view(sol2[ipsi, :], subg2)
    psi3    = view(sol3[ipsi, :], subg3); psi4 = view(sol4[ipsi, :], subg4)

    ###########################
    Ec1       = En[regionPero]./q .- psi1; Ec2 = En[regionPero]./q .- psi2
    Ec3       = En[regionPero]./q .- psi3; Ec4 = En[regionPero]./q .- psi4

    Ev1       = Ep[regionPero]./q .- psi1; Ev2 = Ep[regionPero]./q .- psi2
    Ev3       = Ep[regionPero]./q .- psi3; Ev4 = Ep[regionPero]./q .- psi4

    # println(" ")
    # @show minimum(Ec1), maximum(Ec1)
    # @show minimum(Ec2), maximum(Ec2)
    # @show minimum(Ec3), maximum(Ec3)
    # @show minimum(Ec4), maximum(Ec4)

    # println(" ")
    # @show minimum(Ev1), maximum(Ev1)
    # @show minimum(Ev2), maximum(Ev2)
    # @show minimum(Ev3), maximum(Ev3)
    # @show minimum(Ev4), maximum(Ev4)

    #################################################################
    ## Plotting
    #################################################################

    if IVDirection == "forw"  && V == "end"
        ######### uniform #########
        # Electrons: min: 0.12839615000670834  max: 0.2841084898473256
        # Holes:     min: -1.5016038499932929  max: -1.3458915101526756
        ######### Maxwell #########
        # Electrons: min: 0.12840019434854044  max: 0.28407871131014106
        # Holes:     min: -1.5015998056514608  max: -1.3459212886898602
        vminE = 0.13; vmaxE = 0.28
        vminH = -1.5; vmaxH = -1.35
    else
        vminE = 0.13; vmaxE = 0.28
        vminH = -1.5; vmaxH = -1.35
    end

    tripcolor(tridata(subg1)..., vcat(Ec1...), vmin=vminE, vmax=vmaxE, shading="gouraud", cmap = parula_map, rasterized=true) #
    xlabel("\$ x\$ [nm]", fontsize=17)
    ylabel("\$ y\$ [nm]", fontsize=17)
    axis([-20, 770, 20, 800])
    title("Ec -- V = $V")
    cbar = colorbar(orientation = "vertical", label = "\$ E_{\\mathrm{c}}\$ [eV]", extend="both")
    tight_layout()

    if saveFig
        savefig(datadir("2D-planar-Ec-scanrate-$textSR-generation-$generation-$IVDirection-$V.pdf"))
    end

    #####################
    figure()
    tripcolor(tridata(subg2)..., vcat(Ec2...), vmin=vminE, vmax=vmaxE, shading="gouraud", cmap = parula_map, rasterized=true)
    xlabel("\$ x\$ [nm]", fontsize=17)
    ylabel(" \$y\$ [nm]", fontsize=17)
    axis([-20, 770, 20, 800])
    title("Ec -- V = $V")
    cbar = colorbar(orientation = "vertical", label = "\$ E_{\\mathrm{c}}\$ [eV]", extend="both")
    tight_layout()

    if saveFig
        savefig(datadir("2D-nanotextured-ampl-$textampl-Ec-scanrate-$textSR-generation-$generation-$IVDirection-$V.pdf"))
    end

    #####################
    figure()
    tripcolor(tridata(subg3)..., vcat(Ec3...), vmin=vminE, vmax=vmaxE, shading="gouraud", cmap = parula_map, rasterized=true)
    xlabel(" \$x\$ [nm]", fontsize=17)
    ylabel("\$ y\$ [nm]", fontsize=17)
    axis([-20, 770, 20, 800])
    title("Ec -- V = $V")
    cbar = colorbar(orientation = "vertical", label = "\$ E_{\\mathrm{c}}\$ [eV]", extend="both")
    tight_layout()

    if saveFig
        savefig(datadir("2D-nanotextured-ampl-$textampl2-Ec-scanrate-$textSR-generation-$generation-$IVDirection-$V.pdf"))
    end
    #####################
    figure()
    tripcolor(tridata(subg4)..., vcat(Ec4...), vmin=vminE, vmax=vmaxE, shading="gouraud", cmap = parula_map, rasterized=true)
    xlabel("\$ x\$ [nm]", fontsize=17)
    ylabel("\$ y\$ [nm]", fontsize=17)
    axis([-20, 770, 20, 800])
    title("Ec -- V = $V")
    cbar = colorbar(orientation = "vertical", label = "\$ E_{\\mathrm{c}}\$ [eV]", extend="both")
    tight_layout()

    if saveFig
        savefig(datadir("2D-nanotextured-ampl-$textampl3-Ec-scanrate-$textSR-generation-$generation-$IVDirection-$V.pdf"))
    end

    ##############################################################################
    ##############################################################################

    figure()
    tripcolor(tridata(subg1)..., vcat(Ev1...), vmin=vminH, vmax=vmaxH, shading="gouraud", cmap = parula_map, rasterized=true) #
    xlabel(" \$x\$ [nm]", fontsize=17)
    ylabel(" \$y\$ [nm]", fontsize=17)
    axis([-20, 770, 20, 800])
    title("Ev -- V = $V")
    cbar = colorbar(orientation = "vertical", label = "\$ E_{\\mathrm{v}}\$ [eV]", extend="both")
    tight_layout()

    if saveFig
        savefig(datadir("2D-planar-Ev-scanrate-$textSR-generation-$generation-$IVDirection-$V.pdf"))
    end

    #####################
    figure()
    tripcolor(tridata(subg2)..., vcat(Ev2...), vmin=vminH, vmax=vmaxH, shading="gouraud", cmap = parula_map, rasterized=true)
    xlabel(" \$x\$ [nm]", fontsize=17)
    ylabel("\$ y\$ [nm]", fontsize=17)
    axis([-20, 770, 20, 800])
    title("Ev -- V = $V")
    cbar = colorbar(orientation = "vertical", label = "\$ E_{\\mathrm{v}}\$ [eV]", extend="both")
    tight_layout()

    if saveFig
        savefig(datadir("2D-nanotextured-ampl-$textampl-Ev-scanrate-$textSR-generation-$generation-$IVDirection-$V.pdf"))
    end

    #####################
    figure()
    tripcolor(tridata(subg3)..., vcat(Ev3...), vmin=vminH, vmax=vmaxH, shading="gouraud", cmap = parula_map, rasterized=true)
    xlabel(" \$x\$ [nm]", fontsize=17)
    ylabel(" \$y\$ [nm]", fontsize=17)
    axis([-20, 770, 20, 800])
    title("Ev -- V = $V")
    cbar = colorbar(orientation = "vertical", label = "\$ E_{\\mathrm{v}}\$ [eV]", extend="both")
    tight_layout()

    if saveFig
        savefig(datadir("2D-nanotextured-ampl-$textampl2-Ev-scanrate-$textSR-generation-$generation-$IVDirection-$V.pdf"))
    end
    #####################
    figure()
    tripcolor(tridata(subg4)..., vcat(Ev4...), vmin=vminH, vmax=vmaxH, shading="gouraud", cmap = parula_map, rasterized=true)
    xlabel(" \$x\$ [nm]", fontsize=17)
    ylabel("\$ y\$ [nm]", fontsize=17)
    axis([-20, 770, 20, 800])
    title("Ev -- V = $V")
    cbar = colorbar(orientation = "vertical", label = "\$ E_{\\mathrm{v}}\$ [eV]", extend="both")
    tight_layout()

    if saveFig
        savefig(datadir("2D-nanotextured-ampl-$textampl3-Ev-scanrate-$textSR-generation-$generation-$IVDirection-$V.pdf"))
    end

   return nothing
end

end