#=

Code for electric field strength

=#

module FigS2FieldStrength

using PyPlot
using DelimitedFiles
using ChargeTransport
using PyCall
using ExtendableGrids
using VoronoiFVM
using LinearAlgebra
using TexturedPerovskiteSolarCells

# for convenience
datadir = TexturedPerovskiteSolarCells.datadir
scriptsdir = TexturedPerovskiteSolarCells.scriptsdir

# thx https://discourse.julialang.org/t/meshgrid-function-in-julia/48679/4?u=j-fu
function meshgrid(rc)
    nx = length(rc[1])
    ny = length(rc[2])
    xout = zeros(ny, nx)
    yout = zeros(ny, nx)
    for jx in 1:nx
        for ix in 1:ny
            xout[ix, jx] = rc[1][jx]
            yout[ix, jx] = rc[2][ix]
        end
    end
    return xout, yout
end

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
              generationOn      = true,
              IVDirection = "forw", # "rev" #
              V = "inival", # "end"
              saveFig = false,
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

    if generationOn == false
        generation = "none"
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
    subg1        = subgrid(grid1, [regionPero])
    data1        = ctsys1.fvmsys.physics.data
    ########
    grid2, ctsys2 = SingleJunction.main(gridDim = 2, typeGrid = "nanotextured", amplitude = amplitude, generation = true, generationUniform = generationUniform, plotPostProcess = true, printText = false, plotting = false)
    subg2        = subgrid(grid2, [regionPero])
    data2        = ctsys2.fvmsys.physics.data
    ########
    grid3, ctsys3 = SingleJunction.main(gridDim = 2, typeGrid = "nanotextured", amplitude = amplitude2, generation = true, generationUniform = generationUniform, plotPostProcess = true, printText = false, plotting = false)
    subg3        = subgrid(grid3, [regionPero])
    data3        = ctsys3.fvmsys.physics.data
    ########
    grid4, ctsys4 = SingleJunction.main(gridDim = 2, typeGrid = "nanotextured", amplitude = amplitude3, generation = true, generationUniform = generationUniform, plotPostProcess = true, printText = false, plotting = false)
    subg4         = subgrid(grid4, [regionPero])
    data4         = ctsys4.fvmsys.physics.data

    #################################################################
    #################################################################

    coord1  = subg1[Coordinates]; coord2  = subg2[Coordinates]
    coord3  = subg3[Coordinates]; coord4  = subg4[Coordinates]

    coord1  = coord1./nm; coord2  = coord2./nm
    coord3  = coord3./nm; coord4  = coord4./nm

    #################################################################
    ## Plotting
    #################################################################

    subg1[Coordinates] = subg1[Coordinates]./nm; subg2[Coordinates] = subg2[Coordinates]./nm
    subg3[Coordinates] = subg3[Coordinates]./nm; subg4[Coordinates] = subg4[Coordinates]./nm

    # THX at JF!!!
    # https://github.com/j-fu/GridVisualize.jl/blob/1f2b299a436b7750702ccca282fa14152d80ebf9/src/pyplot.jl#L86
    function tridata(grid::ExtendableGrid)
        coordx = grid[Coordinates]
        cellnodes = Matrix(grid[CellNodes])
        coordx[1, :], coordx[2, :], transpose(cellnodes .- 1)
    end

    #################################################################
    ## read in solution
    #################################################################

    sol1   = readdlm(datadir("sol", "$pathSol/Sol-2D-$IVDirection-planar-generation-$generation-reco-$typeReco-$V.dat"))'
    sol2   = readdlm(datadir("sol", "$pathSol/Sol-2D-$IVDirection-nanotextured-ampl-$textampl-generation-$generation-reco-$typeReco-$V.dat"))'
    sol3   = readdlm(datadir("sol", "$pathSol/Sol-2D-$IVDirection-nanotextured-ampl-$textampl2-generation-$generation-reco-$typeReco-$V.dat"))'
    sol4   = readdlm(datadir("sol", "$pathSol/Sol-2D-$IVDirection-nanotextured-ampl-$textampl3-generation-$generation-reco-$typeReco-$V.dat"))'

    nodes1 = subg1[NodeParents]; nodes2 = subg2[NodeParents]
    nodes3 = subg3[NodeParents]; nodes4 = subg4[NodeParents]

    nft1   = VoronoiFVM.nodeflux(ctsys1.fvmsys, sol1)
    jPsi1  = nft1[:, ipsi, nodes1]./ (ε0*εr[regionPero]); jPsi1Abs = norm.(eachcol(jPsi1))

    nft2   = VoronoiFVM.nodeflux(ctsys2.fvmsys, sol2)
    jPsi2  = nft2[:, ipsi, nodes2]./ (ε0*εr[regionPero]); jPsi2Abs = norm.(eachcol(jPsi2))

    nft3   = VoronoiFVM.nodeflux(ctsys3.fvmsys, sol3)
    jPsi3  = nft3[:, ipsi, nodes3]./ (ε0*εr[regionPero]); jPsi3Abs = norm.(eachcol(jPsi3))

    nft4   = VoronoiFVM.nodeflux(ctsys4.fvmsys, sol4)
    jPsi4  = nft4[:, ipsi, nodes4]./ (ε0*εr[regionPero]); jPsi4Abs = norm.(eachcol(jPsi4))

    # println("  ")
    # @show minimum(jPsi1Abs), maximum(jPsi1Abs)
    # @show minimum(jPsi2Abs), maximum(jPsi2Abs)
    # @show minimum(jPsi3Abs), maximum(jPsi3Abs)
    # @show minimum(jPsi4Abs), maximum(jPsi4Abs)

    #################################################################
    ## Plotting
    #################################################################

    if generationOn
        if IVDirection == "forw"
            if V == "inival"
                # uniform
                # min: 102711.25386983373  max: 1.5184528454051398e7
                # Maxwell
                # min: 102991.51867742644  max: 1.5184637602326512e7
                vmin = 1.03e5; vmax = 1.1e7
            elseif V== "end"
                # uniform
                # min: 75.49049284822992  max: 3.834495502745478e6
                # Maxwell
                # min: 74.96420127270625  max: 3.8348852946209046e6
                vmin = 1.03e5; vmax = 1.1e7
            end
        elseif IVDirection == "rev"
            # min: 16.346514954714937; # max: 3.681224035003072e6
            vmin = 1.0e4; vmax = 3.68e6
        end
    else
        if IVDirection == "forw" && V == "inival"
            # min: 102162.94200289267  max: 1.5185368956049463e7
            vmin = 1.03e5; vmax = 1.1e7
        else
            vmin = 1.0e4; vmax = 3.68e6
        end
    end

    tripcolor(tridata(subg1)..., vcat(jPsi1Abs...), norm=matplotlib[:colors][:LogNorm](vmin=vmin, vmax=vmax), shading="gouraud", cmap = parula_map, rasterized=true) #
    xlabel(" x [nm]", fontsize=17)
    ylabel(" y [nm]", fontsize=17)
    axis([-20, 770, 20, 800])
    title("E-Field -- V = $V")
    cbar = colorbar(orientation = "vertical", label = "El. field strength [\$\\frac{V}{\\mathrm{m}}\$]", extend="both")
    tight_layout()

    if saveFig
        savefig(datadir("2D-planar-E-field-scanrate-$textSR-generation-$generation-$IVDirection-$V.pdf"))
    end

    #####################
    figure()
    tripcolor(tridata(subg2)..., vcat(jPsi2Abs...), norm=matplotlib[:colors][:LogNorm](vmin=vmin, vmax=vmax), shading="gouraud", cmap = parula_map, rasterized=true)
    xlabel(" x [nm]", fontsize=17)
    ylabel(" y [nm]", fontsize=17)
    axis([-20, 770, 20, 800])
    title("E-Field -- V = $V")
    cbar = colorbar(orientation = "vertical", label = "El. field strength [\$\\frac{V}{\\mathrm{m}}\$]", extend="both")
    tight_layout()

    if saveFig
        savefig(datadir("2D-nanotextured-ampl-$textampl-E-field-scanrate-$textSR-generation-$generation-$IVDirection-$V.pdf"))
    end

    #####################
    figure()
    tripcolor(tridata(subg3)..., vcat(jPsi3Abs...), norm=matplotlib[:colors][:LogNorm](vmin=vmin, vmax=vmax), shading="gouraud", cmap = parula_map, rasterized=true)
    xlabel(" x [nm]", fontsize=17)
    ylabel(" y [nm]", fontsize=17)
    axis([-20, 770, 20, 800])
    title("E-Field -- V = $V")
    cbar = colorbar(orientation = "vertical", label = "El. field strength [\$\\frac{V}{\\mathrm{m}}\$]", extend="both")
    tight_layout()

    if saveFig
        savefig(datadir("2D-nanotextured-ampl-$textampl2-E-field-scanrate-$textSR-generation-$generation-$IVDirection-$V.pdf"))
    end

    #####################
    figure()
    tripcolor(tridata(subg4)..., vcat(jPsi4Abs...), norm=matplotlib[:colors][:LogNorm](vmin=vmin, vmax=vmax), shading="gouraud", cmap = parula_map, rasterized=true)
    xlabel(" x [nm]", fontsize=17)
    ylabel(" y [nm]", fontsize=17)
    axis([-20, 770, 20, 800])
    title("E-Field -- V = $V")
    cbar = colorbar(orientation = "vertical", label = "El. field strength [\$\\frac{V}{\\mathrm{m}}\$]", extend="both")
    tight_layout()

    if saveFig
        savefig(datadir("2D-nanotextured-ampl-$textampl3-E-field-scanrate-$textSR-generation-$generation-$IVDirection-$V.pdf"))
    end

end

end