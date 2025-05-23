
#=

Code for electron and hole density

=#

module FigS3ElectronAndHoleDensity

using PyPlot
using DelimitedFiles
using ChargeTransport
using PyCall
using ExtendableGrids
using Roots
using GridVisualize
using VoronoiFVM
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
              generationOn      = true,
              IVDirection       = "forw",   # "rev" #
              V                 = "inival", # "end"
              plotElectrons     = true,
              plotHoles         = true,
              printText         = true,
              saveFig           = false,
              parameter_file = scriptsdir("params_single_junction.jl"))

    include(parameter_file)

    if generationUniform
        generation = "uniform"
    else
        generation = "Maxwell"
    end

    if generationOn == false
        generation = "none"
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
    ########
    grid4, ctsys4 = SingleJunction.main(gridDim = 2, typeGrid = "nanotextured", amplitude = amplitude3, generation = true, generationUniform = generationUniform, plotPostProcess = true, printText = false, plotting = false)
    subg4         = subgrid(grid4, [regionPero])
    data4         = ctsys4.fvmsys.physics.data

    ################################################################################
    ################################################################################

    coord1 = subg1[Coordinates]; coord1  = coord1./nm; coord2 = subg2[Coordinates]; coord2  = coord2./nm
    coord3 = subg3[Coordinates]; coord3  = coord3./nm; coord4 = subg4[Coordinates]; coord4  = coord4./nm

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

    #########################################################################################################
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
    ######################
    np1  = Np[regionPero] .* Fcc[iphip].(zp*( q*(view(sol1[iphip, :], subg1).-view(sol1[ipsi, :], subg1)) .+ Ep[regionPero])./(kB*T))
    np2  = Np[regionPero] .* Fcc[iphip].(zp*( q*(view(sol2[iphip, :], subg2).-view(sol2[ipsi, :], subg2)) .+ Ep[regionPero])./(kB*T))
    np3  = Np[regionPero] .* Fcc[iphip].(zp*( q*(view(sol3[iphip, :], subg3).-view(sol3[ipsi, :], subg3)) .+ Ep[regionPero])./(kB*T))
    np4  = Np[regionPero] .* Fcc[iphip].(zp*( q*(view(sol4[iphip, :], subg4).-view(sol4[ipsi, :], subg4)) .+ Ep[regionPero])./(kB*T))

    # println(" ")
    # @show minimum(nn1), maximum(nn1)
    # @show minimum(nn2), maximum(nn2)
    # @show minimum(nn3), maximum(nn3)
    # @show minimum(nn4), maximum(nn4)

    # println(" ")
    # @show minimum(np1), maximum(np1)
    # @show minimum(np2), maximum(np2)
    # @show minimum(np3), maximum(np3)
    # @show minimum(np4), maximum(np4)

    #################################################################
    ## Plotting
    #################################################################

    if generationOn
        if IVDirection == "forw" # V = 0.0
            if V == "inival"
                if generationUniform
                    # Electrons: min: 9.383691713026936e14   max: 4.074824716536279e21
                    # Holes:     min: 1.0689897680120492e15  max: 7.751197277681249e21
                    vmin = 1.0e16; vmax = 6.5e21
                else
                    # Electrons: min: 4.5737717076138444e14   max: 4.0761049484493864e21
                    # Holes:     min: 5.021042521379344e14    max: 7.769599864126911e21
                    vmin = 1.0e16; vmax = 6.5e21
                end
            elseif V == "end"
                if generationUniform
                    # Electrons: min: 3.819303789636114e19 max: 1.5316010835957468e22
                    # Holes:     min: 1.763649499428583e19  max: 7.787049549761347e21
                    vmin = 5.0e19; vmax = 7.5e21
                else
                    # Electrons: min: 4.6326199582600585e19 max: 1.5313525444898185e22
                    # Holes:     min: 1.7456827170702346e19  max: 7.804517294253407e21
                    vmin = 5.0e19; vmax = 7.5e21
                end
            end
        elseif IVDirection == "rev" # V = 1.2
            # Electrons: min: 3.805848508591096e19   max: 1.6637739539610584e22
            # Holes:     min: 1.6117683513696866e19  max: 7.796118401656785e21
            vmin = 4.0e19; vmax = 1.0e22
        end
    else
        if IVDirection == "forw" # V = 0.0
            # Electrons: min: 0.2590123184772801   max: 4.064443583988271e21
            # Holes:     min: 0.49325902188972215  max: 7.740424601084793e21
            vmin = 1.0e-1; vmax = 7.5e21
        elseif IVDirection == "rev" # V = 1.2
            # Electrons: min: 2.7792349630924354e19  max: 1.6602032666764819e22
            # Holes:     min: 1.5247139300870613e19  max: 7.785628282134006e21
            vmin = 1.6e19; vmax = 1.5e22
        end
    end

    if plotElectrons

        tripcolor(tridata(subg1)..., vcat(nn1...), norm=matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax), shading="gouraud", cmap = parula_map, rasterized=true)
        xlabel(" \$x\$ [nm]", fontsize=17)
        ylabel(" \$y\$ [nm]", fontsize=17)
        axis([-20, 770, 20, 800])
        title("nn -- V = $V ($IVDirection)")
        cbar = colorbar(orientation = "vertical", label = " Density [\$\\mathrm{m}^{-3}\$]", extend="both")
        tight_layout()

        if saveFig
            savefig(datadir("2D-planar-nn-scanrate-$textSR-generation-$generation-$IVDirection-$V.pdf"))
        end

        #####################
        figure()
        tripcolor(tridata(subg2)..., vcat(nn2...), norm=matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax), shading="gouraud",  cmap = parula_map, rasterized=true)
        xlabel(" \$x\$ [nm]", fontsize=17)
        ylabel("\$ y\$ [nm]", fontsize=17)
        axis([-20, 770, 20, 800])
        title("nn -- V = $V ($IVDirection)")
        cbar = colorbar(orientation = "vertical", label = " Density [\$\\mathrm{m}^{-3}\$]", extend="both")
        tight_layout()

        if saveFig
            savefig(datadir("2D-nanotextured-ampl-$textampl-nn-scanrate-$textSR-generation-$generation-$IVDirection-$V.pdf"))
        end

        #####################
        figure()
        tripcolor(tridata(subg3)..., vcat(nn3...), norm=matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax), shading="gouraud", cmap = parula_map, rasterized=true)
        xlabel("\$ x\$ [nm]", fontsize=17)
        ylabel("\$ y\$ [nm]", fontsize=17)
        axis([-20, 770, 20, 800])
        title("nn -- V = $V ($IVDirection)")
        cbar = colorbar(orientation = "vertical", label = " Density [\$\\mathrm{m}^{-3}\$]", extend="both")
        tight_layout()

        if saveFig
            savefig(datadir("2D-nanotextured-ampl-$textampl2-nn-scanrate-$textSR-generation-$generation-$IVDirection-$V.pdf"))
        end

        #####################
        figure()
        tripcolor(tridata(subg4)..., vcat(nn4...), norm=matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax), shading="gouraud", cmap = parula_map, rasterized=true)
        xlabel("\$ x\$ [nm]", fontsize=17)
        ylabel("\$ y\$ [nm]", fontsize=17)
        axis([-20, 770, 20, 800])
        title("nn -- V =$V ($IVDirection)")
        cbar = colorbar(orientation = "vertical", label = " Density [\$\\mathrm{m}^{-3}\$]", extend="both")
        tight_layout()

        if saveFig
            savefig(datadir("2D-nanotextured-ampl-$textampl3-nn-scanrate-$textSR-generation-$generation-$IVDirection-$V.pdf"))
        end

    end # plotElectrons

    ####################################################################################################
    ####################################################################################################

    if plotHoles

        figure()
        tripcolor(tridata(subg1)..., vcat(np1...), norm=matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax), shading="gouraud", cmap = parula_map, rasterized=true)
        xlabel(" \$x\$ [nm]", fontsize=17)
        ylabel("\$ y\$ [nm]", fontsize=17)
        axis([-20, 770, 20, 800])
        title("np -- V = $V ($IVDirection)")
        cbar = colorbar(orientation = "vertical", label = " Density [\$\\mathrm{m}^{-3}\$]", extend="both")
        tight_layout()

        if saveFig
            savefig(datadir("2D-planar-np-scanrate-$textSR-generation-$generation-$IVDirection-$V.pdf"))
        end

        #####################
        figure()
        tripcolor(tridata(subg2)..., vcat(np2...), norm=matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax), shading="gouraud",  cmap = parula_map, rasterized=true)
        xlabel("\$ x\$ [nm]", fontsize=17)
        ylabel(" \$y\$ [nm]", fontsize=17)
        axis([-20, 770, 20, 800])
        title("np -- V = $V ($IVDirection)")
        cbar = colorbar(orientation = "vertical", label = " Density [\$\\mathrm{m}^{-3}\$]", extend="both")
        tight_layout()

        if saveFig
            savefig(datadir("2D-nanotextured-ampl-$textampl-np-scanrate-$textSR-generation-$generation-$IVDirection-$V.pdf"))
        end

        #####################
        figure()
        tripcolor(tridata(subg3)..., vcat(np3...), norm=matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax), shading="gouraud", cmap = parula_map, rasterized=true)
        xlabel(" \$x\$ [nm]", fontsize=17)
        ylabel(" \$y\$ [nm]", fontsize=17)
        axis([-20, 770, 20, 800])
        title("np -- V = $V ($IVDirection)")
        cbar = colorbar(orientation = "vertical", label = " Density [\$\\mathrm{m}^{-3}\$]", extend="both")
        tight_layout()

        if saveFig
            savefig(datadir("2D-nanotextured-ampl-$textampl2-np-scanrate-$textSR-generation-$generation-$IVDirection-$V.pdf"))
        end

        #####################
        figure()
        tripcolor(tridata(subg4)..., vcat(np4...), norm=matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax), shading="gouraud", cmap = parula_map, rasterized=true)
        xlabel("\$ x\$ [nm]", fontsize=17)
        ylabel("\$ y\$ [nm]", fontsize=17)
        axis([-20, 770, 20, 800])
        title("np -- V = $V ($IVDirection)")
        cbar = colorbar(orientation = "vertical", label = " Density [\$\\mathrm{m}^{-3}\$]", extend="both")
        tight_layout()

        if saveFig
            savefig(datadir("2D-nanotextured-ampl-$textampl3-np-scanrate-$textSR-generation-$generation-$IVDirection-$V.pdf"))
        end


    end # plotHoles

    ###################################################################
    ###################################################################

    ## calculate the measure of the perovskite region
    mOmega = 0.0
    for icellVol in subg1[CellVolumes]
        mOmega = mOmega + icellVol*nm*nm
    end

    Int1     =  ChargeTransport.integrate(ctsys1, storage!, sol1)./q
    nnAvg1   = -Int1[iphin, regionPero]/mOmega
    npAvg1   =  Int1[iphip, regionPero]/mOmega

    IntProd1 = ChargeTransport.integrate(ctsys1, DensityProduct, sol1)[1, regionPero]/mOmega

    #####################################################

    ## calculate the measure of the perovskite region
    mOmega = 0.0
    for icellVol in subg2[CellVolumes]
        mOmega = mOmega + icellVol*nm*nm
    end

    Int2     = ChargeTransport.integrate(ctsys2, storage!, sol2)./q
    nnAvg2   = -Int2[iphin, regionPero]/mOmega
    npAvg2   =  Int2[iphip, regionPero]/mOmega

    IntProd2 = ChargeTransport.integrate(ctsys2, DensityProduct, sol2)[1, regionPero]/mOmega

    #####################################################

    ## calculate the measure of the perovskite region
    mOmega = 0.0
    for icellVol in subg3[CellVolumes]
        mOmega = mOmega + icellVol*nm*nm
    end

    Int3     = ChargeTransport.integrate(ctsys3, storage!, sol3)./q
    nnAvg3   = -Int3[iphin, regionPero]/mOmega
    npAvg3   =  Int3[iphip, regionPero]/mOmega

    IntProd3 = ChargeTransport.integrate(ctsys3, DensityProduct, sol3)[1, regionPero]/mOmega

    #####################################################

    ## calculate the measure of the perovskite region
    mOmega = 0.0
    for icellVol in subg4[CellVolumes]
        mOmega = mOmega + icellVol*nm*nm
    end

    Int4     = ChargeTransport.integrate(ctsys4, storage!, sol4)./q
    nnAvg4   = -Int4[iphin, regionPero]/mOmega
    npAvg4   =  Int4[iphip, regionPero]/mOmega

    IntProd4 = ChargeTransport.integrate(ctsys4, DensityProduct, sol4)[1, regionPero]/mOmega

    #####################################################
    #####################################################
    #####################################################

    if printText
        println("  ")
        println("Avg nn for planar is:                $(nnAvg1)")
        println("Avg nn for textured ($textampl m) is:   $(nnAvg2)")
        println("Avg nn for textured ($textampl2 m) is:   $(nnAvg3)")
        println("Avg nn for textured ($textampl3 m) is:   $(nnAvg4)")

        println(" ")
        println("Avg np for planar is:                 $(npAvg1)")
        println("Avg np for textured ($textampl m) is:   $(npAvg2)")
        println("Avg np for textured ($textampl2 m) is:   $(npAvg3)")
        println("Avg np for textured ($textampl3 m) is:   $(npAvg4)")

        println(" ")

        println("Sqrt(nn * np) for planar is:               $(sqrt(IntProd1))")
        println("Sqrt(nn * np) for textured ($textampl m) is:  $(sqrt(IntProd2))")
        println("Sqrt(nn * np) for textured ($textampl2 m) is:  $(sqrt(IntProd3))")
        println("Sqrt(nn * np) for textured ($textampl3 m) is:  $(sqrt(IntProd4))")
    end

    return nothing
end

end