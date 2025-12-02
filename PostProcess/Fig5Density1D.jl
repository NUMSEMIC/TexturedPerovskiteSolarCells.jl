#=

Code for visualizing carrier density cross section

=#

module Fig5Density1D

using PyPlot
using DelimitedFiles
using ChargeTransport
using ExtendableGrids
using TexturedPerovskiteSolarCells

# for convenience
datadir = TexturedPerovskiteSolarCells.datadir
scriptsdir = TexturedPerovskiteSolarCells.scriptsdir

include(scriptsdir("SingleJunction.jl"))

## for generating the data in manuscript, the options are:
# vETL = 10,   vHTL = 10,  and then V = "inival" or V = "V-1p2"
# vETL = 2000, vHTL = 500, and then V = "inival" or V = "end" (V = 1.2 V)
function main(;
        vETL = 10,  # in cm/s
        vHTL = 10,  # in cm/s
        V = "V-1p2", # "inival"
        saveFig = false,
        parameter_set = ParamsSingleJunction,
    )

    # use the destructuring operator to extract all the necessary parameters
    (; regionPero, regionETL, regionHTL, Nn, Np, zn, zp, iphin, iphip, ipsi, Fcc, En, Ep, T) = parameter_set()

    pathSol = "vETL=$(vETL)vHTL=$(vHTL)/"

    (; q, k_B) = ChargeTransport.constants

    @local_unitfactors nm

    PyPlot.rc("font", family = "sans-serif", size = 14)
    PyPlot.rc("mathtext", fontset = "dejavusans")
    PyPlot.close("all")

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
    subg1 = subgrid(grid1, [regionPero]); subgn1 = subgrid(grid1, [regionETL]); subgp1 = subgrid(grid1, [regionHTL])
    ########
    grid2, ctsys2 = SingleJunction.main(gridDim = 2, typeGrid = "nanotextured", amplitude = amplitude, generation = true, plotPostProcess = true, printText = false, plotting = false)
    subg2 = subgrid(grid2, [regionPero]); subgn2 = subgrid(grid2, [regionETL]); subgp2 = subgrid(grid2, [regionHTL])
    ########
    grid3, ctsys3 = SingleJunction.main(gridDim = 2, typeGrid = "nanotextured", amplitude = amplitude2, generation = true, plotPostProcess = true, printText = false, plotting = false)
    subg3 = subgrid(grid3, [regionPero]); subgn3 = subgrid(grid3, [regionETL]); subgp3 = subgrid(grid3, [regionHTL])
    ########

    #################################################################
    #################################################################

    coord1 = subg1[Coordinates]
    coord2 = subg2[Coordinates]
    coord3 = subg3[Coordinates]

    coord1 = coord1 ./ nm
    coord2 = coord2 ./ nm
    coord3 = coord3 ./ nm

    triang1 = matplotlib.tri.Triangulation(coord1[1, :], coord1[2, :], transpose(Matrix(subg1[CellNodes]) .- 1))
    triang2 = matplotlib.tri.Triangulation(coord2[1, :], coord2[2, :], transpose(Matrix(subg2[CellNodes]) .- 1))
    triang3 = matplotlib.tri.Triangulation(coord3[1, :], coord3[2, :], transpose(Matrix(subg3[CellNodes]) .- 1))

    ###########################################################
    ###########################################################

    coordn1 = subgn1[Coordinates]
    coordn2 = subgn2[Coordinates]
    coordn3 = subgn3[Coordinates]

    coordn1 = coordn1 ./ nm
    coordn2 = coordn2 ./ nm
    coordn3 = coordn3 ./ nm

    triangn1 = matplotlib.tri.Triangulation(coordn1[1, :], coordn1[2, :], transpose(Matrix(subgn1[CellNodes]) .- 1))
    triangn2 = matplotlib.tri.Triangulation(coordn2[1, :], coordn2[2, :], transpose(Matrix(subgn2[CellNodes]) .- 1))
    triangn3 = matplotlib.tri.Triangulation(coordn3[1, :], coordn3[2, :], transpose(Matrix(subgn3[CellNodes]) .- 1))

    ###########################################################
    ###########################################################

    coordp1 = subgp1[Coordinates]
    coordp2 = subgp2[Coordinates]
    coordp3 = subgp3[Coordinates]

    coordp1 = coordp1 ./ nm
    coordp2 = coordp2 ./ nm
    coordp3 = coordp3 ./ nm

    triangp1 = matplotlib.tri.Triangulation(coordp1[1, :], coordp1[2, :], transpose(Matrix(subgp1[CellNodes]) .- 1))
    triangp2 = matplotlib.tri.Triangulation(coordp2[1, :], coordp2[2, :], transpose(Matrix(subgp2[CellNodes]) .- 1))
    triangp3 = matplotlib.tri.Triangulation(coordp3[1, :], coordp3[2, :], transpose(Matrix(subgp3[CellNodes]) .- 1))

    #################################################################
    ## Plotting
    #################################################################

    subg1[Coordinates] = subg1[Coordinates] ./ nm
    subg2[Coordinates] = subg2[Coordinates] ./ nm
    subg3[Coordinates] = subg3[Coordinates] ./ nm

    #################################################################
    ## read in solution
    #################################################################

    sol1 = readdlm(datadir("sol", "$pathSol/Sol-2D-forw-planar-generation-Maxwell-reco-all-$V.dat"))'
    sol2 = readdlm(datadir("sol", "$pathSol/Sol-2D-forw-nanotextured-ampl-$textampl-generation-Maxwell-reco-all-$V.dat"))'
    sol3 = readdlm(datadir("sol", "$pathSol/Sol-2D-forw-nanotextured-ampl-$textampl2-generation-Maxwell-reco-all-$V.dat"))'

    ###########################################

    nn1 = Nn[regionPero] .* Fcc[iphin].(zn * (q * (view(sol1[iphin, :], subg1) .- view(sol1[ipsi, :], subg1)) .+ En[regionPero]) ./ (k_B * T))
    nn2 = Nn[regionPero] .* Fcc[iphin].(zn * (q * (view(sol2[iphin, :], subg2) .- view(sol2[ipsi, :], subg2)) .+ En[regionPero]) ./ (k_B * T))
    nn3 = Nn[regionPero] .* Fcc[iphin].(zn * (q * (view(sol3[iphin, :], subg3) .- view(sol3[ipsi, :], subg3)) .+ En[regionPero]) ./ (k_B * T))

    np1 = Np[regionPero] .* Fcc[iphip].(zp * (q * (view(sol1[iphip, :], subg1) .- view(sol1[ipsi, :], subg1)) .+ Ep[regionPero]) ./ (k_B * T))
    np2 = Np[regionPero] .* Fcc[iphip].(zp * (q * (view(sol2[iphip, :], subg2) .- view(sol2[ipsi, :], subg2)) .+ Ep[regionPero]) ./ (k_B * T))
    np3 = Np[regionPero] .* Fcc[iphip].(zp * (q * (view(sol3[iphip, :], subg3) .- view(sol3[ipsi, :], subg3)) .+ Ep[regionPero]) ./ (k_B * T))

    ######
    nn1itp = matplotlib.tri.LinearTriInterpolator(triang1, nn1)
    nn2itp = matplotlib.tri.LinearTriInterpolator(triang2, nn2)
    nn3itp = matplotlib.tri.LinearTriInterpolator(triang3, nn3)

    np1itp = matplotlib.tri.LinearTriInterpolator(triang1, np1)
    np2itp = matplotlib.tri.LinearTriInterpolator(triang2, np2)
    np3itp = matplotlib.tri.LinearTriInterpolator(triang3, np3)

    ##########################################
    nnn1 = Nn[regionETL] .* Fcc[iphin].(zn * (q * (view(sol1[iphin, :], subgn1) .- view(sol1[ipsi, :], subgn1)) .+ En[regionETL]) ./ (k_B * T))
    nnn2 = Nn[regionETL] .* Fcc[iphin].(zn * (q * (view(sol2[iphin, :], subgn2) .- view(sol2[ipsi, :], subgn2)) .+ En[regionETL]) ./ (k_B * T))
    nnn3 = Nn[regionETL] .* Fcc[iphin].(zn * (q * (view(sol3[iphin, :], subgn3) .- view(sol3[ipsi, :], subgn3)) .+ En[regionETL]) ./ (k_B * T))

    npn1 = Np[regionETL] .* Fcc[iphip].(zp * (q * (view(sol1[iphip, :], subgn1) .- view(sol1[ipsi, :], subgn1)) .+ Ep[regionETL]) ./ (k_B * T))
    npn2 = Np[regionETL] .* Fcc[iphip].(zp * (q * (view(sol2[iphip, :], subgn2) .- view(sol2[ipsi, :], subgn2)) .+ Ep[regionETL]) ./ (k_B * T))
    npn3 = Np[regionETL] .* Fcc[iphip].(zp * (q * (view(sol3[iphip, :], subgn3) .- view(sol3[ipsi, :], subgn3)) .+ Ep[regionETL]) ./ (k_B * T))

    ######
    nnn1itp = matplotlib.tri.LinearTriInterpolator(triangn1, nnn1)
    nnn2itp = matplotlib.tri.LinearTriInterpolator(triangn2, nnn2)
    nnn3itp = matplotlib.tri.LinearTriInterpolator(triangn3, nnn3)

    npn1itp = matplotlib.tri.LinearTriInterpolator(triangn1, npn1)
    npn2itp = matplotlib.tri.LinearTriInterpolator(triangn2, npn2)
    npn3itp = matplotlib.tri.LinearTriInterpolator(triangn3, npn3)

    ##########################################

    nnp1 = Nn[regionHTL] .* Fcc[iphin].(zn * (q * (view(sol1[iphin, :], subgp1) .- view(sol1[ipsi, :], subgp1)) .+ En[regionHTL]) ./ (k_B * T))
    nnp2 = Nn[regionHTL] .* Fcc[iphin].(zn * (q * (view(sol2[iphin, :], subgp2) .- view(sol2[ipsi, :], subgp2)) .+ En[regionHTL]) ./ (k_B * T))
    nnp3 = Nn[regionHTL] .* Fcc[iphin].(zn * (q * (view(sol3[iphin, :], subgp3) .- view(sol3[ipsi, :], subgp3)) .+ En[regionHTL]) ./ (k_B * T))

    npp1 = Np[regionHTL] .* Fcc[iphip].(zp * (q * (view(sol1[iphip, :], subgp1) .- view(sol1[ipsi, :], subgp1)) .+ Ep[regionHTL]) ./ (k_B * T))
    npp2 = Np[regionHTL] .* Fcc[iphip].(zp * (q * (view(sol2[iphip, :], subgp2) .- view(sol2[ipsi, :], subgp2)) .+ Ep[regionHTL]) ./ (k_B * T))
    npp3 = Np[regionHTL] .* Fcc[iphip].(zp * (q * (view(sol3[iphip, :], subgp3) .- view(sol3[ipsi, :], subgp3)) .+ Ep[regionHTL]) ./ (k_B * T))

    ######
    nnp1itp = matplotlib.tri.LinearTriInterpolator(triangp1, nnp1)
    nnp2itp = matplotlib.tri.LinearTriInterpolator(triangp2, nnp2)
    nnp3itp = matplotlib.tri.LinearTriInterpolator(triangp3, nnp3)

    npp1itp = matplotlib.tri.LinearTriInterpolator(triangp1, npp1)
    npp2itp = matplotlib.tri.LinearTriInterpolator(triangp2, npp2)
    npp3itp = matplotlib.tri.LinearTriInterpolator(triangp3, npp3)

    #################################################################
    ## Plotting
    #################################################################

    fCos(x, ampl) = 0.5 * ampl .* cos.(pi .+ 2 .* pi .* x ./ 750.0) .+ 0.5 * ampl .+ 30.0 + 400 - ampl / 2

    Blues = get_cmap(:Blues); Oranges = get_cmap(:Oranges)

    XVal = 187.0; nodes = 600
    XX = XVal .* ones(nodes)
    YY = collect(range(0.0, 450.0, length = nodes))

    LS = "--"

    ######################################################################################################
    ######################################################################################################

    #########################################################
    ###### varying texture height -- density
    figure(figsize = (5.2, 5.6))
    fontsize = 25

    semilogy(YY, nnn1itp(XX, YY), linewidth = 5, color = Blues(241))
    semilogy(YY, nn1itp(XX, YY), linewidth = 5, color = Blues(241))
    semilogy(YY, nnp1itp(XX, YY), linewidth = 5, color = Blues(241))

    semilogy(YY, nnn2itp(XX, YY), linewidth = 5, color = Blues(201))
    semilogy(YY, nn2itp(XX, YY), linewidth = 5, color = Blues(201))
    semilogy(YY, nnp2itp(XX, YY), linewidth = 5, color = Blues(201))

    semilogy(YY, nnn3itp(XX, YY), linewidth = 5, color = Blues(161))
    semilogy(YY, nn3itp(XX, YY), linewidth = 5, color = Blues(161))
    semilogy(YY, nnp3itp(XX, YY), linewidth = 5, color = Blues(161))

    ################
    semilogy(YY, npn1itp(XX, YY), linewidth = 5, color = Oranges(241))
    semilogy(YY, np1itp(XX, YY), linewidth = 5, color = Oranges(241))
    semilogy(YY, npp1itp(XX, YY), linewidth = 5, color = Oranges(241))

    semilogy(YY, npn2itp(XX, YY), linewidth = 5, color = Oranges(201))
    semilogy(YY, np2itp(XX, YY), linewidth = 5, color = Oranges(201))
    semilogy(YY, npp2itp(XX, YY), linewidth = 5, color = Oranges(201))

    semilogy(YY, npn3itp(XX, YY), linewidth = 5, color = Oranges(161))
    semilogy(YY, np3itp(XX, YY), linewidth = 5, color = Oranges(161))
    semilogy(YY, npp3itp(XX, YY), linewidth = 5, color = Oranges(161))

    PyPlot.xlim(0, 440)
    PyPlot.xticks([0, 200, 400])
    PyPlot.xlabel("\$y\$ [nm]", fontsize = fontsize)
    PyPlot.ylim(1.0e14, 1.0e25)
    PyPlot.ylabel("Density [\$\\mathrm{m}^{-3}\$]", fontsize = fontsize)
    PyPlot.title("x = $XVal nm")
    PyPlot.tick_params(which = "both", labelsize = fontsize)
    PyPlot.tight_layout()

    if saveFig
        savefig(datadir("2D-density-vertical-line-x-$XVal-generation-Maxwell-forw-$V.pdf"))
    end

    return nothing
end

end
