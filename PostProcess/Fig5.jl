#=

Code for visualizing carrier densities, band diagram and VOC, injection barrier study

=#

module Fig5

using PyPlot
using DelimitedFiles
using ChargeTransport
using ExtendableGrids
using TexturedPerovskiteSolarCells

# for convenience
datadir = TexturedPerovskiteSolarCells.datadir
scriptsdir = TexturedPerovskiteSolarCells.scriptsdir

include(scriptsdir("SingleJunction.jl"))

function main(;scanrate         = 1000.0,  # "10p0" # "0p001"
              typeReco          = "all",   # "radiative"
              generationUniform = false,
              IVDirection       = "forw",  # "rev" #
              V                 = "end",   # "inival", # "V-1p15"
              printText = true, saveFig = false,
              parameter_file = scriptsdir("params_single_junction.jl")
              )

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
    ########
    grid2, ctsys2 = SingleJunction.main(gridDim = 2, typeGrid = "nanotextured", amplitude = amplitude, generation = true, generationUniform = generationUniform, plotPostProcess = true, printText = false, plotting = false)
    subg2        = subgrid(grid2, [regionPero]); subgn2 = subgrid(grid2, [regionETL1]); subgp2 = subgrid(grid2, [regionHTL])
    ########
    grid3, ctsys3 = SingleJunction.main(gridDim = 2, typeGrid = "nanotextured", amplitude = amplitude2, generation = true, generationUniform = generationUniform, plotPostProcess = true, printText = false, plotting = false)
    subg3        = subgrid(grid3, [regionPero]); subgn3 = subgrid(grid3, [regionETL1]); subgp3 = subgrid(grid3, [regionHTL])
    ########
    ########
    grid4, ctsys4 = SingleJunction.main(gridDim = 2, typeGrid = "nanotextured", amplitude = amplitude3, generation = true, generationUniform = generationUniform, plotPostProcess = true, printText = false, plotting = false)
    subg4         = subgrid(grid4, [regionPero]); subgn4 = subgrid(grid4, [regionETL1]); subgp4 = subgrid(grid4, [regionHTL])
    #################################################################
    #################################################################

    coord1  = subg1[Coordinates]; coord2 = subg2[Coordinates]
    coord3  = subg3[Coordinates]; coord4 = subg4[Coordinates]

    coord1  = coord1./nm; coord2 = coord2./nm
    coord3  = coord3./nm; coord4 = coord4./nm

    triang1 = matplotlib.tri.Triangulation(coord1[1,:], coord1[2, :], transpose(Matrix(subg1[CellNodes]) .- 1))
    triang2 = matplotlib.tri.Triangulation(coord2[1,:], coord2[2, :], transpose(Matrix(subg2[CellNodes]) .- 1))
    triang3 = matplotlib.tri.Triangulation(coord3[1,:], coord3[2, :], transpose(Matrix(subg3[CellNodes]) .- 1))
    triang4 = matplotlib.tri.Triangulation(coord4[1,:], coord4[2, :], transpose(Matrix(subg4[CellNodes]) .- 1))

    ###########################################################
    ###########################################################

    coordn1  = subgn1[Coordinates]; coordn2 = subgn2[Coordinates]
    coordn3  = subgn3[Coordinates]; coordn4 = subgn4[Coordinates]

    coordn1  = coordn1./nm; coordn2 = coordn2./nm
    coordn3  = coordn3./nm; coordn4 = coordn4./nm

    triangn1 = matplotlib.tri.Triangulation(coordn1[1,:], coordn1[2, :], transpose(Matrix(subgn1[CellNodes]) .- 1))
    triangn2 = matplotlib.tri.Triangulation(coordn2[1,:], coordn2[2, :], transpose(Matrix(subgn2[CellNodes]) .- 1))
    triangn3 = matplotlib.tri.Triangulation(coordn3[1,:], coordn3[2, :], transpose(Matrix(subgn3[CellNodes]) .- 1))
    triangn4 = matplotlib.tri.Triangulation(coordn4[1,:], coordn4[2, :], transpose(Matrix(subgn4[CellNodes]) .- 1))

    ###########################################################
    ###########################################################

    coordp1  = subgp1[Coordinates]; coordp2 = subgp2[Coordinates]
    coordp3  = subgp3[Coordinates]; coordp4 = subgp4[Coordinates]

    coordp1  = coordp1./nm; coordp2 = coordp2./nm
    coordp3  = coordp3./nm; coordp4 = coordp4./nm

    triangp1 = matplotlib.tri.Triangulation(coordp1[1,:], coordp1[2, :], transpose(Matrix(subgp1[CellNodes]) .- 1))
    triangp2 = matplotlib.tri.Triangulation(coordp2[1,:], coordp2[2, :], transpose(Matrix(subgp2[CellNodes]) .- 1))
    triangp3 = matplotlib.tri.Triangulation(coordp3[1,:], coordp3[2, :], transpose(Matrix(subgp3[CellNodes]) .- 1))
    triangp4 = matplotlib.tri.Triangulation(coordp4[1,:], coordp4[2, :], transpose(Matrix(subgp4[CellNodes]) .- 1))

    #################################################################
    ## Plotting
    #################################################################

    subg1[Coordinates] = subg1[Coordinates]./nm; subg2[Coordinates] = subg2[Coordinates]./nm
    subg3[Coordinates] = subg3[Coordinates]./nm; subg4[Coordinates] = subg4[Coordinates]./nm

    #################################################################
    ## read in solution
    #################################################################

    sol1    = readdlm(datadir("sol", "$pathSol/Sol-2D-$IVDirection-planar-generation-$generation-reco-$typeReco-$V.dat"))'
    sol2    = readdlm(datadir("sol", "$pathSol/Sol-2D-$IVDirection-nanotextured-ampl-$textampl-generation-$generation-reco-$typeReco-$V.dat"))'
    sol3    = readdlm(datadir("sol", "$pathSol/Sol-2D-$IVDirection-nanotextured-ampl-$textampl2-generation-$generation-reco-$typeReco-$V.dat"))'
    sol4    = readdlm(datadir("sol", "$pathSol/Sol-2D-$IVDirection-nanotextured-ampl-$textampl3-generation-$generation-reco-$typeReco-$V.dat"))'

    psi1      = view(sol1[ipsi, :], subg1); psi2 = view(sol2[ipsi, :], subg2)
    psi3      = view(sol3[ipsi, :], subg3); psi4 = view(sol4[ipsi, :], subg4)

    psi1itp   = matplotlib.tri.LinearTriInterpolator(triang1, psi1)
    psi2itp   = matplotlib.tri.LinearTriInterpolator(triang2, psi2)
    psi3itp   = matplotlib.tri.LinearTriInterpolator(triang3, psi3)
    psi4itp   = matplotlib.tri.LinearTriInterpolator(triang4, psi4)

    psin1     = view(sol1[ipsi, :], subgn1); psin2 = view(sol2[ipsi, :], subgn2)
    psin3     = view(sol3[ipsi, :], subgn3); psin4 = view(sol4[ipsi, :], subgn4)

    psin1itp  = matplotlib.tri.LinearTriInterpolator(triangn1, psin1)
    psin2itp  = matplotlib.tri.LinearTriInterpolator(triangn2, psin2)
    psin3itp  = matplotlib.tri.LinearTriInterpolator(triangn3, psin3)
    psin4itp  = matplotlib.tri.LinearTriInterpolator(triangn4, psin4)

    psip1     = view(sol1[ipsi, :], subgp1); psip2 = view(sol2[ipsi, :], subgp2)
    psip3     = view(sol3[ipsi, :], subgp3); psip4 = view(sol4[ipsi, :], subgp4)

    psip1itp  = matplotlib.tri.LinearTriInterpolator(triangp1, psip1)
    psip2itp  = matplotlib.tri.LinearTriInterpolator(triangp2, psip2)
    psip3itp  = matplotlib.tri.LinearTriInterpolator(triangp3, psip3)
    psip4itp  = matplotlib.tri.LinearTriInterpolator(triangp4, psip4)

    ###########################################
    phin1     = view(sol1[iphin, :], subg1); phin2 = view(sol2[iphin, :], subg2)
    phin3     = view(sol3[iphin, :], subg3); phin4 = view(sol4[iphin, :], subg4)

    phip1     = view(sol1[iphip, :], subg1); phip2 = view(sol2[iphip, :], subg2)
    phip3     = view(sol3[iphip, :], subg3); phip4 = view(sol4[iphip, :], subg4)
    ###
    phinn1    = view(sol1[iphin, :], subgn1); phinn2 = view(sol2[iphin, :], subgn2)
    phinn3    = view(sol3[iphin, :], subgn3); phinn4 = view(sol4[iphin, :], subgn4)

    phipn1    = view(sol1[iphip, :], subgn1); phipn2 = view(sol2[iphip, :], subgn2)
    phipn3    = view(sol3[iphip, :], subgn3); phipn4 = view(sol4[iphip, :], subgn4)
    ###
    phinp1    = view(sol1[iphin, :], subgp1); phinp2 = view(sol2[iphin, :], subgp2)
    phinp3    = view(sol3[iphin, :], subgp3); phinp4 = view(sol4[iphin, :], subgp4)

    phipp1    = view(sol1[iphip, :], subgp1); phipp2 = view(sol2[iphip, :], subgp2)
    phipp3    = view(sol3[iphip, :], subgp3); phipp4 = view(sol4[iphip, :], subgp4)

    phin1itp  = matplotlib.tri.LinearTriInterpolator(triang1, phin1)
    phin2itp  = matplotlib.tri.LinearTriInterpolator(triang2, phin2)
    phin3itp  = matplotlib.tri.LinearTriInterpolator(triang3, phin3)
    phin4itp  = matplotlib.tri.LinearTriInterpolator(triang4, phin4)

    phip1itp  = matplotlib.tri.LinearTriInterpolator(triang1, phip1)
    phip2itp  = matplotlib.tri.LinearTriInterpolator(triang2, phip2)
    phip3itp  = matplotlib.tri.LinearTriInterpolator(triang3, phip3)
    phip4itp  = matplotlib.tri.LinearTriInterpolator(triang4, phip4)
    ###
    phinn1itp = matplotlib.tri.LinearTriInterpolator(triangn1, phinn1)
    phinn2itp = matplotlib.tri.LinearTriInterpolator(triangn2, phinn2)
    phinn3itp = matplotlib.tri.LinearTriInterpolator(triangn3, phinn3)
    phinn4itp = matplotlib.tri.LinearTriInterpolator(triangn4, phinn4)

    phipn1itp = matplotlib.tri.LinearTriInterpolator(triangn1, phipn1)
    ###
    phinp1itp = matplotlib.tri.LinearTriInterpolator(triangp1, phinp1)

    phipp1itp = matplotlib.tri.LinearTriInterpolator(triangp1, phipp1)
    phipp2itp = matplotlib.tri.LinearTriInterpolator(triangp2, phipp2)
    phipp3itp = matplotlib.tri.LinearTriInterpolator(triangp3, phipp3)
    phipp4itp = matplotlib.tri.LinearTriInterpolator(triangp4, phipp4)

    ###########################################

    nn1     = Nn[regionPero] .* Fcc[iphin].(zn*( q*(view(sol1[iphin, :], subg1).-view(sol1[ipsi, :], subg1)) .+ En[regionPero])./(kB*T))
    nn2     = Nn[regionPero] .* Fcc[iphin].(zn*( q*(view(sol2[iphin, :], subg2).-view(sol2[ipsi, :], subg2)) .+ En[regionPero])./(kB*T))
    nn3     = Nn[regionPero] .* Fcc[iphin].(zn*( q*(view(sol3[iphin, :], subg3).-view(sol3[ipsi, :], subg3)) .+ En[regionPero])./(kB*T))
    nn4     = Nn[regionPero] .* Fcc[iphin].(zn*( q*(view(sol4[iphin, :], subg4).-view(sol4[ipsi, :], subg4)) .+ En[regionPero])./(kB*T))

    np1     = Np[regionPero] .* Fcc[iphip].(zp*( q*(view(sol1[iphip, :], subg1).-view(sol1[ipsi, :], subg1)) .+ Ep[regionPero])./(kB*T))
    np2     = Np[regionPero] .* Fcc[iphip].(zp*( q*(view(sol2[iphip, :], subg2).-view(sol2[ipsi, :], subg2)) .+ Ep[regionPero])./(kB*T))
    np3     = Np[regionPero] .* Fcc[iphip].(zp*( q*(view(sol3[iphip, :], subg3).-view(sol3[ipsi, :], subg3)) .+ Ep[regionPero])./(kB*T))
    np4     = Np[regionPero] .* Fcc[iphip].(zp*( q*(view(sol4[iphip, :], subg4).-view(sol4[ipsi, :], subg4)) .+ Ep[regionPero])./(kB*T))

    ######
    nn1itp  = matplotlib.tri.LinearTriInterpolator(triang1, nn1)
    nn2itp  = matplotlib.tri.LinearTriInterpolator(triang2, nn2)
    nn3itp  = matplotlib.tri.LinearTriInterpolator(triang3, nn3)
    nn4itp  = matplotlib.tri.LinearTriInterpolator(triang4, nn4)

    np1itp  = matplotlib.tri.LinearTriInterpolator(triang1, np1)
    np2itp  = matplotlib.tri.LinearTriInterpolator(triang2, np2)
    np3itp  = matplotlib.tri.LinearTriInterpolator(triang3, np3)
    np4itp  = matplotlib.tri.LinearTriInterpolator(triang4, np4)

    ##########################################
    nnn1    = Nn[regionETL1] .* Fcc[iphin].(zn*( q*(view(sol1[iphin, :], subgn1).-view(sol1[ipsi, :], subgn1)) .+ En[regionETL1])./(kB*T))
    nnn2    = Nn[regionETL1] .* Fcc[iphin].(zn*( q*(view(sol2[iphin, :], subgn2).-view(sol2[ipsi, :], subgn2)) .+ En[regionETL1])./(kB*T))
    nnn3    = Nn[regionETL1] .* Fcc[iphin].(zn*( q*(view(sol3[iphin, :], subgn3).-view(sol3[ipsi, :], subgn3)) .+ En[regionETL1])./(kB*T))
    nnn4    = Nn[regionETL1] .* Fcc[iphin].(zn*( q*(view(sol4[iphin, :], subgn4).-view(sol4[ipsi, :], subgn4)) .+ En[regionETL1])./(kB*T))

    npn1    = Np[regionETL1] .* Fcc[iphip].(zp*( q*(view(sol1[iphip, :], subgn1).-view(sol1[ipsi, :], subgn1)) .+ Ep[regionETL1])./(kB*T))
    npn2    = Np[regionETL1] .* Fcc[iphip].(zp*( q*(view(sol2[iphip, :], subgn2).-view(sol2[ipsi, :], subgn2)) .+ Ep[regionETL1])./(kB*T))
    npn3    = Np[regionETL1] .* Fcc[iphip].(zp*( q*(view(sol3[iphip, :], subgn3).-view(sol3[ipsi, :], subgn3)) .+ Ep[regionETL1])./(kB*T))
    npn4    = Np[regionETL1] .* Fcc[iphip].(zp*( q*(view(sol4[iphip, :], subgn4).-view(sol4[ipsi, :], subgn4)) .+ Ep[regionETL1])./(kB*T))

    ######
    nnn1itp = matplotlib.tri.LinearTriInterpolator(triangn1, nnn1)
    nnn2itp = matplotlib.tri.LinearTriInterpolator(triangn2, nnn2)
    nnn3itp = matplotlib.tri.LinearTriInterpolator(triangn3, nnn3)
    nnn4itp = matplotlib.tri.LinearTriInterpolator(triangn4, nnn4)

    npn1itp = matplotlib.tri.LinearTriInterpolator(triangn1, npn1)
    npn2itp = matplotlib.tri.LinearTriInterpolator(triangn2, npn2)
    npn3itp = matplotlib.tri.LinearTriInterpolator(triangn3, npn3)
    npn4itp = matplotlib.tri.LinearTriInterpolator(triangn4, npn4)

    ##########################################

    nnp1    = Nn[regionHTL] .* Fcc[iphin].(zn*( q*(view(sol1[iphin, :], subgp1).-view(sol1[ipsi, :], subgp1)) .+ En[regionHTL])./(kB*T))
    nnp2    = Nn[regionHTL] .* Fcc[iphin].(zn*( q*(view(sol2[iphin, :], subgp2).-view(sol2[ipsi, :], subgp2)) .+ En[regionHTL])./(kB*T))
    nnp3    = Nn[regionHTL] .* Fcc[iphin].(zn*( q*(view(sol3[iphin, :], subgp3).-view(sol3[ipsi, :], subgp3)) .+ En[regionHTL])./(kB*T))
    nnp4    = Nn[regionHTL] .* Fcc[iphin].(zn*( q*(view(sol4[iphin, :], subgp4).-view(sol4[ipsi, :], subgp4)) .+ En[regionHTL])./(kB*T))

    npp1    = Np[regionHTL] .* Fcc[iphip].(zp*( q*(view(sol1[iphip, :], subgp1).-view(sol1[ipsi, :], subgp1)) .+ Ep[regionHTL])./(kB*T))
    npp2    = Np[regionHTL] .* Fcc[iphip].(zp*( q*(view(sol2[iphip, :], subgp2).-view(sol2[ipsi, :], subgp2)) .+ Ep[regionHTL])./(kB*T))
    npp3    = Np[regionHTL] .* Fcc[iphip].(zp*( q*(view(sol3[iphip, :], subgp3).-view(sol3[ipsi, :], subgp3)) .+ Ep[regionHTL])./(kB*T))
    npp4    = Np[regionHTL] .* Fcc[iphip].(zp*( q*(view(sol4[iphip, :], subgp4).-view(sol4[ipsi, :], subgp4)) .+ Ep[regionHTL])./(kB*T))

    ######
    nnp1itp = matplotlib.tri.LinearTriInterpolator(triangp1, nnp1)
    nnp2itp = matplotlib.tri.LinearTriInterpolator(triangp2, nnp2)
    nnp3itp = matplotlib.tri.LinearTriInterpolator(triangp3, nnp3)
    nnp4itp = matplotlib.tri.LinearTriInterpolator(triangp4, nnp4)

    npp1itp = matplotlib.tri.LinearTriInterpolator(triangp1, npp1)
    npp2itp = matplotlib.tri.LinearTriInterpolator(triangp2, npp2)
    npp3itp = matplotlib.tri.LinearTriInterpolator(triangp3, npp3)
    npp4itp = matplotlib.tri.LinearTriInterpolator(triangp4, npp4)

    #################################################################
    ## Plotting
    #################################################################

    fCos(x, ampl) = 0.5*ampl .* cos.(pi .+ 2 .* pi .* x./ 750.0) .+ 0.5*ampl .+ 30.0 + 400 - ampl/2

    Blues   = get_cmap(:Blues); Oranges = get_cmap(:Oranges)

    XVal  = 187.0; nodes = 600
    XX    = XVal .* ones(nodes)
    YY    = collect(range(0.0, 450.0, length = nodes))

    LS = "--"

    ######################################################################################################
    ######################################################################################################

    #########################################################
    ###### varying texture height -- density
    semilogy(YY, nnn1itp(XX, YY), linewidth = 5, color = Blues(241))
    semilogy(YY, nn1itp(XX, YY),  linewidth = 5, color = Blues(241))
    semilogy(YY, nnp1itp(XX, YY), linewidth = 5, color = Blues(241))

    semilogy(YY, nnn2itp(XX, YY), linewidth = 5, color = Blues(201))
    semilogy(YY, nn2itp(XX, YY),  linewidth = 5, color = Blues(201))
    semilogy(YY, nnp2itp(XX, YY), linewidth = 5, color = Blues(201))

    semilogy(YY, nnn3itp(XX, YY), linewidth = 5, color = Blues(161))
    semilogy(YY, nn3itp(XX, YY),  linewidth = 5, color = Blues(161))
    semilogy(YY, nnp3itp(XX, YY), linewidth = 5, color = Blues(161))

    semilogy(YY, nnn4itp(XX, YY), linewidth = 5, color = Blues(121))
    semilogy(YY, nn4itp(XX, YY),  linewidth = 5, color = Blues(121))
    semilogy(YY, nnp4itp(XX, YY), linewidth = 5, color = Blues(121))
    ################
    semilogy(YY, npn1itp(XX, YY), linewidth = 5, color = Oranges(241))
    semilogy(YY, np1itp(XX, YY),  linewidth = 5, color = Oranges(241))
    semilogy(YY, npp1itp(XX, YY), linewidth = 5, color = Oranges(241))

    semilogy(YY, npn2itp(XX, YY), linewidth = 5, color = Oranges(201))
    semilogy(YY, np2itp(XX, YY),  linewidth = 5, color = Oranges(201))
    semilogy(YY, npp2itp(XX, YY), linewidth = 5, color = Oranges(201))

    semilogy(YY, npn3itp(XX, YY), linewidth = 5, color = Oranges(161))
    semilogy(YY, np3itp(XX, YY),  linewidth = 5, color = Oranges(161))
    semilogy(YY, npp3itp(XX, YY), linewidth = 5, color = Oranges(161))

    semilogy(YY, npn4itp(XX, YY), linewidth = 5, color = Oranges(121))
    semilogy(YY, np4itp(XX, YY),  linewidth = 5, color = Oranges(121))
    semilogy(YY, npp4itp(XX, YY), linewidth = 5, color = Oranges(121))

    PyPlot.xlim(0, 440)
    PyPlot.xticks([0, 200, 400])
    PyPlot.xlabel("\$y\$ [nm]", fontsize=18)
    PyPlot.ylim(1.0e14, 1.0e25)
    PyPlot.ylabel("Density [\$\\frac{1}{\\mathrm{m}^3}\$]", fontsize=18)
    PyPlot.title("x = $XVal nm")
    PyPlot.tick_params(which ="both", labelsize=18)
    PyPlot.tight_layout()

    if saveFig
        savefig(datadir("2D-density-vertical-line-x-$XVal-scanrate-$scanrate-generation-$generation-$IVDirection-$V.pdf"))
    end

    ###### planar plot -- all energies
    figure()
    plot(YY, En[regionETL1]./q .- psin1itp(XX, YY), linewidth = 5, color = Blues(241))
    plot(YY, En[regionHTL]./q  .- psip1itp(XX, YY), linewidth = 5, color = Blues(241))
    plot(YY, En[regionPero]./q .- psi1itp(XX, YY),  linewidth = 5, color = Blues(241))
    ###
    plot(YY, Ep[regionETL1]./q .- psin1itp(XX, YY), linewidth = 5, color = Oranges(241))
    plot(YY, Ep[regionHTL]./q  .- psip1itp(XX, YY), linewidth = 5, color = Oranges(241))
    plot(YY, Ep[regionPero]./q .- psi1itp(XX, YY),  linewidth = 5, color = Oranges(241))

    plot(YY, -phinn1itp(XX, YY),                    linewidth = 5, linestyle = LS, color = Blues(241))
    plot(YY, -phinp1itp(XX, YY),                    linewidth = 5, linestyle = LS, color = Blues(241))
    plot(YY, -phin1itp(XX, YY),                     linewidth = 5, linestyle = LS, color = Blues(241))
    ###
    plot(YY, -phipn1itp(XX, YY),                    linewidth = 5, linestyle = LS, color = Oranges(201))
    plot(YY, -phipp1itp(XX, YY),                    linewidth = 5, linestyle = LS, color = Oranges(201))
    plot(YY, -phip1itp(XX, YY),                     linewidth = 5, linestyle = LS, color = Oranges(201))

    PyPlot.xlim(0, 440)
    PyPlot.xticks([0, 200, 400])
    PyPlot.xlabel("\$y \$ [nm]", fontsize=18)
    PyPlot.ylim(-3.0, 2.0)
    PyPlot.yticks([-2, -1, 0, 1])
    PyPlot.ylabel("Energy [eV] ", fontsize=18)
    PyPlot.title("x-$XVal")
    PyPlot.tick_params(which ="both", labelsize=18)
    PyPlot.tight_layout()

    if saveFig
        savefig(datadir("2D-band-diagram-planar-vertical-line-x-$XVal-scanrate-$scanrate-generation-$generation-$IVDirection-$V.pdf"))
    end

    #########################################################

    ###### varying texture height -- Ec
    figure()
    plot(YY, En[regionETL1]./q .- psin1itp(XX, YY), linewidth = 5, color = Blues(241))
    plot(YY, En[regionPero]./q .- psi1itp(XX, YY),  linewidth = 5, color = Blues(241))
    plot(YY, En[regionHTL]./q  .- psip1itp(XX, YY), linewidth = 5, color = Blues(241))

    plot(YY, En[regionETL1]./q .- psin2itp(XX, YY), linewidth = 5, color = Blues(201))
    plot(YY, En[regionPero]./q .- psi2itp(XX, YY),  linewidth = 5, color = Blues(201))
    plot(YY, En[regionHTL]./q  .- psip2itp(XX, YY), linewidth = 5, color = Blues(201))

    plot(YY, En[regionETL1]./q .- psin3itp(XX, YY), linewidth = 5, color = Blues(161))
    plot(YY, En[regionPero]./q .- psi3itp(XX, YY),  linewidth = 5, color = Blues(161))
    plot(YY, En[regionHTL]./q  .- psip3itp(XX, YY), linewidth = 5, color = Blues(161))

    plot(YY, En[regionETL1]./q .- psin4itp(XX, YY), linewidth = 5, color = Blues(121))
    plot(YY, En[regionPero]./q .- psi4itp(XX, YY),  linewidth = 5, color = Blues(121))
    plot(YY, En[regionHTL]./q  .- psip4itp(XX, YY), linewidth = 5, color = Blues(121))

    PyPlot.xlim(0, 430)
    PyPlot.xticks([0, 200, 400])
    PyPlot.xlabel("\$y\$ [nm]", fontsize=18)
    PyPlot.ylim(0.07, 0.3)
    PyPlot.yticks([0.1, 0.2, 0.3])
    PyPlot.ylabel("Energy [eV] ", fontsize=18)
    PyPlot.title("x-$XVal")
    PyPlot.tick_params(which ="both", labelsize=18)
    PyPlot.tight_layout()

    if saveFig
        savefig(datadir("2D-Ec-vertical-line-x-$XVal-scanrate-$scanrate-generation-$generation-$IVDirection-$V.pdf"))
    end

    #########################################################
    ###### varying texture height -- EFn
    figure()
    plot(YY, -phinn1itp(XX, YY), linewidth = 5, linestyle = LS, color = Blues(241))
    plot(YY, -phinn2itp(XX, YY), linewidth = 5, linestyle = LS, color = Blues(201))
    plot(YY, -phinn3itp(XX, YY), linewidth = 5, linestyle = LS, color = Blues(161))
    plot(YY, -phinn4itp(XX, YY), linewidth = 5, linestyle = LS, color = Blues(121))

    plot(YY, -phin1itp(XX, YY), linewidth = 5, linestyle = LS, color = Blues(241))
    plot(YY, -phin2itp(XX, YY), linewidth = 5, linestyle = LS, color = Blues(201))
    plot(YY, -phin3itp(XX, YY), linewidth = 5, linestyle = LS, color = Blues(161))
    plot(YY, -phin4itp(XX, YY), linewidth = 5, linestyle = LS, color = Blues(121))

    PyPlot.xlim(0, 430)
    PyPlot.xticks([0, 200, 400])
    PyPlot.xlabel("\$y\$ [nm]", fontsize=18)
    PyPlot.ylim(-0.002, 0.005)
    PyPlot.yticks([0.0, 0.002])
    PyPlot.ylabel("Energy [eV] ", fontsize=18)
    PyPlot.title("x-$XVal")
    PyPlot.tick_params(which ="both", labelsize=18)
    PyPlot.tight_layout()

    if saveFig
        savefig(datadir("2D-EFn-vertical-line-x-$XVal-scanrate-$scanrate-generation-$generation-$IVDirection-$V.pdf"))
    end

    #########################################################
    ###### varying texture height -- EFp
    figure()
    plot(YY, -phipp1itp(XX, YY),  linestyle = LS, linewidth = 5, color = Oranges(241))
    plot(YY, -phipp2itp(XX, YY),  linestyle = LS, linewidth = 5, color = Oranges(201))
    plot(YY, -phipp3itp(XX, YY),  linestyle = LS, linewidth = 5, color = Oranges(161))
    plot(YY, -phipp4itp(XX, YY),  linestyle = LS, linewidth = 5, color = Oranges(121))

    plot(YY, -phip1itp(XX, YY),  linestyle = LS, linewidth = 5, color = Oranges(241))
    plot(YY, -phip2itp(XX, YY),  linestyle = LS, linewidth = 5, color = Oranges(201))
    plot(YY, -phip3itp(XX, YY),  linestyle = LS, linewidth = 5, color = Oranges(161))
    plot(YY, -phip4itp(XX, YY),  linestyle = LS, linewidth = 5, color = Oranges(121))

    PyPlot.xlim(30, 440)
    PyPlot.xticks([200, 400])
    PyPlot.xlabel("\$y\$ [nm]", fontsize=18)
    PyPlot.ylim(-1.202, -1.195)
    PyPlot.yticks([-1.2, -1.197])
    PyPlot.ylabel("Energy [eV] ", fontsize=18)
    PyPlot.title("x-$XVal")
    PyPlot.tick_params(which ="both", labelsize=18)
    PyPlot.tight_layout()

    if saveFig
        savefig(datadir("2D-EFp-vertical-line-x-$XVal-scanrate-$scanrate-generation-$generation-$IVDirection-$V.pdf"))
    end

    #########################################################
    ###### varying texture height -- Ev
    figure()
    plot(YY, Ep[regionETL1]./q .- psin1itp(XX, YY), linewidth = 5, color = Oranges(241))
    plot(YY, Ep[regionPero]./q .- psi1itp(XX, YY),  linewidth = 5, color = Oranges(241))
    plot(YY, Ep[regionHTL]./q  .- psip1itp(XX, YY), linewidth = 5, color = Oranges(241))

    plot(YY, Ep[regionETL1]./q .- psin2itp(XX, YY), linewidth = 5, color = Oranges(201))
    plot(YY, Ep[regionPero]./q .- psi2itp(XX, YY),  linewidth = 5, color = Oranges(201))
    plot(YY, Ep[regionHTL]./q  .- psip2itp(XX, YY), linewidth = 5, color = Oranges(201))

    plot(YY, Ep[regionETL1]./q .- psin3itp(XX, YY), linewidth = 5, color = Oranges(161))
    plot(YY, Ep[regionPero]./q .- psi3itp(XX, YY),  linewidth = 5, color = Oranges(161))
    plot(YY, Ep[regionHTL]./q  .- psip3itp(XX, YY), linewidth = 5, color = Oranges(161))

    plot(YY, Ep[regionETL1]./q .- psin4itp(XX, YY), linewidth = 5, color = Oranges(121))
    plot(YY, Ep[regionPero]./q .- psi4itp(XX, YY),  linewidth = 5, color = Oranges(121))
    plot(YY, Ep[regionHTL]./q  .- psip4itp(XX, YY), linewidth = 5, color = Oranges(121))

    PyPlot.xlim(30, 440)
    PyPlot.xticks([200, 400])
    PyPlot.xlabel("\$y\$ [nm]", fontsize=18)
    PyPlot.ylim(-1.52, -1.29)
    PyPlot.yticks([-1.5, -1.4])
    PyPlot.ylabel("Energy [eV] ", fontsize=18)
    PyPlot.title("x-$XVal")
    PyPlot.tick_params(which ="both", labelsize=18)
    PyPlot.tight_layout()

    if saveFig
        savefig(datadir("2D-Ep-vertical-line-x-$XVal-scanrate-$scanrate-generation-$generation-$IVDirection-$V.pdf"))
    end

    #########################################################
    ###### varying texture height -- Injection barriers
    figure()
    MKsize = 10
    Diffn  = readdlm(datadir("Ec-Efn-difference-params-$paramsname-generation-$generation-scanrate-$textSR.dat"))
    Diffp  = readdlm(datadir("Ev-Efp-difference-params-$paramsname-generation-$generation-scanrate-$textSR.dat"))

    Coln = [143/255, 187/255, 217/255]
    Colp = [234/255, 147/255, 147/255]

    plot(Diffn[:, 1], Diffn[:, 2], marker ="o", markersize = MKsize, color = Coln, label = "\$  | \\Phi_{\\mathrm{n}}  | \$")
    plot(Diffp[:, 1], Diffp[:, 2], marker ="o", markersize = MKsize, color = Colp, label = "\$   | \\Phi_{\\mathrm{p}} | \$")

    PyPlot.xlim(-20.0, 770)
    PyPlot.xlabel("Texture height [nm]", fontsize=18)
    PyPlot.ylim(0.15, 0.3)
    PyPlot.ylabel(" Energy [eV]", fontsize=18)
    PyPlot.yticks([0.2, 0.3])
    PyPlot.tick_params(which ="both", labelsize=18)
    PyPlot.legend(fontsize = 18)
    PyPlot.tight_layout()

    if saveFig
        savefig(datadir("Injection-barrier-params-$paramsname-generation-$generation-scanrate-$scanrate.pdf"))
    end

    #########################################################
    ###### varying texture height -- VOC & QFL
    figure()
    MKsize = 10
    VocIV  = readdlm(datadir("VOC-params-$paramsname-generation-$generation-scanrate-$textSR.dat"))
    QFLS   = readdlm(datadir("Ef-difference-params-$paramsname-generation-$generation-scanrate-$textSR.dat"))

    ColVoc = [129/255, 120/255, 213/255]
    ColQFL = [172/255, 235/255, 180/255]

    plot(VocIV[:, 1], VocIV[:, 2], marker ="o",  markersize = MKsize, color = ColVoc, label = "\$ V_{\\mathrm{OC}} \$")
    plot(QFLS[:, 1],  QFLS[:, 2],  marker ="o",  markersize = MKsize, color = ColQFL, label = "\$ \\Delta E_{\\mathrm{F}} \$")

    PyPlot.xlim(-20.0, 770)
    PyPlot.ylim(1.165, 1.2)
    PyPlot.yticks([1.17, 1.18, 1.19])
    PyPlot.xlabel("Texture height [nm]", fontsize=18)
    PyPlot.ylabel("Voltage [V]", fontsize=18)
    PyPlot.tick_params(which ="both", labelsize=18)
    PyPlot.legend(fontsize = 18)
    PyPlot.tight_layout()

    if saveFig
        savefig(datadir("Energy-difference-params-$paramsname-generation-$generation-scanrate-$scanrate.pdf"))
    end

    return nothing
end

end