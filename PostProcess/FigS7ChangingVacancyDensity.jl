#=

Code for visualizing impact of vacancy density on electrons and holes

=#
module FigS7ChangingVacancyDensity

using PyPlot
using DelimitedFiles
using ChargeTransport
using ExtendableGrids
using TexturedPerovskiteSolarCells

# for convenience
datadir = TexturedPerovskiteSolarCells.datadir
scriptsdir = TexturedPerovskiteSolarCells.scriptsdir

include(scriptsdir("SingleJunction.jl"))

function main(;scanrate   = 1000.0,   # "10p0" # "0p001"
              typeReco    = "all",    # "radiative"
              IVDirection = "forw", #"rev", #
              V           = "end", #"inival", # "V-1p15"
              printText = true, saveFig = false,
              parameter_file = scriptsdir("params_single_junction.jl")
              )

    include(parameter_file)

    PyPlot.rc("font", family="sans-serif", size=14)
    PyPlot.rc("mathtext", fontset="dejavusans")
    PyPlot.close("all")

    helpSR  = collect(string(scanrate));  helpSR[ findall(x -> x == '.', helpSR)[1] ] = 'p'
    textSR  = join(helpSR)

    pathSol = "parameter-$paramsname/$textSR"

    #################################################################
    ## grid
    #################################################################

    grid, ctsys = SingleJunction.main(gridDim = 1, typeGrid = "planar", generation = true, plotPostProcess = true, printText = false, plotting = false)
    subg        = subgrid(grid, [regionPero]); subgn = subgrid(grid, [regionETL1]); subgp = subgrid(grid, [regionHTL])
    data        = ctsys.fvmsys.physics.data

    coord       = subg[Coordinates];  coord    = coord./nm
    coordn      = subgn[Coordinates]; coordn   = coordn./nm
    coordp      = subgp[Coordinates]; coordp   = coordp./nm
    coordALL    = grid[Coordinates];  coordALL = coordALL./nm

    #########################################################################################################
    #########################################################################################################
    #########################################################################################################

    Ca1 = 1.0e21 / (m^3); Ea1 = -5.614 * eV
    Ca2 = 1.0e22 / (m^3); Ea2 = -5.532 * eV
    Ca4 = 1.0e23 / (m^3); Ea4 = -5.283 * eV

    sol1 = readdlm(datadir("sol", "$pathSol/Sol-1D-$IVDirection-planar-Ca-1p0e21-generation-Maxwell-reco-$typeReco-$V.dat"))'
    sol2 = readdlm(datadir("sol", "$pathSol/Sol-1D-$IVDirection-planar-Ca-1p0e22-generation-Maxwell-reco-$typeReco-$V.dat"))'
    sol4 = readdlm(datadir("sol", "$pathSol/Sol-1D-$IVDirection-planar-Ca-1p0e23-generation-Maxwell-reco-$typeReco-$V.dat"))'

    nnn1 = Nn[regionETL1] .* Fcc[iphin].(zn*( q*(view(sol1[iphin, :], subgn).-view(sol1[ipsi, :], subgn)) .+ En[regionETL1])./(kB*T))
    nnn2 = Nn[regionETL1] .* Fcc[iphin].(zn*( q*(view(sol2[iphin, :], subgn).-view(sol2[ipsi, :], subgn)) .+ En[regionETL1])./(kB*T))
    nnn4 = Nn[regionETL1] .* Fcc[iphin].(zn*( q*(view(sol4[iphin, :], subgn).-view(sol4[ipsi, :], subgn)) .+ En[regionETL1])./(kB*T))

    nn1  = Nn[regionPero] .* Fcc[iphin].(zn*( q*(view(sol1[iphin, :], subg).-view(sol1[ipsi, :], subg)) .+ En[regionPero])./(kB*T))
    nn2  = Nn[regionPero] .* Fcc[iphin].(zn*( q*(view(sol2[iphin, :], subg).-view(sol2[ipsi, :], subg)) .+ En[regionPero])./(kB*T))
    nn4  = Nn[regionPero] .* Fcc[iphin].(zn*( q*(view(sol4[iphin, :], subg).-view(sol4[ipsi, :], subg)) .+ En[regionPero])./(kB*T))

    nnp1 = Nn[regionHTL] .* Fcc[iphin].(zn*( q*(view(sol1[iphin, :], subgp).-view(sol1[ipsi, :], subgp)) .+ En[regionHTL])./(kB*T))
    nnp2 = Nn[regionHTL] .* Fcc[iphin].(zn*( q*(view(sol2[iphin, :], subgp).-view(sol2[ipsi, :], subgp)) .+ En[regionHTL])./(kB*T))
    nnp4 = Nn[regionHTL] .* Fcc[iphin].(zn*( q*(view(sol4[iphin, :], subgp).-view(sol4[ipsi, :], subgp)) .+ En[regionHTL])./(kB*T))

    npn1 = Np[regionETL1] .* Fcc[iphip].(zp*( q*(view(sol1[iphip, :], subgn).-view(sol1[ipsi, :], subgn)) .+ Ep[regionETL1])./(kB*T))
    npn2 = Np[regionETL1] .* Fcc[iphip].(zp*( q*(view(sol2[iphip, :], subgn).-view(sol2[ipsi, :], subgn)) .+ Ep[regionETL1])./(kB*T))
    npn4 = Np[regionETL1] .* Fcc[iphip].(zp*( q*(view(sol4[iphip, :], subgn).-view(sol4[ipsi, :], subgn)) .+ Ep[regionETL1])./(kB*T))

    np1  = Np[regionPero] .* Fcc[iphip].(zp*( q*(view(sol1[iphip, :], subg).-view(sol1[ipsi, :], subg)) .+ Ep[regionPero])./(kB*T))
    np2  = Np[regionPero] .* Fcc[iphip].(zp*( q*(view(sol2[iphip, :], subg).-view(sol2[ipsi, :], subg)) .+ Ep[regionPero])./(kB*T))
    np4  = Np[regionPero] .* Fcc[iphip].(zp*( q*(view(sol4[iphip, :], subg).-view(sol4[ipsi, :], subg)) .+ Ep[regionPero])./(kB*T))

    npp1 = Np[regionHTL] .* Fcc[iphip].(zp*( q*(view(sol1[iphip, :], subgp).-view(sol1[ipsi, :], subgp)) .+ Ep[regionHTL])./(kB*T))
    npp2 = Np[regionHTL] .* Fcc[iphip].(zp*( q*(view(sol2[iphip, :], subgp).-view(sol2[ipsi, :], subgp)) .+ Ep[regionHTL])./(kB*T))
    npp4 = Np[regionHTL] .* Fcc[iphip].(zp*( q*(view(sol4[iphip, :], subgp).-view(sol4[ipsi, :], subgp)) .+ Ep[regionHTL])./(kB*T))

    na1  = Na[regionPero] .* Fcc[iphia].(za*( q*(view(sol1[iphia, :], subg).-view(sol1[ipsi, :], subg)) .+ Ea1)./(kB*T))
    na2  = Na[regionPero] .* Fcc[iphia].(za*( q*(view(sol2[iphia, :], subg).-view(sol2[ipsi, :], subg)) .+ Ea2)./(kB*T))
    na4  = Na[regionPero] .* Fcc[iphia].(za*( q*(view(sol4[iphia, :], subg).-view(sol4[ipsi, :], subg)) .+ Ea4)./(kB*T))

    #########################################################################################################
    #########################################################################################################
    #########################################################################################################

    Blues   = get_cmap(:Blues)
    Oranges = get_cmap(:Oranges)
    Gold    = get_cmap(:Wistia)

    #########################################

    semilogy(coordn', nnn1, linewidth = 5, color = Blues(201))
    semilogy(coord',  nn1,  linewidth = 5, color = Blues(201))
    semilogy(coordp', nnp1, linewidth = 5, color = Blues(201))

    semilogy(coordn', npn1, linewidth = 5, color = Oranges(201))
    semilogy(coord',  np1,  linewidth = 5, color = Oranges(201))
    semilogy(coordp', npp1, linewidth = 5, color = Oranges(201))

    semilogy(coord',  na1,  linewidth = 5, color = Gold(201))

    PyPlot.xlabel("\$ y\$  [nm]", fontsize=18)
    PyPlot.ylabel("Density [\$\\frac{1}{\\mathrm{m}^3}\$]", fontsize=18)
    PyPlot.xlim(-5, 455)
    PyPlot.ylim(1.0e14, 5.0e24)
    PyPlot.xticks([0.0, 200, 400])
    PyPlot.yticks([1.0e14, 1.0e17, 1.0e20, 1.0e23])
    PyPlot.tick_params(which ="both", labelsize=18)
    PyPlot.tight_layout()

    if saveFig
        savefig(datadir("1D-dens-Ca-1p0e21-scanrate-$textSR-generation-Maxwell-$IVDirection-$V.pdf"))
    end

    #########################################
    figure()

    semilogy(coordn', nnn2, linewidth = 5, color = Blues(201))
    semilogy(coord',  nn2,  linewidth = 5, color = Blues(201))
    semilogy(coordp', nnp2, linewidth = 5, color = Blues(201))

    semilogy(coordn', npn2, linewidth = 5, color = Oranges(201))
    semilogy(coord',  np2,  linewidth = 5, color = Oranges(201))
    semilogy(coordp', npp2, linewidth = 5, color = Oranges(201))

    semilogy(coord',  na2,  linewidth = 5, color = Gold(201))

    PyPlot.xlabel("\$ y\$ [nm]", fontsize=18)
    PyPlot.ylabel("Density [\$\\frac{1}{\\mathrm{m}^3}\$]", fontsize=18)
    PyPlot.xlim(-5, 455)
    PyPlot.ylim(1.0e14, 5.0e24)
    PyPlot.xticks([0.0, 200, 400])
    PyPlot.yticks([1.0e14, 1.0e17, 1.0e20, 1.0e23])
    PyPlot.tick_params(which ="both", labelsize=18)
    PyPlot.tight_layout()

    if saveFig
        savefig(datadir("1D-dens-Ca-1p0e22-scanrate-$textSR-generation-Maxwell-$IVDirection-$V.pdf"))
    end

    #########################################
    figure()

    semilogy(coordn', nnn4, linewidth = 5, color = Blues(201))
    semilogy(coord', nn4, linewidth = 5, color = Blues(201))
    semilogy(coordp', nnp4, linewidth = 5, color = Blues(201))

    semilogy(coordn', npn4, linewidth = 5, color = Oranges(201))
    semilogy(coord', np4, linewidth = 5, color = Oranges(201))
    semilogy(coordp', npp4, linewidth = 5, color = Oranges(201))

    semilogy(coord', na4, linewidth = 5, color = Gold(201))

    PyPlot.xlabel("\$ y\$  [nm]", fontsize=18)
    PyPlot.ylabel("Density [\$\\frac{1}{\\mathrm{m}^3}\$]", fontsize=18)
    PyPlot.xlim(-5, 455)
    PyPlot.ylim(1.0e14, 5.0e24)
    PyPlot.xticks([0.0, 200, 400])
    PyPlot.yticks([1.0e14, 1.0e17, 1.0e20, 1.0e23])
    PyPlot.tick_params(which ="both", labelsize=18)
    PyPlot.tight_layout()

    if saveFig
        savefig(datadir( "1D-dens-Ca-1p0e23-scanrate-$textSR-generation-Maxwell-$IVDirection-$V.pdf"))
    end

    #########################################################################################
    #########################################################################################

    ## calculate the measure of the perovskite region
    mOmega = 0.0
    for icellVol in subg[CellVolumes]
        mOmega = mOmega + icellVol
    end

    #######################

    ctsys.fvmsys.physics.data.params.bandEdgeEnergy[iphia, regionPero] = Ea1

    Int1     =  ChargeTransport.integrate(ctsys, storage!, sol1)./q
    nnAvg1   = -Int1[iphin, regionPero]/mOmega
    npAvg1   =  Int1[iphip, regionPero]/mOmega
    naAvg1   =  Int1[iphia, regionPero]/mOmega

    ######################################################################

    ctsys.fvmsys.physics.data.params.bandEdgeEnergy[iphia, regionPero] = Ea2

    Int2     =  ChargeTransport.integrate(ctsys, storage!, sol2)./q
    nnAvg2   = -Int2[iphin, regionPero]/mOmega
    npAvg2   =  Int2[iphip, regionPero]/mOmega
    naAvg2   =  Int2[iphia, regionPero]/mOmega

    ######################################################################

    ctsys.fvmsys.physics.data.params.bandEdgeEnergy[iphia, regionPero] = Ea4

    Int4     =  ChargeTransport.integrate(ctsys, storage!, sol4)./q
    nnAvg4   = -Int4[iphin, regionPero]/mOmega
    npAvg4   =  Int4[iphip, regionPero]/mOmega
    naAvg4   =  Int4[iphia, regionPero]/mOmega

   ######################################################################

    if printText
       println("  ")
       println("Avg nn for Ca = $Ca1 is:                 $(nnAvg1)")
       println("Avg nn for Ca = $Ca2 is:                 $(nnAvg2)")
       println("Avg nn for Ca = $Ca4 is:                 $(nnAvg4)")

       println("  ")
       println("Avg np for Ca = $Ca1 is:                 $(npAvg1)")
       println("Avg np for Ca = $Ca2 is:                 $(npAvg2)")
       println("Avg np for Ca = $Ca4 is:                 $(npAvg4)")

       println("  ")
       println("Avg na for Ca = $Ca1 is:                 $(naAvg1)")
       println("Avg na for Ca = $Ca2 is:                 $(naAvg2)")
       println("Avg na for Ca = $Ca4 is:                 $(naAvg4)")
    end

    #########################################################################################
    #########################################################################################

    Enn1 = En[regionETL1]./q .- view(sol1[ipsi, :], subgn)
    Enn2 = En[regionETL1]./q .- view(sol2[ipsi, :], subgn)
    Enn4 = En[regionETL1]./q .- view(sol4[ipsi, :], subgn)

    En1 = En[regionPero]./q .- view(sol1[ipsi, :], subg)
    En2 = En[regionPero]./q .- view(sol2[ipsi, :], subg)
    En4 = En[regionPero]./q .- view(sol4[ipsi, :], subg)

    Enp1 = En[regionHTL]./q .- view(sol1[ipsi, :], subgp)
    Enp2 = En[regionHTL]./q .- view(sol2[ipsi, :], subgp)
    Enp4 = En[regionHTL]./q .- view(sol4[ipsi, :], subgp)

    Epn1 = Ep[regionETL1]./q .- view(sol1[ipsi, :], subgn)
    Epn2 = Ep[regionETL1]./q .- view(sol2[ipsi, :], subgn)
    Epn4 = Ep[regionETL1]./q .- view(sol4[ipsi, :], subgn)

    Ep1 = Ep[regionPero]./q .- view(sol1[ipsi, :], subg)
    Ep2 = Ep[regionPero]./q .- view(sol2[ipsi, :], subg)
    Ep4 = Ep[regionPero]./q .- view(sol4[ipsi, :], subg)

    Epp1 = Ep[regionHTL]./q .- view(sol1[ipsi, :], subgp)
    Epp2 = Ep[regionHTL]./q .- view(sol2[ipsi, :], subgp)
    Epp4 = Ep[regionHTL]./q .- view(sol4[ipsi, :], subgp)

    #########################################################################################
    #########################################################################################

    LS = "--"

    figure()
    plot(coordn',    Enn1, linewidth = 5, color = Blues(241))
    plot(coord',     En1,  linewidth = 5, color = Blues(241))
    plot(coordp',    Enp1, linewidth = 5, color = Blues(241))

    plot(coordn',    Epn1, linewidth = 5, color = Oranges(241))
    plot(coord',     Ep1,  linewidth = 5, color = Oranges(241))
    plot(coordp',    Epp1, linewidth = 5, color = Oranges(241))

    plot(coordALL', -sol1[iphin, :], linewidth = 5, linestyle = LS, color = Blues(201))
    plot(coordALL', -sol1[iphip, :], linewidth = 5, linestyle = LS, color = Oranges(201))

    PyPlot.xlabel("\$ y\$  [nm]",   fontsize=18)
    PyPlot.ylabel(" Energy [eV] ", fontsize=18)
    PyPlot.xlim(-5.0, 445)
    PyPlot.ylim(-2.4, 1.8)
    PyPlot.xticks([0.0, 200, 400])
    PyPlot.tick_params(which ="both", labelsize=18)
    PyPlot.tight_layout()

    if saveFig
        savefig(datadir("1D-energy-Ca-1p0e21-scanrate-$textSR-generation-Maxwell-$IVDirection-$V.pdf"))
    end

    #########################################
    figure()
    plot(coordn',    Enn2, linewidth = 5, color = Blues(241))
    plot(coord',     En2,  linewidth = 5, color = Blues(241))
    plot(coordp',    Enp2, linewidth = 5, color = Blues(241))

    plot(coordn',    Epn2, linewidth = 5, color = Oranges(241))
    plot(coord',     Ep2,  linewidth = 5, color = Oranges(241))
    plot(coordp',    Epp2, linewidth = 5, color = Oranges(241))

    plot(coordALL', -sol2[iphin, :], linewidth = 5, linestyle = LS, color = Blues(201))
    plot(coordALL', -sol2[iphip, :], linewidth = 5, linestyle = LS, color = Oranges(201))

    PyPlot.xlabel("\$ y\$  [nm]",   fontsize=18)
    PyPlot.ylabel(" Energy [eV] ", fontsize=18)
    PyPlot.xlim(-5.0, 445)
    PyPlot.ylim(-2.4, 1.8)
    PyPlot.xticks([0.0, 200, 400])
    PyPlot.tick_params(which ="both", labelsize=18)
    PyPlot.tight_layout()

    if saveFig
        savefig(datadir("1D-energy-Ca-1p0e22-scanrate-$textSR-generation-Maxwell-$IVDirection-$V.pdf"))
    end

    #########################################
    figure()
    plot(coordn',   Enn4, linewidth = 5, color = Blues(241))
    plot(coord',    En4,  linewidth = 5, color = Blues(241))
    plot(coordp',   Enp4, linewidth = 5, color = Blues(241))

    plot(coordn',   Epn4, linewidth = 5, color = Oranges(241))
    plot(coord',    Ep4,  linewidth = 5, color = Oranges(241))
    plot(coordp',   Epp4, linewidth = 5, color = Oranges(241))

    plot(coordALL', -sol4[iphin, :], linewidth = 5, linestyle = LS, color = Blues(201))
    plot(coordALL', -sol4[iphip, :], linewidth = 5, linestyle = LS, color = Oranges(201))

    PyPlot.xlabel("\$ y\$  [nm]",   fontsize=18)
    PyPlot.ylabel(" Energy [eV] ", fontsize=18)
    PyPlot.xlim(-5.0, 445)
    PyPlot.ylim(-2.4, 1.8)
    PyPlot.xticks([0.0, 200, 400])
    PyPlot.tick_params(which ="both", labelsize=18)
    PyPlot.tight_layout()

    if saveFig
        savefig(datadir("1D-energy-Ca-1p0e23-scanrate-$textSR-generation-Maxwell-$IVDirection-$V.pdf"))
    end

    return nothing
end

end