

#=

Simulation of a three-layer perovskite solar cell in 2D with textured interface.

=#

ENV["LC_NUMERIC"]="C"

module SingleJunction

using TexturedPerovskiteSolarCells
using ChargeTransport
using ExtendableGrids
using PyPlot
using LinearAlgebra
using DelimitedFiles
using LinearSolve
using SimplexGridFactory
using Triangulate
using VoronoiFVM
using Roots

# for convenience
datadir = TexturedPerovskiteSolarCells.datadir
scriptsdir = TexturedPerovskiteSolarCells.scriptsdir

###########################################################
###########################################################

function main(;plotting = false, printText = true,
               ########################
               parameter_file = scriptsdir("params_single_junction.jl"),
               ########################
               gridDim = 1,
               ########################
               typeGrid  = "planar",  # "nanotextured", #
               amplitude = 4.0e-7,    # amplitude of nanotexture
               ########################
               typeReco =  "all" ,    # "radiative", #"off", # "radiative", # "SR", # "bulk", # "all"   (SR for surface reco)
               ########################
               generation = true, generationUniform = false, MaxwellSol = true,
               ########################
               scanrate = 1.0e3  * V/s, # 1.0e-3, # 1.0e1, #
               ########################
               CalculateEa = false, EaLoop = -4.0, # for calculating Ea of φa
               ########################
               enableIons = true, saveData = false,
               ########################
               demo_run = false,    # do calculations on coarser mesh
               ########################
               plotPostProcess = false,
               ########################
               )

    PyPlot.rc("font", family="sans-serif", size=14)
    PyPlot.rc("mathtext", fontset="dejavusans")
    PyPlot.close("all")

    # we need this, since methods coincide with VoronoiFVM methods
    SolverControl()                      = ChargeTransport.SolverControl()
    unknowns(a)                          = ChargeTransport.unknowns(a)
    TestFunctionFactory(a)               = ChargeTransport.TestFunctionFactory(a)
    testfunction(a,b,c)                  = ChargeTransport.testfunction(a,b,c)
    integrate(a, b, c, d, e; kwarges...) = ChargeTransport.integrate(a, b, c, d, e; kwarges...)
    integrate(a, b, c; kwarges...)       = ChargeTransport.integrate(a, b, c; kwarges...)

    if gridDim == 1
        typeGrid = "planar"
    end
    ################################################################################
    if printText
        println("Set up grid and regions")
    end
    ################################################################################

    include(parameter_file)

    grid = generate_grid(gridDim = gridDim, type = typeGrid, amplitude = amplitude, parameter_file = parameter_file, demo_run = demo_run)

    if plotting
        gridplot(grid, Plotter= PyPlot, resolution=(600,400), linewidth=0.5, legend=:rc)
    end

    ## planar:                       47 122 nodes
    ## nanotexture (ampl = 2.0e-7):  52 597 nodes
    ## nanotexture (ampl = 4.0e-7):  71 657 nodes
    ## nanotexture (ampl = 6.5e-7): 118 557 nodes
    ## nanotexture (ampl = 8.0e-7): 170 590 nodes

    if printText
        println("*** done\n")
    end

    ################################################################################
    if printText
        println("Define physical parameters and model")
    end
    ################################################################################

    ########## primary data for I-V scan protocol ##########

    if typeReco == "all" # radiative until 1.4 V and reco = all until 1.2 V
        endVoltage       = 1.2   * V
    else
        endVoltage       = 1.4   * V
    end
    tPrecond             = 5 * s  * 1.0e-1/scanrate

    tend                 = endVoltage/scanrate

    ## Define scan protocol function
    function scanProtocol(t)
        if data.calculationType == InEquilibrium
            biasVal = 0.0
        else
            if        0.0      <= t && t <= tPrecond
                    biasVal = endVoltage
            elseif    tPrecond <= t  && t <= tPrecond + tend
                biasVal = endVoltage .- scanrate * (t - tPrecond)
            elseif     tend + tPrecond < t  && t <= tPrecond + 2*tend
                biasVal = 0.0 .+ scanrate * (t - tend- tPrecond)
            else
                biasVal = 0.0
            end
        end
        return biasVal
    end

    # Apply zero voltage on left boundary and a linear scan protocol on right boundary
    contactVoltageFunction = [zeroVoltage, scanProtocol]

    ########################################

    if gridDim == 1
        if generationUniform
            Ea = EaPlanar
        else
            Ea = Ea1D
        end
    end

    ## DA: as we adjust the geometry, the average vacancy density also changes, when not adjusting the intrinsic energy level for vacancies.
    if gridDim == 2
        if typeGrid == "planar"
            Ea = EaPlanar
        elseif typeGrid == "nanotextured"
            if     amplitude == 0.5e-7
                Ea = EaAmpl0p5e7
            elseif amplitude == 1.0e-7
                Ea = EaAmpl1p0e7
            elseif amplitude == 1.5e-7
                Ea = EaAmpl1p5e7
            elseif amplitude == 2.0e-7
                Ea = EaAmpl2p0e7
            elseif amplitude == 2.5e-7
                Ea = EaAmpl2p5e7
            elseif amplitude == 3.0e-7
                Ea = EaAmpl3p0e7
            elseif amplitude == 3.5e-7
                Ea = EaAmpl3p5e7
            elseif amplitude == 4.0e-7
                Ea = EaAmpl4p0e7
            elseif amplitude == 4.5e-7
                Ea = EaAmpl4p5e7
            elseif amplitude == 5.0e-7
                Ea = EaAmpl5p0e7
            elseif amplitude == 5.5e-7
                Ea = EaAmpl5p5e7
            elseif amplitude == 6.0e-7
                Ea = EaAmpl6p0e7
            elseif amplitude == 6.5e-7
                Ea = EaAmpl6p5e7
            elseif amplitude == 7.0e-7
                Ea = EaAmpl7p0e7
            elseif amplitude == 7.5e-7
                Ea = EaAmpl7p5e7
            elseif amplitude == 8.0e-7
                Ea = EaAmpl8p0e7
            end
        end
    end

    ## DA: Here, we adjusted a bit the average vacancy density to see what is the origin of the imbalance between electrons and holes.
    if Ca == 1.0e21
        #int Ea cal: 9.441715265086726e20  #(int - Ca)/Ca = -0.05582847349132742
        Ea = [0.0,  -5.614,   0.0] .*  eV
    elseif Ca == 1.0e22
        # int Ea cal: 9.910724685807422e21 # (int - Ca) / Ca = -0.0089275314192578
        Ea = [0.0,  -5.532,   0.0] .*  eV
    elseif Ca == 1.0e23
        # int Ea cal: 9.92995947626143e22  # (int - Ca) / Ca = -0.0070040523738568675
        Ea = [0.0,  -5.283,   0.0] .*  eV
    elseif Ca == 1.0e24
        # int Ea cal: 9.601509074576119e23 # (int - Ca) / Ca = -0.03984909254238813
        Ea = [0.0,  -5.3,   0.0] .*  eV
    end

    if CalculateEa
        Ea[regionPero] =  EaLoop * eV
    end

    # For all other parameters, we refer to the parameters template
    if printText
        println("*** done\n")
    end

    helpSR = collect(string(scanrate));  helpSR[ findall(x -> x == '.', helpSR)[1] ] = 'p'
    textSR = join(helpSR)

    ################################################################################
    if printText
        println("Define System and fill in information about model")
    end
    ################################################################################

    if enableIons
        numberOfCarriers2 = numberOfCarriers
    else
        numberOfCarriers2 = 2
    end

    # initialize Data instance and fill in data
    if generation
        if MaxwellSol && generationUniform == false

            generationData = MaxwellPhotogeneration(gridDim = gridDim, typeGrid = typeGrid, amplitude = amplitude, parameter_file = parameter_file, demo_run = demo_run)

        elseif MaxwellSol== false && generationUniform == false
            Fph     = incidentPhotonFlux[regionPero]; ag  = absorption[regionPero]
            genPeak = generationPeak;                 inv = invertedIllumination

            function BeerLamb(x)
                if heightLayersPL[1] <= x <= heightLayersPL[2]
                    Fph .* ag .* exp.( - inv .* ag .* (x .- genPeak))
                else
                    0.0
                end
            end

            if gridDim == 1
                coord  = vec(grid[Coordinates])
            elseif gridDim == 2
                coord  = vec(grid[Coordinates][2, :])
            end

            generationData           = BeerLamb.(coord)

        elseif generationUniform

            generationData           = zeros(length(grid[Coordinates][1,:]))
            subg                     = subgrid(grid, [regionPero])
            iNode                    = subg[NodeParents]

            generationData[iNode]   .= 3.35e27 * 1/(m^3 * s) # average Photogeneration for planar

        end

        data                         = Data(grid, numberOfCarriers2, contactVoltageFunction = contactVoltageFunction, generationData = generationData)
    else
        data                         = Data(grid, numberOfCarriers2, contactVoltageFunction = contactVoltageFunction)
    end

    data.modelType                   = Transient

    if typeReco == "off"
        bulk_recomb_radiative = false; bulk_recomb_SRH = false
    elseif typeReco == "radiative"
        bulk_recomb_radiative = true; bulk_recomb_SRH = false
    elseif typeReco == "SRH"
        bulk_recomb_radiative = false; bulk_recomb_SRH = true
    elseif typeReco == "SR"
        bulk_recomb_radiative = false; bulk_recomb_SRH = false
        data.boundaryType[bregionJ1] = InterfaceRecombination
        data.boundaryType[bregionJ2] = InterfaceRecombination
    elseif typeReco == "bulk"
        bulk_recomb_radiative = true;  bulk_recomb_SRH = true
    elseif typeReco == "all"
        bulk_recomb_radiative = true;  bulk_recomb_SRH = true
        data.boundaryType[bregionJ1] = InterfaceRecombination
        data.boundaryType[bregionJ2] = InterfaceRecombination
    end

    data.bulkRecombination           = set_bulk_recombination(;iphin = iphin, iphip = iphip, bulk_recomb_Auger = false,
                                                              bulk_recomb_radiative = bulk_recomb_radiative,
                                                              bulk_recomb_SRH = bulk_recomb_SRH)

    data.boundaryType[bregionLeft]   = MixedOhmicSchottkyContact
    data.boundaryType[bregionRight]  = MixedOhmicSchottkyContact

    if enableIons
        data.F                       = Fcc
        enable_ionic_carrier!(data, ionicCarrier = iphia, regions = [regionPero])
    else
        data.F                       = [Boltzmann, Boltzmann]
    end

    if generation
        if MaxwellSol
            generationType = "Maxwell"
        else
            generationType = "Beer-Lambert"
        end
    else
        generationType     = "none"
    end

    if generation && generationUniform
        generationType     = "uniform"
    end

    if generation
        data.generationModel = GenerationUserDefined
    end

    if printText
        println("*** done\n")
    end

    ################################################################################
    if printText
        println("Define Params and fill in physical parameters")
    end
    ################################################################################

    params                                                    = Params(grid, numberOfCarriers2)

    params.temperature                                        = T
    params.UT                                                 = (kB * params.temperature) / q
    params.chargeNumbers[iphin]                               = zn
    params.chargeNumbers[iphip]                               = zp
    if enableIons
        params.chargeNumbers[iphia]                           = za
    end

    for ireg in 1:numberOfRegions # interior region data

        params.dielectricConstant[ireg]                       = εr[ireg] * ε0

        ## effective DOS, band edge energy and mobilities
        params.densityOfStates[iphin, ireg]                   = Nn[ireg]
        params.densityOfStates[iphip, ireg]                   = Np[ireg]

        params.bandEdgeEnergy[iphin, ireg]                    = En[ireg]
        params.bandEdgeEnergy[iphip, ireg]                    = Ep[ireg]

        params.mobility[iphin, ireg]                          = μn[ireg]
        params.mobility[iphip, ireg]                          = μp[ireg]

        ## vacancy parameters
        if enableIons
            params.densityOfStates[iphia, ireg]               = Na[ireg]
            params.bandEdgeEnergy[iphia, ireg]                = Ea[ireg]
            params.mobility[iphia, ireg]                      = μa[ireg]

        end

        ## recombination parameters
        params.recombinationRadiative[ireg]                   = r0[ireg]
        params.recombinationSRHLifetime[iphin, ireg]          = τn[ireg]
        params.recombinationSRHLifetime[iphip, ireg]          = τp[ireg]
        params.recombinationSRHTrapDensity[iphin, ireg]       = nτ[ireg]
        params.recombinationSRHTrapDensity[iphip, ireg]       = pτ[ireg]

    end

    params.SchottkyBarrier[bregionLeft]                       =  0.1 * eV
    params.SchottkyBarrier[bregionRight]                      = -0.1 * eV + En[3] - Ep[3]

    ##############################################################
    ## inner boundary region data (we choose the intrinsic values)
    if typeReco == "SR" || typeReco == "all"
        params.bDensityOfStates[iphin, bregionJ1]             = Nn[regionETL1]
        params.bDensityOfStates[iphip, bregionJ1]             = Np[regionPero]

        params.bDensityOfStates[iphin, bregionJ2]             = Nn[regionPero]
        params.bDensityOfStates[iphip, bregionJ2]             = Np[regionHTL]

        params.bBandEdgeEnergy[iphin, bregionJ1]              = En[regionETL1]
        params.bBandEdgeEnergy[iphip, bregionJ1]              = Ep[regionPero]

        params.bBandEdgeEnergy[iphin, bregionJ2]              = En[regionPero]
        params.bBandEdgeEnergy[iphip, bregionJ2]              = Ep[regionHTL]

        ## for surface recombination
        params.recombinationSRHvelocity[iphin, bregionJ1]     = SRHvelocityETLn
        params.recombinationSRHvelocity[iphip, bregionJ1]     = SRHvelocityETLp

        params.bRecombinationSRHTrapDensity[iphin, bregionJ1] = params.recombinationSRHTrapDensity[iphin, regionETL1]
        params.bRecombinationSRHTrapDensity[iphip, bregionJ1] = params.recombinationSRHTrapDensity[iphip, regionPero]

        params.recombinationSRHvelocity[iphin, bregionJ2]     = SRHvelocityHTLn
        params.recombinationSRHvelocity[iphip, bregionJ2]     = SRHvelocityHTLp

        params.bRecombinationSRHTrapDensity[iphin, bregionJ2] = params.recombinationSRHTrapDensity[iphin, regionPero]
        params.bRecombinationSRHTrapDensity[iphip, bregionJ2] = params.recombinationSRHTrapDensity[iphip, regionHTL]
    end

    ##############################################################

    ## interior doping
    params.doping[iphip, regionHTL]                           = Cp
    params.doping[iphin, regionETL1]                          = Cn1

    # doping
    if enableIons
        params.doping[iphia, regionPero]                      = Ca
    end

    data.params = params
    ctsys       = System(grid, data, unknown_storage=:sparse)

    ipsi        = data.index_psi

    if plotPostProcess
        ctsys.fvmsys.physics.data.calculationType = OutOfEquilibrium
        return grid, ctsys
    end


    if printText
        println("*** done\n")
    end

    ################################################################################
    if printText
        println("Define control parameters for Newton solver")
    end
    ################################################################################

    control              = SolverControl()
    if printText
        control.verbose  = "e"
    else
        control.verbose  = "d"
    end
    control.method_linear = UMFPACKFactorization()
    control.damp_initial = damp_initial
    control.damp_growth  = damp_growth
    control.max_round    = max_round
    control.maxiters     = maxiters

    control.abstol       = abstol
    control.reltol       = reltol
    control.tol_round    = tol_round

    control.Δu_opt       = Inf

    if demo_run
        # precond: ~ 60 time steps (1 sec); reverse: ~ 150 time steps (~8 sec); forward: ~ 150 time steps (~8 sec)
        control.Δt       = 3.0e-4/scanrate
        control.Δt_min   = 6.0e-3/scanrate
        control.Δt_max   = 8.0e-3/scanrate
        control.Δt_grow  = 1.01
    else
        # precond: ~ 150 time steps; reverse: ~ 460 time steps (~80 sec); forward: ~ 460 time steps (~80 sec)
        control.Δt       = 7.0e-4/scanrate
        control.Δt_min   = 1.0e-4/scanrate
        control.Δt_max   = 5.0e-3/scanrate
        control.Δt_grow  = 1.005
    end

    if printText
        println("*** done\n")
    end

    ################################################################################
    if printText
        println("Compute solution in thermodynamic equilibrium")
    end
    ################################################################################

    # initialize solution and starting vectors
    inival  = unknowns(ctsys); initialCond  = unknowns(ctsys)
    sol     = unknowns(ctsys); sol0p0       = unknowns(ctsys)
    inival .= 0.0;             initialCond .= 0.0

    solEQ   = equilibrium_solve!(ctsys, control = control, nonlinear_steps = 0)

    inival  = solEQ

    if demo_run && gridDim == 2 && generation && MaxwellSol
        if typeGrid == "planar"
            return testval = sum(filter(!isnan, inival))/length(inival) # when using sparse storage, we get NaN values in solution
        elseif typeGrid == "nanotextured" && amplitude == 2.0e-7
            return testval = sum(filter(!isnan, inival))/length(inival) # when using sparse storage, we get NaN values in solution
        end
    end

    if plotting

        label_solution, label_density, label_energy, label_BEE = set_plotting_labels(data)

        if enableIons
            ## add labels for anion vacancy
            label_energy[1, iphia] = "\$E_a-q\\psi\$"; label_energy[2, iphia] = "\$ - q \\varphi_a\$"; label_BEE[iphia] = "\$E_a\$"
            label_density[iphia]   = "\$ n_a \$";      label_solution[iphia]  = "\$ \\varphi_a\$"
        end

    end

    if printText
        println("*** done\n")
    end

    ################################################################################
     if printText  && generation
        println("Loop for generation")
    end
    ################################################################################

    if generation

        control.damp_initial = 0.01
        control.damp_growth  = 1.21 # >= 1
        control.max_round    = 0
        control.abstol       = 1.0e-3
        control.reltol       = 1.0e-3
        control.tol_round    = 1.0e-3

        if amplitude == 8.0e-7 && typeReco == "all"
            control.max_round    = 0
            control.damp_initial = 0.01
            control.abstol       = 1.0e-4
            control.reltol       = 1.0e-4
            control.tol_round    = 1.0e-4
        end

        # these values are needed for putting the generation slightly on
        I      = collect(20:-1:0.0)
        LAMBDA = 10 .^ (-I)

        set_contact!(ctsys, bregionRight, Δu = -endVoltage)

        for istep = 1:length(I)-1

            ## turn slowly generation on
            ctsys.data.λ2   = LAMBDA[istep + 1]

            if printText
                println("increase generation with λ2 = $(data.λ2)")
            end

            solEQ  = ChargeTransport.solve(ctsys, inival = inival, control = control)
            inival = solEQ

        end # generation loop

        inival  = solEQ

        if printText
            println("*** done\n")
        end

    end

    if CalculateEa
        return ctsys, inival
    end

    ################################################################################
    if printText
        println("Loop for increasing bias")
    end
    ################################################################################

    control.damp_initial = 0.5
    control.max_round    = 5
    control.abstol       = abstol
    control.reltol       = reltol
    control.tol_round    = tol_round

    biasLoop = collect(reverse(-range(0.0, endVoltage, length = 31)))

    if amplitude == 8.0e-7 && typeReco == "all"
        I2      = collect(31:-1:0.0)
        LAMBDA2 = 10 .^ (-I2)
    end

    for istep = 1:length(biasLoop)

        ## turn slowly voltage on
        set_contact!(ctsys, bregionRight, Δu = biasLoop[istep])

        if printText
            println(" bias =" , biasLoop[istep])
        end

        if amplitude == 8.0e-7 && typeReco == "all"

            ctsys.data.params.recombinationSRHvelocity[iphin, bregionJ1] = LAMBDA2[istep] * SRHvelocityETLn
            ctsys.data.params.recombinationSRHvelocity[iphip, bregionJ2] = LAMBDA2[istep] * SRHvelocityHTLp

            @show ctsys.data.params.recombinationSRHvelocity[iphin, bregionJ1]
            @show ctsys.data.params.recombinationSRHvelocity[iphip, bregionJ2]
            println(" ")
        end

        sol0p0 = ChargeTransport.solve(ctsys, inival = inival, control = control, tstep = 1.0)
        inival = sol0p0

    end # bias loop

    if plotting

        if gridDim == 1
            figure()
            plot_densities(PyPlot, ctsys, sol0p0, "Initial condition", label_density)
            figure()
            plot_solution(PyPlot, ctsys,  sol0p0, "Initial condition", label_solution)
        end

    end

    if printText
        println("*** done\n")
    end

    ################################################################################
    if printText
        println("Preconditioning")
    end
    ################################################################################

    if demo_run
        control.Δt_grow = 1.03
    else
        control.Δt_grow = 1.02
    end
    solPrecond      = ChargeTransport.solve(ctsys, inival = inival, times=(0.0, tPrecond), control = control)

    if printText
        println("*** done\n")
    end

    ################################################################################
    if printText
        println("Reverse IV Measurement loop")
    end
    ################################################################################

    if demo_run
        control.Δt_grow = 1.01
    else
        control.Δt_grow = 1.005
    end
    solRev          = ChargeTransport.solve(ctsys, inival = solPrecond.u[end], times=(tPrecond, tPrecond + tend),   control = control)

    if printText
        println("*** done\n")
    end

    ################################################################################
    if printText
        println("Reverse IV curve calculation")
    end
    ################################################################################

    IVRev         = zeros(0) # for saving I-V data
    IVRevn        = zeros(0); IVRevp  = zeros(0)
    IVReva        = zeros(0); IVRevψ  = zeros(0)

    ######################
    IRevSRHn      = zeros(0); IRevRadn = zeros(0)
    IRevSRHp      = zeros(0); IRevRadp = zeros(0)
    IRevSRnL      = zeros(0); IRevSRnR = zeros(0)
    IRevSRpL      = zeros(0); IRevSRpR = zeros(0)

    tvaluesRev    = solRev.t
    number_tsteps = length(tvaluesRev)
    biasValuesRev = scanProtocol.(tvaluesRev[2:end])

    factory       = TestFunctionFactory(ctsys)
    tf            = testfunction(factory, [bregionLeft], [bregionRight])

    for istep = 2:number_tsteps

        Δt       = tvaluesRev[istep] - tvaluesRev[istep-1] # Time step size
        inival   = solRev.u[istep-1]
        solution = solRev.u[istep]

        I        = integrate(ctsys, tf, solution, inival, Δt)

        current = 0.0
        for ii = 1:numberOfCarriers2+1
            current = current + I[ii]
        end

        push!(IVRev,   current); push!(IVRevn, I[iphin])
        push!(IVRevp, I[iphip]); push!(IVRevψ, I[ipsi])

        if enableIons
            push!(IVReva, I[iphia])
        end

        IntSRH         = integrate(ctsys, SRHRecombination!, solution)
        IntRad         = integrate(ctsys, RadiativeRecombination!, solution)
        IntSR          = integrate(ctsys, SRRecombination!, solution, boundary = true)

        IntSRHnSum     = 0.0; IntRadnSum = 0.0
        IntSRHpSum     = 0.0; IntRadpSum = 0.0

        for ii = 1:numberOfRegions
            IntSRHnSum = IntSRHnSum - IntSRH[iphin, ii]
            IntRadnSum = IntRadnSum - IntRad[iphin, ii]

            IntSRHpSum = IntSRHpSum + IntSRH[iphip, ii]
            IntRadpSum = IntRadpSum + IntRad[iphip, ii]
        end

        IntSRnL        = - IntSR[iphin, bregionJ1]; IntSRnR  = - IntSR[iphin, bregionJ2]
        IntSRpL        = IntSR[iphip, bregionJ1]; IntSRpR  = IntSR[iphip, bregionJ2]

        push!(IRevSRHn, IntSRHnSum); push!(IRevSRHp, IntSRHpSum)
        push!(IRevRadn, IntRadnSum); push!(IRevRadp, IntRadpSum)
        push!(IRevSRnL, IntSRnL);    push!(IRevSRnR, IntSRnR)
        push!(IRevSRpL, IntSRpL);    push!(IRevSRpR, IntSRpR)

    end

    if printText
        println("*** done\n")
    end

    ################################################################################
    if printText
        println("Forward IV Measurement loop")
    end
    ################################################################################

    solForw = ChargeTransport.solve(ctsys, inival = solRev.u[end], times=(tPrecond + tend, tPrecond + 2*tend),   control = control)

    subg = subgrid(grid, [regionPero])
    ## calculate the integral of each carrier for each region
    intncc = ChargeTransport.integrate(ctsys, storage!, solForw.u[end])./q

    ## calculate the measure of the perovskite region
    mOmega = 0.0
    for icellVol in subg[CellVolumes]
        mOmega = mOmega + icellVol
    end

    if printText
        println("Average vacancy density ", intncc[iphia, regionPero]/mOmega)
    end

    if plotting

        if gridDim == 1
            figure()
            plot_densities(PyPlot, ctsys, solForw[end], "End time", label_density)
            figure()
            plot_solution(PyPlot, ctsys,  solForw[end], "End time", label_solution)
        end

    end

    if printText
        println("*** done\n")
    end

    ################################################################################
    if printText
        println("Forward IV curve calculation")
    end
    ################################################################################

    IV            = zeros(0) # for saving I-V data
    IVn           = zeros(0); IVp  = zeros(0)
    IVa           = zeros(0); IVψ  = zeros(0)

    ######################
    ISRHn         = zeros(0); IRadn = zeros(0)
    ISRHp         = zeros(0); IRadp = zeros(0)
    ISRnL         = zeros(0); ISRnR = zeros(0)
    ISRpL         = zeros(0); ISRpR = zeros(0)
    IGen          = zeros(0);


    tvalues       = solForw.t
    number_tsteps = length(tvalues)
    biasValues    = scanProtocol.(tvalues[2:end])

    factory       = TestFunctionFactory(ctsys)
    tf            = testfunction(factory, [bregionLeft], [bregionRight])

    for istep = 2:number_tsteps

        Δt       = tvalues[istep] - tvalues[istep-1] # Time step size
        inival   = solForw.u[istep-1]
        solution = solForw.u[istep]

        I        = integrate(ctsys, tf, solution, inival, Δt)

        current = 0.0
        for ii = 1:numberOfCarriers2+1
            current = current + I[ii]
        end

        push!(IV,   current); push!(IVn, I[iphin])
        push!(IVp, I[iphip]); push!(IVψ, I[ipsi])

        if enableIons
            push!(IVa, I[iphia])
        end

        IntSRH         = integrate(ctsys, SRHRecombination!, solution)
        IntRad         = integrate(ctsys, RadiativeRecombination!, solution)
        IntSR          = integrate(ctsys, SRRecombination!, solution, boundary = true)
        IntGen         = integrate(ctsys, Photogeneration!, solution)

        IntSRHnSum     = 0.0; IntRadnSum = 0.0
        IntSRHpSum     = 0.0; IntRadpSum = 0.0

        for ii = 1:numberOfRegions
            IntSRHnSum = IntSRHnSum - IntSRH[iphin, ii]
            IntRadnSum = IntRadnSum - IntRad[iphin, ii]

            IntSRHpSum = IntSRHpSum + IntSRH[iphip, ii]
            IntRadpSum = IntRadpSum + IntRad[iphip, ii]
        end

        IntGenp        = IntGen[iphip, regionPero]
        IntSRnL        = - IntSR[iphin, bregionJ1]; IntSRnR  = - IntSR[iphin, bregionJ2]
        IntSRpL        = IntSR[iphip, bregionJ1]; IntSRpR  = IntSR[iphip, bregionJ2]

        push!(ISRHn, IntSRHnSum); push!(ISRHp, IntSRHpSum)
        push!(IRadn, IntRadnSum); push!(IRadp, IntRadpSum)
        push!(ISRnL, IntSRnL);    push!(ISRnR, IntSRnR)
        push!(ISRpL, IntSRpL);    push!(ISRpR, IntSRpR)
        push!(IGen, IntGenp)

    end

    if printText
        println("*** done\n")
    end

    ################################################################################
    if printText
        println("Plotting and saving")
    end
    ################################################################################

    if plotting && gridDim == 1

        figure()
        plot(biasValues,    -IV.*(cm^2).*1.0e3,    linewidth = 5, color = "blue",  label ="forward"
        )
        plot(biasValuesRev, -IVRev.*(cm^2).*1.0e3, linewidth = 5, color = "green", linestyle = "--", label ="reverse")

        PyPlot.grid()
        PyPlot.legend()
        PyPlot.xlabel("bias [V]", fontsize=17)
        PyPlot.ylabel("current density [mAcm\$^{-2} \$]", fontsize=17)
        PyPlot.tick_params(which ="both", labelsize=18)
        PyPlot.tight_layout()

        #####################################
        figure()
        semilogy(biasValues, ISRHn.*(cm^2).*1.0e3, linewidth = 5, color = "blue",  label ="SRH")
        semilogy(biasValues, IRadn.*(cm^2).*1.0e3, linewidth = 5, color = "red",   label ="rad")
        semilogy(biasValues, IGen.*(cm^2).*1.0e3, linewidth = 5, color = "gold",   label ="Gen")

        semilogy(biasValues, ISRnL.*(cm^2).*1.0e3, linewidth = 5, color = "black", label = "SR, n, L")
        semilogy(biasValues, ISRpL.*(cm^2).*1.0e3, linewidth = 5, color = "gray",  linestyle=":",  label = "SR, p, L")

        semilogy(biasValues, ISRnR.*(cm^2).*1.0e3, linewidth = 5, color = "darkgreen",  label = "SR, n, R")
        semilogy(biasValues, ISRpR.*(cm^2).*1.0e3, linewidth = 5, color = "green",  linestyle=":", label = "SR, p, R")

        PyPlot.grid()
        PyPlot.legend()
        PyPlot.xlabel("bias [V]", fontsize=17)
        PyPlot.ylabel("current density [mAcm\$^{-2} \$]", fontsize=17)
        PyPlot.tick_params(which ="both", labelsize=18)
        PyPlot.tight_layout()
    end

    if generation
        IV                = -IV
        bias              = biasValues

        powerDensity      = bias .* (IV)           # power density function
        MaxPD, indexPD    = findmax(powerDensity)

        open_circuit      = compute_open_circuit_voltage(bias, IV)

        IncLightPowerDens = 1000.0 * W/m^2

        fillfactor        = 100 * (bias[indexPD] * IV[indexPD]) / (IV[1] * open_circuit)

        efficiency        = 100 * bias[indexPD] * (IV[indexPD]) / (IncLightPowerDens)

        if gridDim == 1
            JSC = IV[1].*(0.01)^(2).*1.0e3
        elseif gridDim == 2
            JSC = IV[1]./heightDev.*(0.01)^(2).*1.0e3
        end

        if printText
            println(" ")
            println("The JSC                  is $JSC mAcm^{-2}.")
            println("The fill factor          is $fillfactor %.")
            println("The efficiency           is $efficiency %.")
            println("The open circuit voltage is $open_circuit V.")
            println(" ")
        end

        IV = -IV

        # For reverse:

        IVR               = reverse(-IVRev)
        bias              = reverse(biasValuesRev)

        powerDensity      = bias .* (IVR)           # power density function
        MaxPD, indexPD    = findmax(powerDensity)

        open_circuit      = compute_open_circuit_voltage(bias, IVR)

        IncLightPowerDens = 1000.0 * W/m^2

        fillfactor        = 100 * (bias[indexPD] * IVR[indexPD]) / (IVR[1] * open_circuit)

        efficiency        = 100 * bias[indexPD] * (IVR[indexPD]) / (IncLightPowerDens)

        if gridDim == 1
            JSC = IVR[1].*(0.01)^(2).*1.0e3
        elseif gridDim == 2
            JSC = IVR[1]./heightDev.*(0.01)^(2).*1.0e3
        end

        if printText
            println(" ")
            println("The JSC                  is $JSC mAcm^{-2}.")
            println("The fill factor          is $fillfactor %.")
            println("The efficiency           is $efficiency %.")
            println("The open circuit voltage is $open_circuit V.")
            println(" ")
        end

    end

    if saveData

        helpampl = collect(string(amplitude));  helpampl[ findall(x -> x == '.', helpampl)[1] ] = 'p'
        textampl = join(helpampl)

        helpCa   = collect(string(Ca));  helpCa[ findall(x -> x == '.', helpCa)[1] ] = 'p'
        textCa   = join(helpCa)

        if typeGrid == "nanotextured"
            typeGridText = "nanotextured-ampl-$textampl"
        elseif typeGrid == "planar"
            typeGridText = "planar"
        end

        if generation

            openCircuitRev       = compute_open_circuit_voltage(reverse(biasValuesRev), reverse(IVRev))

            gRev1(t)             = scanProtocol(t) - openCircuitRev
            tOpenCircuitRev      = find_zero(gRev1, tPrecond)

            powerDensityRev      = - biasValuesRev .* IVRev    # power density function
            MaxPDRev, indexPDRev = findmax(powerDensityRev)
            gRev2(t)             = scanProtocol(t) - biasValuesRev[indexPDRev]
            tMPPRev              = find_zero(gRev2, tPrecond)

            ## reverse solution
            if enableIons
                if gridDim == 1
                    writedlm(datadir("sol", "parameter-$paramsname/$textSR/Sol-$(gridDim)D-rev-$(typeGridText)-Ca-$textCa-generation-$(generationType)-reco-$(typeReco)-inival.dat"), solRev(tPrecond)')
                    writedlm(datadir("sol", "parameter-$paramsname/$textSR/Sol-$(gridDim)D-rev-$(typeGridText)-Ca-$textCa-generation-$(generationType)-reco-$(typeReco)-MPP.dat"),    solRev(tMPPRev)')
                    writedlm(datadir("sol", "parameter-$paramsname/$textSR/Sol-$(gridDim)D-rev-$(typeGridText)-Ca-$textCa-generation-$(generationType)-reco-$(typeReco)-OCV.dat"),    solRev(tOpenCircuitRev)')
                else
                    writedlm(datadir("sol", "parameter-$paramsname/$textSR/Sol-$(gridDim)D-rev-$(typeGridText)-generation-$(generationType)-reco-$(typeReco)-inival.dat"), solRev(tPrecond)')
                    writedlm(datadir("sol", "parameter-$paramsname/$textSR/Sol-$(gridDim)D-rev-$(typeGridText)-generation-$(generationType)-reco-$(typeReco)-MPP.dat"),    solRev(tMPPRev)')
                    writedlm(datadir("sol", "parameter-$paramsname/$textSR/Sol-$(gridDim)D-rev-$(typeGridText)-generation-$(generationType)-reco-$(typeReco)-OCV.dat"),    solRev(tOpenCircuitRev)')
                end
            else
                writedlm(datadir("sol", "parameter-$paramsname/$textSR/Sol-$(gridDim)D-rev-$(typeGridText)-generation-$(generationType)-reco-$(typeReco)-inival-enableIons-false.dat"), solRev(tPrecond)')
                writedlm(datadir("sol", "parameter-$paramsname/$textSR/Sol-$(gridDim)D-rev-$(typeGridText)-generation-$(generationType)-reco-$(typeReco)-MPP-enableIons-false.dat"),    solRev(tMPPRev)')
                writedlm(datadir("sol", "parameter-$paramsname/$textSR/Sol-$(gridDim)D-rev-$(typeGridText)-generation-$(generationType)-reco-$(typeReco)-OCV-enableIons-false.dat"),    solRev(tOpenCircuitRev)')
            end

            ####################
            openCircuit         = compute_open_circuit_voltage(biasValues, IV)
            g1(t)               = scanProtocol(t) - openCircuit

            tOpenCircuit        = find_zero(g1, tPrecond+1.5*tend)

            powerDensity        = - biasValues .* IV    # power density function
            MaxPD, indexPD      = findmax(powerDensity)
            g2(t)               = scanProtocol(t) - biasValues[indexPD]
            tMPP                = find_zero(g2, tPrecond+1.5*tend)

            if enableIons
                if gridDim == 1
                    writedlm(datadir("sol", "parameter-$paramsname/$textSR/Sol-$(gridDim)D-forw-$(typeGridText)-Ca-$textCa-generation-$(generationType)-reco-$(typeReco)-inival.dat"), solForw(tPrecond+tend)')
                    writedlm(datadir("sol", "parameter-$paramsname/$textSR/Sol-$(gridDim)D-forw-$(typeGridText)-Ca-$textCa-generation-$(generationType)-reco-$(typeReco)-MPP.dat"),    solForw(tMPP)')
                    writedlm(datadir("sol", "parameter-$paramsname/$textSR/Sol-$(gridDim)D-forw-$(typeGridText)-Ca-$textCa-generation-$(generationType)-reco-$(typeReco)-OCV.dat"),    solForw(tOpenCircuit)')
                    writedlm(datadir("sol", "parameter-$paramsname/$textSR/Sol-$(gridDim)D-forw-$(typeGridText)-Ca-$textCa-generation-$(generationType)-reco-$(typeReco)-end.dat"),    solForw(tPrecond+2*tend)')
                else
                    writedlm(datadir("sol", "parameter-$paramsname/$textSR/Sol-$(gridDim)D-forw-$(typeGridText)-generation-$(generationType)-reco-$(typeReco)-inival.dat"), solForw(tPrecond+tend)')
                    writedlm(datadir("sol", "parameter-$paramsname/$textSR/Sol-$(gridDim)D-forw-$(typeGridText)-generation-$(generationType)-reco-$(typeReco)-MPP.dat"),    solForw(tMPP)')
                    writedlm(datadir("sol", "parameter-$paramsname/$textSR/Sol-$(gridDim)D-forw-$(typeGridText)-generation-$(generationType)-reco-$(typeReco)-OCV.dat"),    solForw(tOpenCircuit)')
                    writedlm(datadir("sol", "parameter-$paramsname/$textSR/Sol-$(gridDim)D-forw-$(typeGridText)-generation-$(generationType)-reco-$(typeReco)-end.dat"),    solForw(tPrecond+2*tend)')
                end
            else
                writedlm(datadir("sol", "parameter-$paramsname/$textSR/Sol-$(gridDim)D-forw-$(typeGridText)-generation-$(generationType)-reco-$(typeReco)-inival-enableIons-false.dat"), solForw(tPrecond+tend)')
                writedlm(datadir("sol", "parameter-$paramsname/$textSR/Sol-$(gridDim)D-forw-$(typeGridText)-generation-$(generationType)-reco-$(typeReco)-MPP-enableIons-false.dat"),    solForw(tMPP)')
                writedlm(datadir("sol", "parameter-$paramsname/$textSR/Sol-$(gridDim)D-forw-$(typeGridText)-generation-$(generationType)-reco-$(typeReco)-OCV-enableIons-false.dat"),    solForw(tOpenCircuit)')
                writedlm(datadir("sol", "parameter-$paramsname/$textSR/Sol-$(gridDim)D-forw-$(typeGridText)-generation-$(generationType)-reco-$(typeReco)-end-enableIons-false.dat"),    solForw(tPrecond+2*tend)')
            end

        end

        ## reverse I-V
        if enableIons
            if Ca == 6.0e22 # dont save when average vacancy dens is not equal to the original val
                writedlm(datadir("IV", "parameter-$paramsname/$textSR/IV-$(gridDim)D-rev-$(typeGridText)-generation-$(generationType)-reco-$(typeReco).dat"), [biasValuesRev IVRev])
                ############
                writedlm(datadir("IV", "parameter-$paramsname/$textSR/JSRH-$(gridDim)D-rev-$(typeGridText)-generation-$(generationType)-reco-$(typeReco).dat"),  [biasValuesRev IRevSRHn])
                writedlm(datadir("IV", "parameter-$paramsname/$textSR/JRad-$(gridDim)D-rev-$(typeGridText)-generation-$(generationType)-reco-$(typeReco).dat"),  [biasValuesRev IRevRadn])
                writedlm(datadir("IV", "parameter-$paramsname/$textSR/JSRnL-$(gridDim)D-rev-$(typeGridText)-generation-$(generationType)-reco-$(typeReco).dat"), [biasValuesRev IRevSRnL])
                writedlm(datadir("IV", "parameter-$paramsname/$textSR/JSRpL-$(gridDim)D-rev-$(typeGridText)-generation-$(generationType)-reco-$(typeReco).dat"), [biasValuesRev IRevSRpL])
                writedlm(datadir("IV", "parameter-$paramsname/$textSR/JSRnR-$(gridDim)D-rev-$(typeGridText)-generation-$(generationType)-reco-$(typeReco).dat"), [biasValuesRev IRevSRnR])
                writedlm(datadir("IV", "parameter-$paramsname/$textSR/JSRpR-$(gridDim)D-rev-$(typeGridText)-generation-$(generationType)-reco-$(typeReco).dat"), [biasValuesRev IRevSRpR])
            end
        else
            writedlm(datadir("IV", "parameter-$paramsname/$textSR/IV-$(gridDim)D-rev-$(typeGridText)-generation-$(generationType)-reco-$(typeReco)-enableIons-false.dat"), [biasValuesRev IVRev])
        end

        ## forward I-V
        if enableIons
            if Ca == 6.0e22 # dont save when average vacancy dens is not equal to the original val
                writedlm(datadir("IV", "parameter-$paramsname/$textSR/IV-$(gridDim)D-forw-$(typeGridText)-generation-$(generationType)-reco-$(typeReco).dat"), [biasValues IV])
                ############
                writedlm(datadir("IV", "parameter-$paramsname/$textSR/JSRH-$(gridDim)D-forw-$(typeGridText)-generation-$(generationType)-reco-$(typeReco).dat"),  [biasValues ISRHn])
                writedlm(datadir("IV", "parameter-$paramsname/$textSR/JRad-$(gridDim)D-forw-$(typeGridText)-generation-$(generationType)-reco-$(typeReco).dat"),  [biasValues IRadn])
                writedlm(datadir("IV", "parameter-$paramsname/$textSR/JGen-$(gridDim)D-forw-$(typeGridText)-generation-$(generationType)-reco-$(typeReco).dat"),  [biasValues IGen])
                writedlm(datadir("IV", "parameter-$paramsname/$textSR/JSRnL-$(gridDim)D-forw-$(typeGridText)-generation-$(generationType)-reco-$(typeReco).dat"), [biasValues ISRnL])
                writedlm(datadir("IV", "parameter-$paramsname/$textSR/JSRpL-$(gridDim)D-forw-$(typeGridText)-generation-$(generationType)-reco-$(typeReco).dat"), [biasValues ISRpL])
                writedlm(datadir("IV", "parameter-$paramsname/$textSR/JSRnR-$(gridDim)D-forw-$(typeGridText)-generation-$(generationType)-reco-$(typeReco).dat"), [biasValues ISRnR])
                writedlm(datadir("IV", "parameter-$paramsname/$textSR/JSRpR-$(gridDim)D-forw-$(typeGridText)-generation-$(generationType)-reco-$(typeReco).dat"), [biasValues ISRpR])
            end

        else
            writedlm(datadir("IV", "parameter-$paramsname/$textSR/IV-$(gridDim)D-forw-$(typeGridText)-generation-$(generationType)-reco-$(typeReco)-enableIons-false.dat"), [biasValues IV])
        end

    end

    if printText
        println("*** done\n")
    end

    testval = sum(filter(!isnan, sol0p0))/length(sol0p0) # when using sparse storage, we get NaN values in solution
    return testval

end #  main

function test(;gridDim=1, typeGrid = "planar", amplitude = 2.0e-7, demo_run = false)
    if demo_run
        if gridDim == 1
            testval = -0.6236278584408826 # all reco
        elseif gridDim == 2
            if typeGrid == "planar"
                testval = -1.1597613291016444 # all reco
            elseif typeGrid == "nanotextured" && amplitude == 2.0e-7
                testval = -1.1880359433451766 # all reco
            end
        end
    else
        if gridDim == 1
            testval = -0.6162695397900051 # all reco
        end
    end
    result = main(gridDim = gridDim, typeGrid = typeGrid, amplitude = amplitude, printText = false, demo_run = demo_run, generation = true, generationUniform = false, MaxwellSol = true)
    @info "result  = $result"
    @info "testval = $testval"
    return abs(result - testval) < 1e-15
end

end # module


# demo_run: false (Average vacancy density 5.985820758704952e22)
# The JSC                  is 21.401388798480134 mAcm^{-2}.
# The fill factor          is 83.16990537331397 %.
# The efficiency           is 20.858910215594932 %.
# The open circuit voltage is 1.1718808313367537 V.


# Reverse curve:
# The JSC                  is 21.487861289295004 mAcm^{-2}.
# The fill factor          is 83.6435670039751 %.
# The efficiency           is 21.04063950025071 %.
# The open circuit voltage is 1.1706665209604834 V.

###################################################################
###################################################################

# demo_run: true (Average vacancy density 6.0123145755738964e22)
# The JSC                  is 21.402636816225165 mAcm^{-2}.
# The fill factor          is 83.22201406248647 %.
# The efficiency           is 20.857266019977082 %.
# The open circuit voltage is 1.170986468002918 V.


# Reverse curve:
# The JSC                  is 21.490317966799385 mAcm^{-2}.
# The fill factor          is 83.68497373607075 %.
# The efficiency           is 21.043217928526097 %.
# The open circuit voltage is 1.170096896416928 V.