#=

Simulation of a three-layer perovskite solar cell in 2D with textured interface.
We consider the composition C60 -- triple cation -- PTAA
The parameters follow: https://doi.org/10.25446/oxford.24359959
from the publication Thiesbrummel et al., nature energy, 2024.

=#

ENV["LC_NUMERIC"] = "C"

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

    include("parameter.jl")

    # for convenience
    datadir = TexturedPerovskiteSolarCells.datadir
    scriptsdir = TexturedPerovskiteSolarCells.scriptsdir

    ###########################################################
    ###########################################################

    function main(;
            plotting = false, printText = true,
            ########################
            parameter_set = ParamsSingleJunction,
            ########################
            gridDim = 1, typeGrid = "planar",  # "nanotextured", #
            ########################
            amplitude = 4.0e-7,    # amplitude of nanotexture in base unit
            ########################
            typeReco = "all",    # "radiative", #"off", # "radiative", # "SR", # "bulk", # "all"   (SR for surface reco)
            ########################
            generation = true, generationUniform = false, MaxwellSol = true,
            ########################
            scanrate = 1.0e3, # 1.0e-3, # 1.0e1, # in base unit
            ########################
            vacancyEnergyCalculation = false, # for calculating Ea of φa
            ########################
            enableIons = true, saveData = false,
            ########################
            demo_run = false,    # do calculations on coarser mesh, if true
            ########################
            plotPostProcess = false,
            ########################
        )

        PyPlot.rc("font", family = "sans-serif", size = 14)
        PyPlot.rc("mathtext", fontset = "dejavusans")
        PyPlot.close("all")

        @local_unitfactors V cm m s W

        (; q, ε_0) = ChargeTransport.constants
        eV = q * V

        # we need this, since methods coincide with VoronoiFVM methods
        SolverControl() = ChargeTransport.SolverControl()
        unknowns(a) = ChargeTransport.unknowns(a)
        TestFunctionFactory(a) = ChargeTransport.TestFunctionFactory(a)
        testfunction(a, b, c) = ChargeTransport.testfunction(a, b, c)
        integrate(a, b, c, d, e; kwarges...) = ChargeTransport.integrate(a, b, c, d, e; kwarges...)
        integrate(a, b, c; kwarges...) = ChargeTransport.integrate(a, b, c; kwarges...)

        if gridDim == 1
            typeGrid = "planar"
        end

        ################################################################################
        if printText
            println("Set up grid and regions")
        end
        ################################################################################

        # use the destructuring operator to extract all the necessary parameters
        (;
            Ea1D, Ca, numberOfCarriers, bregionJ1, bregionJ2, iphin, iphip, bregionLeft, bregionRight, heightLayersPL, Fcc, regionPero, iphia, T,
            zn, zp, za, numberOfRegions, εr, Nn, Np, Na, En, Ep, ΔEFLeft, ΔEFRight, μn, μp, μa, r0, τn, τp, nτ, pτ, regionETL, regionHTL, Cp, Ca, Cn, damp_initial, damp_growth, maxiters, max_round, abstol, reltol,
            tol_round, EaPlanar, EaAmpl0p5e7, EaAmpl1p0e7, EaAmpl1p5e7, EaAmpl2p0e7, EaAmpl2p5e7, EaAmpl3p0e7, EaAmpl3p5e7, EaAmpl4p0e7, EaAmpl4p5e7,
            EaAmpl5p0e7, EaAmpl5p5e7, EaAmpl6p0e7, EaAmpl6p5e7, EaAmpl7p0e7, EaAmpl7p5e7, incidentPhotonFlux, absorption, generationPeak, invertedIllumination,
            heightDev,
        ) = parameter_set()

        grid = generate_grid(gridDim = gridDim, type = typeGrid, amplitude = amplitude, parameter_set = parameter_set, demo_run = demo_run)

        if plotting
            gridplot(grid, Plotter = PyPlot, resolution = (600, 400), linewidth = 0.5, legend = :rc)
        end

        if printText
            println("*** done\n")
        end

        ################################################################################
        if printText
            println("Define physical parameters and model")
        end
        ################################################################################

        ## recombination velocities
        # ## correct values:
        SRHvelocityETLn = 2000 * cm / s
        SRHvelocityETLp = 2000 * cm / s

        SRHvelocityHTLn = 500 * cm / s
        SRHvelocityHTLp = 500 * cm / s

        ########## primary data for I-V scan protocol ##########
        if typeReco == "all" # radiative until 1.4 V and reco = all until 1.2 V
            endVoltage = 1.2 * V
        else
            endVoltage = 1.4 * V
        end

        if SRHvelocityETLn <= 200 * cm / s
            endVoltage = 1.27 * V
        end

        if SRHvelocityETLn == 2000 * cm / s && SRHvelocityHTLn <= 100 * cm / s
            endVoltage = 1.2 * V
        end

        tPrecond = 0.001
        tend = endVoltage / scanrate

        ## Define scan protocol function
        function scanProtocolPrecond(t)
            if data.calculationType == InEquilibrium
                biasVal = 0.0
            else
                if 0.0 <= t && t <= tPrecond
                    biasVal = endVoltage
                elseif tPrecond <= t  && t <= tPrecond + tend
                    biasVal = endVoltage .- scanrate * (t - tPrecond)
                elseif tend + tPrecond < t  && t <= tPrecond + 2 * tend
                    biasVal = 0.0 .+ scanrate * (t - tend - tPrecond)
                else
                    biasVal = 0.0
                end
            end

            # we need this such that during the computation of the equilibrium solution, we are safe that no bias is applied.
            if !data.generationComplete
                biasVal = 0.0
            end
            return biasVal
        end

        ## Define scan protocol function
        function scanProtocolLinear(t)

            if 0.0 <= t && t <= tend
                biasVal = scanrate * t
            elseif tend <= t  && t <= 2 * tend
                biasVal = endVoltage .- scanrate * (t - tend)
            else
                biasVal = 0.0
            end

            # we need this such that during the computation of the equilibrium solution, we are safe that no bias is applied.
            if !data.generationComplete
                biasVal = 0.0
            end

            return biasVal
        end

        # Apply zero voltage on left boundary and a linear scan protocol on right boundary
        contactVoltageFunction = [zeroVoltage, scanProtocolPrecond]

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
                if amplitude == 0.5e-7
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
                end
            end
        end

        # For all other parameters, we refer to the parameters template
        if printText
            println("*** done\n")
        end

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

                generationData = MaxwellPhotogeneration(gridDim = gridDim, typeGrid = typeGrid, amplitude = amplitude, parameter_set = parameter_set, demo_run = demo_run)

            elseif MaxwellSol == false && generationUniform == false
                Fph = incidentPhotonFlux[regionPero]; ag = absorption[regionPero]
                genPeak = generationPeak;                 inv = invertedIllumination

                function BeerLamb(x)
                    return if heightLayersPL[1] <= x <= heightLayersPL[2]
                        Fph .* ag .* exp.(- inv .* ag .* (x .- genPeak))
                    else
                        0.0
                    end
                end

                if gridDim == 1
                    coord = vec(grid[Coordinates])
                elseif gridDim == 2
                    coord = vec(grid[Coordinates][2, :])
                end

                generationData = BeerLamb.(coord)

            elseif generationUniform

                generationData = zeros(length(grid[Coordinates][1, :]))
                subg = subgrid(grid, [regionPero])
                iNode = subg[NodeParents]

                generationData[iNode] .= 3.22e27 * 1 / (m^3 * s) # average Photogeneration for planar

            end

            data = Data(grid, numberOfCarriers2, contactVoltageFunction = contactVoltageFunction, generationData = generationData)
        else
            data = Data(grid, numberOfCarriers2, contactVoltageFunction = contactVoltageFunction)
        end

        data.modelType = Transient

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

        data.bulkRecombination = set_bulk_recombination(;
            iphin = iphin, iphip = iphip, bulk_recomb_Auger = false,
            bulk_recomb_radiative = bulk_recomb_radiative,
            bulk_recomb_SRH = bulk_recomb_SRH
        )

        data.boundaryType[bregionLeft] = MixedOhmicSchottkyContact
        data.boundaryType[bregionRight] = MixedOhmicSchottkyContact

        if enableIons
            data.F = Fcc
            enable_ionic_carrier!(data, ionicCarrier = iphia, regions = [regionPero])
        else
            data.F = Fcc[1:2]
        end

        if generation
            if MaxwellSol
                generationType = "Maxwell"
            else
                generationType = "Beer-Lambert"
            end
        else
            generationType = "none"
        end

        if generation && generationUniform
            generationType = "uniform"
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

        params = Params(grid[NumCellRegions], grid[NumBFaceRegions], numberOfCarriers2)

        params.temperature = T
        params.chargeNumbers[iphin] = zn
        params.chargeNumbers[iphip] = zp
        if enableIons
            params.chargeNumbers[iphia] = za
        end

        for ireg in 1:numberOfRegions # interior region data

            params.dielectricConstant[ireg] = εr[ireg] * ε_0

            ## effective DOS, band edge energy and mobilities
            params.densityOfStates[iphin, ireg] = Nn[ireg]
            params.densityOfStates[iphip, ireg] = Np[ireg]

            params.bandEdgeEnergy[iphin, ireg] = En[ireg]
            params.bandEdgeEnergy[iphip, ireg] = Ep[ireg]

            params.mobility[iphin, ireg] = μn[ireg]
            params.mobility[iphip, ireg] = μp[ireg]

            ## vacancy parameters
            if enableIons
                params.densityOfStates[iphia, ireg] = Na[ireg]
                params.bandEdgeEnergy[iphia, ireg] = Ea[ireg]
                params.mobility[iphia, ireg] = μa[ireg]

            end

            ## recombination parameters
            params.recombinationRadiative[ireg] = r0[ireg]
            params.recombinationSRHLifetime[iphin, ireg] = τn[ireg]
            params.recombinationSRHLifetime[iphip, ireg] = τp[ireg]
            params.recombinationSRHTrapDensity[iphin, ireg] = nτ[ireg]
            params.recombinationSRHTrapDensity[iphip, ireg] = pτ[ireg]

        end

        params.SchottkyBarrier[bregionLeft] = ΔEFLeft
        params.SchottkyBarrier[bregionRight] = En[3] - Ep[3] - ΔEFRight

        ##############################################################
        ## inner boundary region data (we choose the intrinsic values)
        if typeReco == "SR" || typeReco == "all"
            params.bDensityOfStates[iphin, bregionJ1] = Nn[regionETL]
            params.bDensityOfStates[iphip, bregionJ1] = Np[regionPero]

            params.bDensityOfStates[iphin, bregionJ2] = Nn[regionPero]
            params.bDensityOfStates[iphip, bregionJ2] = Np[regionHTL]

            params.bBandEdgeEnergy[iphin, bregionJ1] = En[regionETL]
            params.bBandEdgeEnergy[iphip, bregionJ1] = Ep[regionPero]

            params.bBandEdgeEnergy[iphin, bregionJ2] = En[regionPero]
            params.bBandEdgeEnergy[iphip, bregionJ2] = Ep[regionHTL]

            ## for surface recombination
            params.recombinationSRHvelocity[iphin, bregionJ1] = SRHvelocityETLn
            params.recombinationSRHvelocity[iphip, bregionJ1] = SRHvelocityETLp

            params.bRecombinationSRHTrapDensity[iphin, bregionJ1] = params.recombinationSRHTrapDensity[iphin, regionETL]
            params.bRecombinationSRHTrapDensity[iphip, bregionJ1] = params.recombinationSRHTrapDensity[iphip, regionPero]

            params.recombinationSRHvelocity[iphin, bregionJ2] = SRHvelocityHTLn
            params.recombinationSRHvelocity[iphip, bregionJ2] = SRHvelocityHTLp

            params.bRecombinationSRHTrapDensity[iphin, bregionJ2] = params.recombinationSRHTrapDensity[iphin, regionPero]
            params.bRecombinationSRHTrapDensity[iphip, bregionJ2] = params.recombinationSRHTrapDensity[iphip, regionHTL]
        end

        ##############################################################

        ## interior doping
        params.doping[iphip, regionHTL] = Cp
        params.doping[iphin, regionETL] = Cn

        # doping
        if enableIons
            params.doping[iphia, regionPero] = Ca
        end

        data.params = params
        ctsys = System(grid, data, unknown_storage = :sparse)

        ipsi = data.index_psi

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

        control = SolverControl()
        if printText
            control.verbose = "e"
        else
            control.verbose = false
        end
        control.method_linear = UMFPACKFactorization()
        control.damp_initial = damp_initial
        control.damp_growth = damp_growth
        control.max_round = max_round
        control.maxiters = maxiters

        control.abstol = abstol
        control.reltol = reltol
        control.tol_round = tol_round

        control.Δu_opt = Inf

        if printText
            println("*** done\n")
        end

        ################################################################################
        if printText
            println("Compute solution in thermodynamic equilibrium")
        end
        ################################################################################

        # initialize solution and starting vectors
        inival = unknowns(ctsys); initialCond = unknowns(ctsys)
        sol0p0 = unknowns(ctsys)
        inival .= 0.0;             initialCond .= 0.0

        solEQ = equilibrium_solve!(ctsys, control = control, vacancyEnergyCalculation = vacancyEnergyCalculation)

        inival = solEQ

        if plotting
            label_solution, label_density, label_energy, label_BEE = set_plotting_labels(data)

            if enableIons
                ## add labels for anion vacancy
                label_energy[1, iphia] = "\$E_a-q\\psi\$"; label_energy[2, iphia] = "\$ - q \\varphi_a\$"; label_BEE[iphia] = "\$E_a\$"
                label_density[iphia] = "\$ n_a \$";      label_solution[iphia] = "\$ \\varphi_a\$"
            end
        end

        integral = integrated_density(ctsys, sol = solEQ, icc = iphia, ireg = regionPero)
        mOmega = data.regionVolumes[regionPero]

        println("Calculated average vacancy density is: ", integral / mOmega)
        println(" ")
        vacancyEnergy = data.params.bandEdgeEnergy[iphia, regionPero] / q
        println("Value for vacancy energy is: ", vacancyEnergy, " eV. Save this value for later use.")
        println(" ")

        if demo_run && gridDim == 2 && generation && MaxwellSol
            if typeGrid == "planar"
                return testval = sum(filter(!isnan, inival)) / length(inival) # when using sparse storage, we get NaN values in solution
            elseif typeGrid == "nanotextured" && amplitude == 2.0e-7
                return testval = sum(filter(!isnan, inival)) / length(inival) # when using sparse storage, we get NaN values in solution
            end
        end

        if printText
            println("*** done\n")
        end

        ################################################################################
        if printText
            println("Loop for increasing bias (to start at V near VOC)")
        end
        ################################################################################

        control.max_round = 5
        control.abstol = abstol
        control.reltol = reltol
        control.tol_round = tol_round

        if typeGrid == "nanotextured" && gridDim == 2 && amplitude == 7.5e-7 && SRHvelocityETLn <= 50 * cm / s
            control.damp_initial = 0.5
            biasLoop = collect(reverse(-range(0.0, endVoltage, length = 31)))
        else
            control.damp_initial = 0.5
            biasLoop = collect(reverse(-range(0.0, endVoltage, length = 21)))
        end

        for istep in 1:length(biasLoop)

            ## turn slowly voltage on
            set_contact!(ctsys, bregionRight, Δu = biasLoop[istep])

            if printText
                println(" bias =", biasLoop[istep])
            end

            ## keep some big time stepping here, so we preserve the vacancy density
            sol0p0 = ChargeTransport.solve(ctsys, inival = inival, control = control, tstep = 1.0e3)
            inival = sol0p0

        end # bias loop

        if plotting

            if gridDim == 1
                figure()
                plot_densities(PyPlot, ctsys, sol0p0, "Initial condition", label_density)
                figure()
                plot_solution(PyPlot, ctsys, sol0p0, "Initial condition", label_solution)
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

        control.damp_initial = 0.5
        control.Δt = 1.0e-4
        control.Δt_min = 1.0e-5
        control.Δt_max = 1.0e-4
        control.Δt_grow = 1.03

        solPrecond = ChargeTransport.solve(ctsys, inival = inival, times = (0.0, tPrecond), control = control)

        if printText
            println("*** done\n")
        end

        ################################################################################
        if printText
            println("Reverse IV Measurement loop")
        end
        ################################################################################

        control.Δt = 6.0e-4 / scanrate
        control.Δt_min = 6.0e-3 / scanrate
        control.Δt_max = 8.0e-3 / scanrate
        control.Δt_grow = 1.03
        solRev = ChargeTransport.solve(ctsys, inival = solPrecond.u[end], times = (tPrecond, tPrecond + tend), control = control)

        if printText
            println("*** done\n")
        end

        ################################################################################
        if printText
            println("Forward IV Measurement loop")
        end
        ################################################################################

        solForw = ChargeTransport.solve(ctsys, inival = solRev.u[end], times = (tPrecond + tend, tPrecond + 2 * tend), control = control)

        if plotting

            if gridDim == 1
                figure()
                plot_densities(PyPlot, ctsys, solForw[end], "End time", label_density)
                figure()
                plot_solution(PyPlot, ctsys, solForw[end], "End time", label_solution)
            end

        end

        if printText
            println("*** done\n")
        end

        ################################################################################
        if printText
            println("Reverse IV curve calculation")
        end
        ################################################################################

        IVRev = zeros(0) # for saving I-V data
        IVRevn = zeros(0); IVRevp = zeros(0)
        IVReva = zeros(0); IVRevψ = zeros(0)

        ######################
        IRevSRHn = zeros(0); IRevRadn = zeros(0)
        IRevSRHp = zeros(0); IRevRadp = zeros(0)
        IRevSRnL = zeros(0); IRevSRnR = zeros(0)
        IRevSRpL = zeros(0); IRevSRpR = zeros(0)

        tvaluesRev = solRev.t
        number_tsteps = length(tvaluesRev)
        biasValuesRev = data.contactVoltageFunction[2].(tvaluesRev[2:end])

        factory = TestFunctionFactory(ctsys)
        tf = testfunction(factory, [bregionLeft], [bregionRight])

        for istep in 2:number_tsteps

            Δt = tvaluesRev[istep] - tvaluesRev[istep - 1] # Time step size
            inival = solRev.u[istep - 1]
            solution = solRev.u[istep]

            I = integrate(ctsys, tf, solution, inival, Δt)

            current = 0.0
            for ii in 1:(numberOfCarriers2 + 1)
                current = current + I[ii]
            end

            push!(IVRev, current); push!(IVRevn, I[iphin])
            push!(IVRevp, I[iphip]); push!(IVRevψ, I[ipsi])

            if enableIons
                push!(IVReva, I[iphia])
            end

            IntSRH = integrate(ctsys, SRHRecombination!, solution)
            IntRad = integrate(ctsys, RadiativeRecombination!, solution)
            IntSR = integrate(ctsys, SRRecombination!, solution, boundary = true)

            IntSRHnSum = 0.0; IntRadnSum = 0.0
            IntSRHpSum = 0.0; IntRadpSum = 0.0

            for ii in 1:numberOfRegions
                IntSRHnSum = IntSRHnSum - IntSRH[iphin, ii]
                IntRadnSum = IntRadnSum - IntRad[iphin, ii]

                IntSRHpSum = IntSRHpSum + IntSRH[iphip, ii]
                IntRadpSum = IntRadpSum + IntRad[iphip, ii]
            end

            IntSRnL = - IntSR[iphin, bregionJ1]; IntSRnR = - IntSR[iphin, bregionJ2]
            IntSRpL = IntSR[iphip, bregionJ1]; IntSRpR = IntSR[iphip, bregionJ2]

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
            println("Forward IV curve calculation")
        end
        ################################################################################

        IV = zeros(0) # for saving I-V data
        IVn = zeros(0); IVp = zeros(0)
        IVa = zeros(0); IVψ = zeros(0)

        ######################
        ISRHn = zeros(0); IRadn = zeros(0)
        ISRHp = zeros(0); IRadp = zeros(0)
        ISRnL = zeros(0); ISRnR = zeros(0)
        ISRpL = zeros(0); ISRpR = zeros(0)
        IGen = zeros(0)

        tvalues = solForw.t
        number_tsteps = length(tvalues)
        biasValues = data.contactVoltageFunction[2].(tvalues[2:end])

        factory = TestFunctionFactory(ctsys)
        tf = testfunction(factory, [bregionLeft], [bregionRight])

        for istep in 2:number_tsteps

            Δt = tvalues[istep] - tvalues[istep - 1] # Time step size
            inival = solForw.u[istep - 1]
            solution = solForw.u[istep]

            I = integrate(ctsys, tf, solution, inival, Δt)

            current = 0.0
            for ii in 1:(numberOfCarriers2 + 1)
                current = current + I[ii]
            end

            push!(IV, current); push!(IVn, I[iphin])
            push!(IVp, I[iphip]); push!(IVψ, I[ipsi])

            if enableIons
                push!(IVa, I[iphia])
            end

            IntSRH = integrate(ctsys, SRHRecombination!, solution)
            IntRad = integrate(ctsys, RadiativeRecombination!, solution)
            IntSR = integrate(ctsys, SRRecombination!, solution, boundary = true)
            IntGen = integrate(ctsys, Photogeneration!, solution)

            IntSRHnSum = 0.0; IntRadnSum = 0.0
            IntSRHpSum = 0.0; IntRadpSum = 0.0

            for ii in 1:numberOfRegions
                IntSRHnSum = IntSRHnSum - IntSRH[iphin, ii]
                IntRadnSum = IntRadnSum - IntRad[iphin, ii]

                IntSRHpSum = IntSRHpSum + IntSRH[iphip, ii]
                IntRadpSum = IntRadpSum + IntRad[iphip, ii]
            end

            IntGenp = IntGen[iphip, regionPero]
            IntSRnL = - IntSR[iphin, bregionJ1]; IntSRnR = - IntSR[iphin, bregionJ2]
            IntSRpL = IntSR[iphip, bregionJ1]; IntSRpR = IntSR[iphip, bregionJ2]

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
            tEnd = tPrecond + 2 * tend

            tt = collect(range(0.0, tEnd, length = 201))
            T = data.contactVoltageFunction[2]
            plot(tt, T.(tt), marker = "o")
            xlabel("time [s]")
            ylabel("voltage [V]")
            PyPlot.tight_layout()

            figure()
            plot(biasValues, -IV .* (cm^2) .* 1.0e3, linewidth = 5, color = "blue", label = "forward")
            plot(biasValuesRev, -IVRev .* (cm^2) .* 1.0e3, linewidth = 5, color = "green", linestyle = "--", label = "reverse")

            PyPlot.grid()
            PyPlot.legend()
            PyPlot.xlabel("bias [V]", fontsize = 17)
            PyPlot.ylabel("current density [mAcm\$^{-2} \$]", fontsize = 17)
            PyPlot.tick_params(which = "both", labelsize = 18)
            PyPlot.tight_layout()

            #####################################
            figure()
            semilogy(biasValues, ISRHn .* (cm^2) .* 1.0e3, linewidth = 5, color = "blue", label = "SRH")
            semilogy(biasValues, IRadn .* (cm^2) .* 1.0e3, linewidth = 5, color = "red", label = "rad")
            semilogy(biasValues, IGen .* (cm^2) .* 1.0e3, linewidth = 5, color = "gold", label = "Gen")

            semilogy(biasValues, ISRnL .* (cm^2) .* 1.0e3, linewidth = 5, color = "black", label = "SR, n, L")
            semilogy(biasValues, ISRpL .* (cm^2) .* 1.0e3, linewidth = 5, color = "gray", linestyle = ":", label = "SR, p, L")

            semilogy(biasValues, ISRnR .* (cm^2) .* 1.0e3, linewidth = 5, color = "darkgreen", label = "SR, n, R")
            semilogy(biasValues, ISRpR .* (cm^2) .* 1.0e3, linewidth = 5, color = "green", linestyle = ":", label = "SR, p, R")

            PyPlot.grid()
            PyPlot.legend()
            PyPlot.xlabel("bias [V]", fontsize = 17)
            PyPlot.ylabel("current density [mAcm\$^{-2} \$]", fontsize = 17)
            PyPlot.tick_params(which = "both", labelsize = 18)
            PyPlot.tight_layout()
        end

        if generation
            IV = -IV
            bias = biasValues

            powerDensity = bias .* (IV)           # power density function
            MaxPD, indexPD = findmax(powerDensity)

            open_circuit = compute_open_circuit_voltage(bias, IV)

            IncLightPowerDens = 1000.0 * W / m^2

            fillfactor = 100 * (bias[indexPD] * IV[indexPD]) / (IV[1] * open_circuit)

            if gridDim == 1
                efficiency = 100 * bias[indexPD] * (IV[indexPD]) / (IncLightPowerDens)
                JSC = IV[1] .* (0.01)^(2) .* 1.0e3
            elseif gridDim == 2
                efficiency = 100 * bias[indexPD] * (IV[indexPD] ./ heightDev) / (IncLightPowerDens)
                JSC = IV[1] ./ heightDev .* (0.01)^(2) .* 1.0e3
            end

            if printText
                println(" ")
                println("The JSC                  is $(round(JSC, digits = 6)) mAcm^{-2}.")
                println("The fill factor          is $(round(fillfactor, digits = 6)) %.")
                println("The efficiency           is $(round(efficiency, digits = 6)) %.")
                println("The open circuit voltage is $(round(open_circuit, digits = 6)) V.")
                println(" ")
            end

            IV = -IV

            # For reverse:
            IVR = reverse(-IVRev)
            bias = reverse(biasValuesRev)

            powerDensity = bias .* (IVR)           # power density function
            MaxPD, indexPD = findmax(powerDensity)

            open_circuit = compute_open_circuit_voltage(bias, IVR)

            IncLightPowerDens = 1000.0 * W / m^2

            fillfactor = 100 * (bias[indexPD] * IVR[indexPD]) / (IVR[1] * open_circuit)


            if gridDim == 1
                efficiency = 100 * bias[indexPD] * (IVR[indexPD]) / (IncLightPowerDens)
                JSC = IVR[1] .* (0.01)^(2) .* 1.0e3
            elseif gridDim == 2
                efficiency = 100 * bias[indexPD] * (IVR[indexPD] ./ heightDev) / (IncLightPowerDens)
                JSC = IVR[1] ./ heightDev .* (0.01)^(2) .* 1.0e3
            end

            if printText
                println(" ")
                println("The JSC                  is $(round(JSC, digits = 6)) mAcm^{-2}.")
                println("The fill factor          is $(round(fillfactor, digits = 6)) %.")
                println("The efficiency           is $(round(efficiency, digits = 6)) %.")
                println("The open circuit voltage is $(round(open_circuit, digits = 6)) V.")
                println(" ")
            end

        end

        if saveData

            helpampl = collect(string(amplitude));  helpampl[findall(x -> x == '.', helpampl)[1]] = 'p'
            textampl = join(helpampl)

            helpCa = collect(string(Ca));  helpCa[findall(x -> x == '.', helpCa)[1]] = 'p'
            textCa = join(helpCa)

            if typeGrid == "nanotextured"
                typeGridText = "nanotextured-ampl-$textampl"
            elseif typeGrid == "planar"
                typeGridText = "planar"
            end

            if generation

                ## reverse scan
                tGuess = tPrecond

                powerDensityRev = - biasValuesRev .* IVRev    # power density function
                MaxPDRev, indexPDRev = findmax(powerDensityRev)
                gRev2(t) = data.contactVoltageFunction[2](t) - biasValuesRev[indexPDRev]
                tMPPRev = find_zero(gRev2, tGuess)

                inivalRev = tPrecond
                endRev = tPrecond + tend

                if enableIons
                    if gridDim == 1
                        writedlm(datadir("sol", "Sol-$(gridDim)D-rev-$(typeGridText)-Ca-$textCa-generation-$(generationType)-reco-$(typeReco)-inival.dat"), solRev(inivalRev)')
                        writedlm(datadir("sol", "Sol-$(gridDim)D-rev-$(typeGridText)-Ca-$textCa-generation-$(generationType)-reco-$(typeReco)-MPP.dat"), solRev(tMPPRev)')
                        writedlm(datadir("sol", "Sol-$(gridDim)D-rev-$(typeGridText)-Ca-$textCa-generation-$(generationType)-reco-$(typeReco)-end.dat"), solRev(endRev)')
                    else
                        writedlm(datadir("sol", "Sol-$(gridDim)D-rev-$(typeGridText)-generation-$(generationType)-reco-$(typeReco)-inival.dat"), solRev(inivalRev)')
                        writedlm(datadir("sol", "Sol-$(gridDim)D-rev-$(typeGridText)-generation-$(generationType)-reco-$(typeReco)-MPP.dat"), solRev(tMPPRev)')
                        writedlm(datadir("sol", "Sol-$(gridDim)D-rev-$(typeGridText)-generation-$(generationType)-reco-$(typeReco)-end.dat"), solRev(endRev)')
                    end
                else
                    writedlm(datadir("sol", "Sol-$(gridDim)D-rev-$(typeGridText)-generation-$(generationType)-reco-$(typeReco)-inival-enableIons-false.dat"), solRev(inivalRev)')
                    writedlm(datadir("sol", "Sol-$(gridDim)D-rev-$(typeGridText)-generation-$(generationType)-reco-$(typeReco)-MPP-enableIons-false.dat"), solRev(tMPPRev)')
                    writedlm(datadir("sol", "Sol-$(gridDim)D-rev-$(typeGridText)-generation-$(generationType)-reco-$(typeReco)-end-enableIons-false.dat"), solRev(endRev)')
                end

                ####################
                ## forward scan
                tGuess = tPrecond + 1.5 * tend

                powerDensity = - biasValues .* IV    # power density function
                MaxPD, indexPD = findmax(powerDensity)
                g2(t) = data.contactVoltageFunction[2](t) - biasValues[indexPD]
                tMPP = find_zero(g2, tGuess)

                g3(t) = data.contactVoltageFunction[2](t) - 1.0
                t1p0 = find_zero(g3, tGuess)

                if SRHvelocityETLn <= 200 * cm / s
                    g4(t) = data.contactVoltageFunction[2](t) - 1.2
                    t1p2 = find_zero(g4, tGuess)
                    ###########
                    g5(t) = data.contactVoltageFunction[2](t) - 1.25
                    t1p25 = find_zero(g5, tGuess)
                end

                inivalForw = tPrecond + tend
                endForw = tPrecond + 2 * tend

                if enableIons
                    if gridDim == 1
                        writedlm(datadir("sol", "Sol-$(gridDim)D-forw-$(typeGridText)-Ca-$textCa-generation-$(generationType)-reco-$(typeReco)-inival.dat"), solForw(inivalForw)')
                        writedlm(datadir("sol", "Sol-$(gridDim)D-forw-$(typeGridText)-Ca-$textCa-generation-$(generationType)-reco-$(typeReco)-MPP.dat"), solForw(tMPP)')
                        writedlm(datadir("sol", "Sol-$(gridDim)D-forw-$(typeGridText)-Ca-$textCa-generation-$(generationType)-reco-$(typeReco)-end.dat"), solForw(endForw)')
                    else
                        writedlm(datadir("sol", "Sol-$(gridDim)D-forw-$(typeGridText)-generation-$(generationType)-reco-$(typeReco)-inival.dat"), solForw(inivalForw)')
                        writedlm(datadir("sol", "Sol-$(gridDim)D-forw-$(typeGridText)-generation-$(generationType)-reco-$(typeReco)-MPP.dat"), solForw(tMPP)')
                        writedlm(datadir("sol", "Sol-$(gridDim)D-forw-$(typeGridText)-generation-$(generationType)-reco-$(typeReco)-end.dat"), solForw(endForw)')
                        ############
                        writedlm(datadir("sol", "Sol-$(gridDim)D-forw-$(typeGridText)-generation-$(generationType)-reco-$(typeReco)-V-1p0.dat"), solForw(t1p0)')
                        if SRHvelocityETLn <= 200 * cm / s
                            writedlm(datadir("sol", "Sol-$(gridDim)D-forw-$(typeGridText)-generation-$(generationType)-reco-$(typeReco)-V-1p2.dat"), solForw(t1p2)')
                            writedlm(datadir("sol", "Sol-$(gridDim)D-forw-$(typeGridText)-generation-$(generationType)-reco-$(typeReco)-V-1p25.dat"), solForw(t1p25)')
                        end
                    end
                else
                    writedlm(datadir("sol", "Sol-$(gridDim)D-forw-$(typeGridText)-generation-$(generationType)-reco-$(typeReco)-inival-enableIons-false.dat"), solForw(inivalForw)')
                    writedlm(datadir("sol", "Sol-$(gridDim)D-forw-$(typeGridText)-generation-$(generationType)-reco-$(typeReco)-MPP-enableIons-false.dat"), solForw(tMPP)')
                    writedlm(datadir("sol", "Sol-$(gridDim)D-forw-$(typeGridText)-generation-$(generationType)-reco-$(typeReco)-end-enableIons-false.dat"), solForw(endForw)')
                end

            end

            ## reverse I-V
            if enableIons
                if Ca == 6.0e22 # dont save when average vacancy dens is not equal to the original val
                    writedlm(datadir("IV", "IV-$(gridDim)D-rev-$(typeGridText)-generation-$(generationType)-reco-$(typeReco).dat"), [biasValuesRev IVRev])
                    ############
                    writedlm(datadir("IV", "JSRH-$(gridDim)D-rev-$(typeGridText)-generation-$(generationType)-reco-$(typeReco).dat"), [biasValuesRev IRevSRHn])
                    writedlm(datadir("IV", "JRad-$(gridDim)D-rev-$(typeGridText)-generation-$(generationType)-reco-$(typeReco).dat"), [biasValuesRev IRevRadn])
                    writedlm(datadir("IV", "JSRnL-$(gridDim)D-rev-$(typeGridText)-generation-$(generationType)-reco-$(typeReco).dat"), [biasValuesRev IRevSRnL])
                    writedlm(datadir("IV", "JSRpL-$(gridDim)D-rev-$(typeGridText)-generation-$(generationType)-reco-$(typeReco).dat"), [biasValuesRev IRevSRpL])
                    writedlm(datadir("IV", "JSRnR-$(gridDim)D-rev-$(typeGridText)-generation-$(generationType)-reco-$(typeReco).dat"), [biasValuesRev IRevSRnR])
                    writedlm(datadir("IV", "JSRpR-$(gridDim)D-rev-$(typeGridText)-generation-$(generationType)-reco-$(typeReco).dat"), [biasValuesRev IRevSRpR])
                end
            else
                writedlm(datadir("IV", "IV-$(gridDim)D-rev-$(typeGridText)-generation-$(generationType)-reco-$(typeReco)-enableIons-false.dat"), [biasValuesRev IVRev])
            end

            ## forward I-V
            if enableIons
                if Ca == 6.0e22 # dont save when average vacancy dens is not equal to the original val
                    writedlm(datadir("IV", "IV-$(gridDim)D-forw-$(typeGridText)-generation-$(generationType)-reco-$(typeReco).dat"), [biasValues IV])
                    ############
                    writedlm(datadir("IV", "JSRH-$(gridDim)D-forw-$(typeGridText)-generation-$(generationType)-reco-$(typeReco).dat"), [biasValues ISRHn])
                    writedlm(datadir("IV", "JRad-$(gridDim)D-forw-$(typeGridText)-generation-$(generationType)-reco-$(typeReco).dat"), [biasValues IRadn])
                    writedlm(datadir("IV", "JGen-$(gridDim)D-forw-$(typeGridText)-generation-$(generationType)-reco-$(typeReco).dat"), [biasValues IGen])
                    writedlm(datadir("IV", "JSRnL-$(gridDim)D-forw-$(typeGridText)-generation-$(generationType)-reco-$(typeReco).dat"), [biasValues ISRnL])
                    writedlm(datadir("IV", "JSRpL-$(gridDim)D-forw-$(typeGridText)-generation-$(generationType)-reco-$(typeReco).dat"), [biasValues ISRpL])
                    writedlm(datadir("IV", "JSRnR-$(gridDim)D-forw-$(typeGridText)-generation-$(generationType)-reco-$(typeReco).dat"), [biasValues ISRnR])
                    writedlm(datadir("IV", "JSRpR-$(gridDim)D-forw-$(typeGridText)-generation-$(generationType)-reco-$(typeReco).dat"), [biasValues ISRpR])
                end

            else
                writedlm(datadir("IV", "IV-$(gridDim)D-forw-$(typeGridText)-generation-$(generationType)-reco-$(typeReco)-enableIons-false.dat"), [biasValues IV])
                ############
                writedlm(datadir("IV", "JSRH-$(gridDim)D-forw-$(typeGridText)-generation-$(generationType)-reco-$(typeReco)-enableIons-false.dat"), [biasValues ISRHn])
                writedlm(datadir("IV", "JRad-$(gridDim)D-forw-$(typeGridText)-generation-$(generationType)-reco-$(typeReco)-enableIons-false.dat"), [biasValues IRadn])
                writedlm(datadir("IV", "JGen-$(gridDim)D-forw-$(typeGridText)-generation-$(generationType)-reco-$(typeReco)-enableIons-false.dat"), [biasValues IGen])
                writedlm(datadir("IV", "JSRnL-$(gridDim)D-forw-$(typeGridText)-generation-$(generationType)-reco-$(typeReco)-enableIons-false.dat"), [biasValues ISRnL])
                writedlm(datadir("IV", "JSRpL-$(gridDim)D-forw-$(typeGridText)-generation-$(generationType)-reco-$(typeReco)-enableIons-false.dat"), [biasValues ISRpL])
                writedlm(datadir("IV", "JSRnR-$(gridDim)D-forw-$(typeGridText)-generation-$(generationType)-reco-$(typeReco)-enableIons-false.dat"), [biasValues ISRnR])
                writedlm(datadir("IV", "JSRpR-$(gridDim)D-forw-$(typeGridText)-generation-$(generationType)-reco-$(typeReco)-enableIons-false.dat"), [biasValues ISRpR])
            end

        end

        if printText
            println("*** done\n")
        end

        integral = integrated_density(ctsys, sol = sol0p0, icc = iphia, ireg = regionPero)
        mOmega = data.regionVolumes[regionPero]

        println("Calculated average vacancy density is: ", integral / mOmega)
        println(" ")
        vacancyEnergy = data.params.bandEdgeEnergy[iphia, regionPero] / q
        println("Value for vacancy energy is: ", vacancyEnergy, " eV. Save this value for later use.")
        println(" ")


        testval = sum(filter(!isnan, sol0p0)) / length(sol0p0) # when using sparse storage, we get NaN values in solution
        return testval

    end #  main

    function test(; gridDim = 1, typeGrid = "planar", amplitude = 2.0e-7, demo_run = false)
        if demo_run
            if gridDim == 1
                testval = -0.6419535732692843 # all reco
            elseif gridDim == 2
                if typeGrid == "planar"
                    testval = -1.3136375582633124 # all reco
                elseif typeGrid == "nanotextured" && amplitude == 2.0e-7
                    testval = -1.3771899183248726 # all reco
                end
            end
        else
            if gridDim == 1
                testval = -0.6354871107155567 # all reco
            end
        end
        result = main(gridDim = gridDim, typeGrid = typeGrid, amplitude = amplitude, printText = false, demo_run = demo_run, generation = true, generationUniform = false, MaxwellSol = true)
        @info "result  = $result"
        @info "testval = $testval"
        return abs(result - testval) < 7.0e-3
    end

end # module
