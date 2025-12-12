# Parameters for a single-junction perovskite solar cell
# with the composition C60 -- triple cation -- PTAA
# The parameters follow: https://doi.org/10.25446/oxford.24359959
# from the publication Thiesbrummel et al., nature energy, 2024.

#####################################################################
############################ parameters ############################

@kwdef struct ParamsSingleJunction

    # unit factors
    cm = ufac"cm"
    nm = ufac"nm"
    K = ufac"K"
    m = ufac"m"
    s = ufac"s"
    V = ufac"V"

    q = ChargeTransport.constants.q
    k_B = ChargeTransport.constants.k_B

    eV = q * V

    ## set indices of the quasi Fermi potentials
    iphin = 1 # electron quasi Fermi potential
    iphip = 2 # hole quasi Fermi potential
    iphia = 3 # anion vacancy quasi Fermi potential
    ipsi = 4

    numberOfCarriers = 3 # electrons, holes and anion vacancies

    ########## device geometry ##########
    # region numbers
    regionETL = 1
    regionPero = 2
    regionHTL = 3

    regions = [regionETL, regionPero, regionHTL]
    numberOfRegions = length(regions)

    # boundary region numbers
    bregionLeft = 1
    bregionRight = 2
    bregionJ1 = 3
    bregionJ2 = 4
    bregionNoFlux = 5

    heightDev = 750.0 * nm

    ## length of regions
    h_ETL = 30.0 * nm # C60
    h_activePL = 400.0 * nm # perovskite
    h_HTL = 10.0 * nm # PTAA

    heightLayersPL = [h_ETL, h_ETL + h_activePL, h_ETL + h_activePL + h_HTL]
    h_totalPL = heightLayersPL[end]

    ########## physical values ##########

    ## charge numbers
    zn = -1
    zp = 1
    za = 1

    ## temperature
    T = 300.0 * K

    ## relative dielectric permittivity
    εr = [5.0, 22.0, 3.5] .* 1.0

    ΔEFLeft = 0.05 * eV # offset between metal and ETL
    ΔEFRight = 0.05 * eV # offset between metal and HTL

    ## band edge energies
    En = [-3.9, -3.9, -2.5] .* eV
    Ep = [-5.9, -5.53, -5.5] .* eV

    Ea1D = [0.0, -5.267, 0.0] .* eV
    ###################################################################
    ## for different amplitude lengths, save the energy for faster calculations (also post-processing more convenient with them)
    EaPlanar = [0.0, -5.28, 0.0] .* eV
    EaAmpl0p5e7 = [0.0, -5.286, 0.0] .* eV
    EaAmpl1p0e7 = [0.0, -5.303, 0.0] .* eV
    EaAmpl1p5e7 = [0.0, -5.327, 0.0] .* eV
    EaAmpl2p0e7 = [0.0, -5.355, 0.0] .* eV
    EaAmpl2p5e7 = [0.0, -5.384, 0.0] .* eV
    EaAmpl3p0e7 = [0.0, -5.412, 0.0] .* eV
    EaAmpl3p5e7 = [0.0, -5.436, 0.0] .* eV
    EaAmpl4p0e7 = [0.0, -5.459, 0.0] .* eV
    EaAmpl4p5e7 = [0.0, -5.478, 0.0] .* eV
    EaAmpl5p0e7 = [0.0, -5.495, 0.0] .* eV
    EaAmpl5p5e7 = [0.0, -5.51, 0.0] .* eV
    EaAmpl6p0e7 = [0.0, -5.523, 0.0] .* eV
    EaAmpl6p5e7 = [0.0, -5.535, 0.0] .* eV
    EaAmpl7p0e7 = [0.0, -5.547, 0.0] .* eV
    EaAmpl7p5e7 = [0.0, -5.558, 0.0] .* eV

    ## effective densities of density of states
    Nn1 = 1.0e26 / (m^3)
    Nn2 = 2.2e24 / (m^3)
    Nn3 = 1.0e26 / (m^3)

    Np1 = 1.0e26 / (m^3)
    Np2 = 2.2e24 / (m^3)
    Np3 = 1.0e26 / (m^3)

    Na2 = 1.0e27 / (m^3)

    Nn = [Nn1, Nn2, Nn3]
    Np = [Np1, Np2, Np3]
    Na = [0.0, Na2, 0.0]

    Da = 5.0e-14
    ## mobilities
    μn = [1.0e-6, 1.0e-4, 1.0e-8] .* (m^2) / (V * s)
    μp = [1.0e-6, 1.0e-4, 1.0e-8] .* (m^2) / (V * s)
    μa = [0.0, Da / (k_B * T / q), 0.0] .* (m^2) / (V * s) # 1.9e-12

    ## statistics functions
    Fcc = [Boltzmann, Boltzmann, FermiDiracMinusOne]

    ## radiative recombination
    r0 = [0.0, 3.0e-17, 0.0] .* m^3 / s

    ## life times
    τn = [1.0e100, 2.0e-7, 1.0e100] .* s
    τp = [1.0e100, 2.0e-7, 1.0e100] .* s

    ## trap densities
    ni1 = sqrt(Nn1 * Np1 * exp(- (En[1] - Ep[1]) / (k_B * T)))
    ni2 = sqrt(Nn2 * Np2 * exp(- (En[2] - Ep[2]) / (k_B * T)))
    ni3 = sqrt(Nn3 * Np3 * exp(- (En[3] - Ep[3]) / (k_B * T)))

    nτ = [ni1, ni2, ni3]
    pτ = [ni1, ni2, ni3]

    ## doping
    ## doping of transport layers are just for numerical stability, they are so small, they count as undoped.
    Cn = 1.0e20 / (m^3)
    Ca = 6.0e22 / (m^3)
    Cp = 1.0e20 / (m^3)

    ## generation
    incidentPhotonFlux = [0.0, 1.42e21, 0.0] ./ (m^2 * s)
    absorption = [0.0, 6.34e6, 0.0] ./ m
    generationPeak = h_ETL + h_activePL
    invertedIllumination = -1

    #####################################################################
    #####################      Newton Parameter      ####################

    damp_initial = 0.5
    damp_growth = 1.61 # >= 1
    max_round = 5
    maxiters = 1000

    abstol = 1.0e-7
    reltol = 1.0e-7
    tol_round = 1.0e-7
    Δu_opt = Inf

end
