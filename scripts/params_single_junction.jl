

#####################################################################
############################ parameters ############################

paramsname           = "2"

## set indices of the quasi Fermi potentials
iphin                = 1 # electron quasi Fermi potential
iphip                = 2 # hole quasi Fermi potential
iphia                = 3 # anion vacancy quasi Fermi potential
ipsi                 = 4

numberOfCarriers     = 3 # electrons, holes and anion vacancies

########## device geometry ##########
# region numbers
regionETL1           = 1
regionPero           = 2
regionHTL            = 3

regions              = [regionETL1, regionPero, regionHTL]
numberOfRegions      = length(regions)

# boundary region numbers
bregionLeft          = 1
bregionRight         = 2
bregionJ1            = 3
bregionJ2            = 4
bregionNoFlux        = 5

heightDev            = 750.0 * nm

## length of regions
h_ETL1               = 30.0  * nm # C60
h_activePL           = 400.0 * nm # perovskite
h_HTL                = 10.0  * nm # PTAA

heightLayersPL       = [h_ETL1, h_ETL1 + h_activePL, h_ETL1 +  h_activePL + h_HTL]
h_totalPL            = heightLayersPL[end]

########## physical values ##########

## charge numbers
zn                   = -1
zp                   = 1
za                   = 1

## temperature
T                    = 300.0                                     *  K
UT                   = kB * T / q

## relative dielectric permittivity
εr                   = [5.0, 22.0, 3.5]                         .* 1.0

## band edge energies
En                   = [-3.9, -3.9,  -2.5]                      .*  eV
Ep                   = [-5.9, -5.53, -5.5]                      .*  eV

Ea1D                 = [0.0,  -5.361, 0.0]                      .*  eV

###################################################################
## for different amplitude lengths
# EaLoop = -5.365 # int Ea cal: 5.997578454738973e22 # (int - Ca) / Ca = -0.0004035908768378476
EaPlanar             = [0.0,  -5.365 , 0.0]                     .*  eV
# EaLoop = -5.368 # int Ea cal: 5.9982616585419165e22 # (int - Ca) / Ca = -0.0002897235763472455
EaAmpl0p5e7          = [0.0,  -5.368 , 0.0]                     .*  eV
# EaLoop = -5.378 # int Ea cal: 5.997072291536053e22 # (int - Ca) / Ca = -0.0004879514106577901
EaAmpl1p0e7          = [0.0,  -5.378 , 0.0]                     .*  eV
# EaLoop = -5.392 # int Ea cal: 6.003544241339036e22 # (int - Ca) / Ca = 0.0005907068898393456
EaAmpl1p5e7          = [0.0,  -5.392 , 0.0]                     .*  eV
# EaLoop = -5.41 # int Ea cal: 5.998555990281471e22 # (int - Ca) / Ca = -0.00024066828642154797
EaAmpl2p0e7          = [0.0,  -5.410 , 0.0]                     .*  eV
# EaLoop = -5.425 # int Ea cal: 6.004910764089687e22 # (int - Ca) / Ca = 0.0008184606816145332
EaAmpl2p5e7          = [0.0,  -5.425,   0.0]                    .*  eV
# EaLoop = -5.441 # int Ea cal: 6.006807857887035e22 # (int - Ca) / Ca = 0.0011346429811724648
EaAmpl3p0e7          = [0.0,  -5.441,   0.0]                    .*  eV
# EaLoop = -5.457 # int Ea cal: 5.998931157924821e22 # (int - Ca) / Ca = -0.00017814034586320678
EaAmpl3p5e7          = [0.0,  -5.457,   0.0]                    .*  eV
# EaLoop = -5.471 # int Ea cal: 5.996443836847383e22 # (int - Ca) / Ca = -0.0005926938587694455
EaAmpl4p0e7          = [0.0,  -5.471,   0.0]                    .*  eV
# EaLoop = -5.482 # int Ea cal: 6.0043888554251045e22 # (int - Ca) / Ca = 0.000731475904184084
EaAmpl4p5e7          = [0.0,  -5.482 , 0.0]                     .*  eV
# EaLoop = -5.493 # int Ea cal: 6.0002047196495195e22 # (int - Ca) / Ca = 3.4119941586577e-5
EaAmpl5p0e7          = [0.0,  -5.493 , 0.0]                     .*  eV
# EaLoop = -5.503 # int Ea cal: 6.006123717439521e22 # (int - Ca) / Ca = 0.0010206195732534614
EaAmpl5p5e7          = [0.0,  -5.503 , 0.0]                     .*  eV
# EaLoop = -5.513 # int Ea cal: 6.009538163397028e22 # (int - Ca) / Ca = 0.0015896938995047028
EaAmpl6p0e7          = [0.0,  -5.513 , 0.0]                     .*  eV
# EaLoop = -5.523 # int Ea cal: 6.004118957467857e22 # (int - Ca) / Ca = 0.0006864929113094751
EaAmpl6p5e7          = [0.0,  -5.523 , 0.0]                     .*  eV
# EaLoop = -5.533 # int Ea cal: 5.989365556878181e22 # (int - Ca) / Ca = -0.001772407186969784
EaAmpl7p0e7          = [0.0,  -5.533 , 0.0]                     .*  eV
# EaLoop = -5.541 # int Ea cal: 6.010879247846272e22 # (int - Ca) / Ca = 0.0018132079743787336
EaAmpl7p5e7          = [0.0,  -5.541 , 0.0]                     .*  eV
# EaLoop = -5.55 # int Ea cal: 5.988885681337782e22 # (int - Ca) / Ca = -0.0018523864437030037
EaAmpl8p0e7          = [0.0,  -5.55, 0.0]                     .*  eV

## effective densities of density of states
Nn1                  = 1.0e26                                   / (m^3)
Nn2                  = 2.2e24                                   / (m^3)
Nn3                  = 1.0e26                                   / (m^3)

Np1                  = 1.0e26                                   / (m^3)
Np2                  = 2.2e24                                   / (m^3)
Np3                  = 1.0e26                                   / (m^3)

Na2                  = 1.0e27                                   / (m^3)

Nn                   = [Nn1, Nn2, Nn3]
Np                   = [Np1, Np2, Np3]
Na                   = [0.0, Na2, 0.0]

## mobilities
μn                   = [1.0e-6, 5.0e-4,  1.0e-8]               .* (m^2) / (V * s)
μp                   = [1.0e-6, 5.0e-4,  1.0e-8]               .* (m^2) / (V * s)
μa                   = [0.0,    1.9e-12, 0.0]                  .* (m^2) / (V * s)

## statistics functions
Fcc                  = [Boltzmann, Boltzmann, FermiDiracMinusOne]

## radiative recombination
r0                   = [0.0, 3.0e-17, 0.0]                     .* m^3 / s

## recombination velocities
SRHvelocityETLn      = 1.0e5                                    * m / s
SRHvelocityETLp      = 20.0                                     * m / s

SRHvelocityHTLn      = 1.0                                      * m / s
SRHvelocityHTLp      = 1.0e5                                    * m / s

## life times
τn                   = [1.0e100, 4.0e-8, 1.0e100]              .* s
τp                   = [1.0e100, 4.0e-8, 1.0e100]              .* s

## trap densities
ni1                  = sqrt( Nn1 * Np1 * exp( - (En[1] - Ep[1])/(kB * T)) )
ni2                  = sqrt( Nn2 * Np2 * exp( - (En[2] - Ep[2])/(kB * T)) )
ni3                  = sqrt( Nn3 * Np3 * exp( - (En[3] - Ep[3])/(kB * T)) )

nτ                   = [ni1, ni2, ni3]
pτ                   = [ni1, ni2, ni3]

## doping
Cn1                  = 2.09e24                                  / (m^3)
# DA: Ea values available for 1.0e21, 1.0e22, 1.0e23, 1.0e24
Ca                   = 6.0e22                                   / (m^3)
Cp                   = 2.09e24                                  / (m^3)

## generation
incidentPhotonFlux   = [0.0, 1.4e21, 0.0]                      ./ (m^2 * s)
absorption           = [0.0, 1.3e7,  0.0]                      ./ m
generationPeak       = h_ETL1 + h_activePL
invertedIllumination = -1

#####################################################################
#####################      Newton Parameter      ####################

damp_initial = 0.5
damp_growth  = 1.61 # >= 1
max_round    = 5
maxiters     = 1000

abstol       = 1.0e-7
reltol       = 1.0e-7
tol_round    = 1.0e-7
Δu_opt       = Inf
