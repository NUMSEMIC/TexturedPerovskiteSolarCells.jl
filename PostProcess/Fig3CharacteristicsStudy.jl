#=

Code for generating the solar cell characteristics parameter study

=#

module Fig3CharacteristicsStudy

using PyPlot
using DelimitedFiles
using ChargeTransport
using TexturedPerovskiteSolarCells

# for convenience
datadir = TexturedPerovskiteSolarCells.datadir
scriptsdir = TexturedPerovskiteSolarCells.scriptsdir

function main(;
        printText = true, saveFig = false,
        parameter_set = ParamsSingleJunction,
    )

    PyPlot.rc("font", family = "sans-serif", size = 16)
    PyPlot.rc("mathtext", fontset = "dejavusans")
    PyPlot.close("all")

    figsizeIV = (6.2, 5.7)
    figsize = (5.2, 6)
    fontsize = 25
    # use the destructuring operator to extract all the necessary parameters
    (; heightDev) = parameter_set()

    @local_unitfactors W m nm cm

    Textampl = [
        "5p0e-8", "1p0e-7", "1p5e-7", "2p0e-7", "2p5e-7", "3p0e-7", "3p5e-7", "4p0e-7",
        "4p5e-7", "5p0e-7", "5p5e-7", "6p0e-7", "6p5e-7", "7p0e-7", "7p5e-7",
    ]
    ##################################################################

    vETL = 2000
    vHTL = 500

    if printText
        println("Calculate for realistic cell (vETL = $(vETL), vHTL = $(vHTL))")
    end

    path = "vETL=$(vETL)vHTL=$(vHTL)/"

    ####################################################################################################
    ampl2 = collect(0.5:0.5:7.5) .* 1.0e-7; ampl = vcat(0.0, ampl2)

    textampl1 = "3p0e-7"
    textampl2 = "6p0e-7"

    Greys = get_cmap(:Greys)
    Greens = get_cmap(:Greens)

    figure(figsize = figsizeIV)
    IVPL = readdlm(datadir("IV", "$path/IV-2D-forw-planar-generation-Maxwell-reco-all.dat"))
    IVNT1 = readdlm(datadir("IV", "$path/IV-2D-forw-nanotextured-ampl-$textampl1-generation-Maxwell-reco-all.dat"))
    IVNT2 = readdlm(datadir("IV", "$path/IV-2D-forw-nanotextured-ampl-$textampl2-generation-Maxwell-reco-all.dat"))

    PyPlot.plot(IVPL[:, 1], -IVPL[:, 2] .* (cm^2) .* 1.0e3 ./ heightDev, linewidth = 5, color = Greys(251), label = "planar")
    PyPlot.plot(IVNT1[:, 1], -IVNT1[:, 2] .* (cm^2) .* 1.0e3 ./ heightDev, linewidth = 5, color = Greys(201), label = "textured (300 nm)")
    PyPlot.plot(IVNT2[:, 1], -IVNT2[:, 2] .* (cm^2) .* 1.0e3 ./ heightDev, linewidth = 5, color = Greys(131), label = "textured (600 nm)")

    PyPlot.xlim(0.1, 1.27)
    PyPlot.xticks([0.0, 0.4, 0.8, 1.2])
    PyPlot.ylim(-1, 23)
    yticks([0, 5, 10, 15, 20])
    PyPlot.legend(fontsize = 12)
    PyPlot.xlabel("Appl. voltage [V]", fontsize = fontsize)
    PyPlot.ylabel("\$ J \$ [mA cm\$^{-2} \$]", fontsize = fontsize)
    PyPlot.tick_params(which = "both", labelsize = fontsize)
    subplots_adjust(left = 0.17, right = 0.97, top = 0.97, bottom = 0.15)

    if saveFig
        savefig(datadir("IV-planar-nanotexture-generation-Maxwell-vETL=$(vETL)vHTL=$(vHTL).pdf"))
    end

    ####################################################################################################

    JSCVec1 = zeros(0); FFVec1 = zeros(0)
    VOCVec1 = zeros(0); EfficiencyVec1 = zeros(0)
    JMP1 = zeros(0); VMP1 = zeros(0)

    IVPL = readdlm(datadir("IV", "$path/IV-2D-forw-planar-generation-Maxwell-reco-all.dat"))

    IV = -IVPL[:, 2]
    bias = IVPL[:, 1]

    powerDensity = bias .* (IV)           # power density function
    MaxPD, indexPD = findmax(powerDensity)

    open_circuit = compute_open_circuit_voltage(bias, IV)
    IncLightPowerDens = 1000.0 * W / m^2
    fillfactor = 100 * (bias[indexPD] * IV[indexPD]) / (IV[1] * open_circuit)
    JSC = IV[1] ./ heightDev .* (0.01)^(2) .* 1.0e3
    efficiency = 100 * bias[indexPD] * (IV[indexPD] ./ heightDev) / (IncLightPowerDens)

    push!(JSCVec1, JSC)
    push!(FFVec1, fillfactor)
    push!(VOCVec1, open_circuit)
    push!(EfficiencyVec1, efficiency)
    push!(JMP1, IV[indexPD] ./ heightDev .* (0.01)^(2) .* 1.0e3)
    push!(VMP1, bias[indexPD])

    ##############

    for textampl in Textampl
        if printText
            println("Texture height = ", textampl, " m")
        end

        IVNT = readdlm(datadir("IV", "$path/IV-2D-forw-nanotextured-ampl-$textampl-generation-Maxwell-reco-all.dat"))

        IV = -IVNT[:, 2]
        bias = IVNT[:, 1]

        powerDensity = bias .* (IV)           # power density function
        MaxPD, indexPD = findmax(powerDensity)

        open_circuit = compute_open_circuit_voltage(bias, IV)
        fillfactor = 100 * (bias[indexPD] * IV[indexPD]) / (IV[1] * open_circuit)
        JSC = IV[1] ./ heightDev .* (0.01)^(2) .* 1.0e3
        efficiency = 100 * bias[indexPD] * (IV[indexPD] ./ heightDev) / (IncLightPowerDens)

        push!(JSCVec1, JSC)
        push!(FFVec1, fillfactor)
        push!(VOCVec1, open_circuit)
        push!(EfficiencyVec1, efficiency)
        push!(JMP1, IV[indexPD] ./ heightDev .* (0.01)^(2) .* 1.0e3)
        push!(VMP1, bias[indexPD])
    end

    ##################################################################
    vETL = 2000
    vHTL = 10

    if printText
        println("Calculate for realistic cell (vETL = $(vETL), vHTL = $(vHTL))")
    end

    path = "vETL=$(vETL)vHTL=$(vHTL)/"

    ####################################################################################################

    Oranges = get_cmap(:Oranges)

    figure(figsize = figsizeIV)
    IVPL = readdlm(datadir("IV", "$path/IV-2D-forw-planar-generation-Maxwell-reco-all.dat"))
    IVNT1 = readdlm(datadir("IV", "$path/IV-2D-forw-nanotextured-ampl-$textampl1-generation-Maxwell-reco-all.dat"))
    IVNT2 = readdlm(datadir("IV", "$path/IV-2D-forw-nanotextured-ampl-$textampl2-generation-Maxwell-reco-all.dat"))

    PyPlot.plot(IVPL[:, 1], -IVPL[:, 2] .* (cm^2) .* 1.0e3 ./ heightDev, linewidth = 5, color = Oranges(251), label = "planar")
    PyPlot.plot(IVNT1[:, 1], -IVNT1[:, 2] .* (cm^2) .* 1.0e3 ./ heightDev, linewidth = 5, color = Oranges(201), label = "textured (300 nm)")
    PyPlot.plot(IVNT2[:, 1], -IVNT2[:, 2] .* (cm^2) .* 1.0e3 ./ heightDev, linewidth = 5, color = Oranges(131), label = "textured (600 nm)")

    PyPlot.xlim(0.1, 1.27)
    PyPlot.xticks([0.0, 0.4, 0.8, 1.2])
    PyPlot.ylim(-1, 23)
    yticks([0, 5, 10, 15, 20])
    PyPlot.legend(fontsize = 12)
    PyPlot.xlabel("Appl. voltage [V]", fontsize = fontsize)
    PyPlot.ylabel("\$ J \$ [mA cm\$^{-2} \$]", fontsize = fontsize)
    PyPlot.tick_params(which = "both", labelsize = fontsize)
    subplots_adjust(left = 0.17, right = 0.97, top = 0.97, bottom = 0.15)

    if saveFig
        savefig(datadir("IV-planar-nanotexture-generation-Maxwell-vETL=$(vETL)vHTL=$(vHTL).pdf"))
    end

    ####################################################################################################

    JSCVec2 = zeros(0); FFVec2 = zeros(0)
    VOCVec2 = zeros(0); EfficiencyVec2 = zeros(0)
    JMP2 = zeros(0); VMP2 = zeros(0)

    IVPL = readdlm(datadir("IV", "$path/IV-2D-forw-planar-generation-Maxwell-reco-all.dat"))

    IV = -IVPL[:, 2]
    bias = IVPL[:, 1]

    powerDensity = bias .* (IV)           # power density function
    MaxPD, indexPD = findmax(powerDensity)

    open_circuit = compute_open_circuit_voltage(bias, IV)
    IncLightPowerDens = 1000.0 * W / m^2
    fillfactor = 100 * (bias[indexPD] * IV[indexPD]) / (IV[1] * open_circuit)
    JSC = IV[1] ./ heightDev .* (0.01)^(2) .* 1.0e3
    efficiency = 100 * bias[indexPD] * (IV[indexPD] ./ heightDev) / (IncLightPowerDens)

    push!(JSCVec2, JSC)
    push!(FFVec2, fillfactor)
    push!(VOCVec2, open_circuit)
    push!(EfficiencyVec2, efficiency)
    push!(JMP2, IV[indexPD] ./ heightDev .* (0.01)^(2) .* 1.0e3)
    push!(VMP2, bias[indexPD])

    ##############
    for textampl in Textampl
        if printText
            println("Texture height = ", textampl, " m")
        end

        IVNT = readdlm(datadir("IV", "$path/IV-2D-forw-nanotextured-ampl-$textampl-generation-Maxwell-reco-all.dat"))

        IV = -IVNT[:, 2]
        bias = IVNT[:, 1]

        powerDensity = bias .* (IV)           # power density function
        MaxPD, indexPD = findmax(powerDensity)

        open_circuit = compute_open_circuit_voltage(bias, IV)
        fillfactor = 100 * (bias[indexPD] * IV[indexPD]) / (IV[1] * open_circuit)
        JSC = IV[1] ./ heightDev .* (0.01)^(2) .* 1.0e3
        efficiency = 100 * bias[indexPD] * (IV[indexPD] ./ heightDev) / (IncLightPowerDens)

        push!(JSCVec2, JSC)
        push!(FFVec2, fillfactor)
        push!(VOCVec2, open_circuit)
        push!(EfficiencyVec2, efficiency)
        push!(JMP2, IV[indexPD] ./ heightDev .* (0.01)^(2) .* 1.0e3)
        push!(VMP2, bias[indexPD])
    end

    ##################################################################
    vETL = 10
    vHTL = 500

    if printText
        println("Calculate for realistic cell (vETL = $(vETL), vHTL = $(vHTL))")
    end

    path = "vETL=$(vETL)vHTL=$(vHTL)/"

    ####################################################################################################

    Greens = get_cmap(:Greens)

    figure(figsize = figsizeIV)
    IVPL = readdlm(datadir("IV", "$path/IV-2D-forw-planar-generation-Maxwell-reco-all.dat"))
    IVNT1 = readdlm(datadir("IV", "$path/IV-2D-forw-nanotextured-ampl-$textampl1-generation-Maxwell-reco-all.dat"))
    IVNT2 = readdlm(datadir("IV", "$path/IV-2D-forw-nanotextured-ampl-$textampl2-generation-Maxwell-reco-all.dat"))

    PyPlot.plot(IVPL[:, 1], -IVPL[:, 2] .* (cm^2) .* 1.0e3 ./ heightDev, linewidth = 5, color = Greens(251), label = "planar")
    PyPlot.plot(IVNT1[:, 1], -IVNT1[:, 2] .* (cm^2) .* 1.0e3 ./ heightDev, linewidth = 5, color = Greens(201), label = "textured (300 nm)")
    PyPlot.plot(IVNT2[:, 1], -IVNT2[:, 2] .* (cm^2) .* 1.0e3 ./ heightDev, linewidth = 5, color = Greens(131), label = "textured (600 nm)")

    PyPlot.xlim(0.1, 1.27)
    PyPlot.xticks([0.0, 0.4, 0.8, 1.2])
    PyPlot.ylim(-1, 23)
    yticks([0, 5, 10, 15, 20])
    PyPlot.legend(fontsize = 12)
    PyPlot.xlabel("Appl. voltage [V]", fontsize = fontsize)
    PyPlot.ylabel("\$ J \$ [mA cm\$^{-2} \$]", fontsize = fontsize)
    PyPlot.tick_params(which = "both", labelsize = fontsize)
    subplots_adjust(left = 0.17, right = 0.97, top = 0.97, bottom = 0.15)

    if saveFig
        savefig(datadir("IV-planar-nanotexture-generation-Maxwell-vETL=$(vETL)vHTL=$(vHTL).pdf"))
    end

    ####################################################################################################

    JSCVec3 = zeros(0); FFVec3 = zeros(0)
    VOCVec3 = zeros(0); EfficiencyVec3 = zeros(0)
    JMP3 = zeros(0); VMP3 = zeros(0)

    IVPL = readdlm(datadir("IV", "$path/IV-2D-forw-planar-generation-Maxwell-reco-all.dat"))

    IV = -IVPL[:, 2]
    bias = IVPL[:, 1]

    powerDensity = bias .* (IV)           # power density function
    MaxPD, indexPD = findmax(powerDensity)

    open_circuit = compute_open_circuit_voltage(bias, IV)
    IncLightPowerDens = 1000.0 * W / m^2
    fillfactor = 100 * (bias[indexPD] * IV[indexPD]) / (IV[1] * open_circuit)
    JSC = IV[1] ./ heightDev .* (0.01)^(2) .* 1.0e3
    efficiency = 100 * bias[indexPD] * (IV[indexPD] ./ heightDev) / (IncLightPowerDens)

    push!(JSCVec3, JSC)
    push!(FFVec3, fillfactor)
    push!(VOCVec3, open_circuit)
    push!(EfficiencyVec3, efficiency)
    push!(JMP3, IV[indexPD] ./ heightDev .* (0.01)^(2) .* 1.0e3)
    push!(VMP3, bias[indexPD])

    ##############
    for textampl in Textampl
        if printText
            println("Texture height = ", textampl, " m")
        end

        IVNT = readdlm(datadir("IV", "$path/IV-2D-forw-nanotextured-ampl-$textampl-generation-Maxwell-reco-all.dat"))

        IV = -IVNT[:, 2]
        bias = IVNT[:, 1]

        powerDensity = bias .* (IV)           # power density function
        MaxPD, indexPD = findmax(powerDensity)

        open_circuit = compute_open_circuit_voltage(bias, IV)
        fillfactor = 100 * (bias[indexPD] * IV[indexPD]) / (IV[1] * open_circuit)
        JSC = IV[1] ./ heightDev .* (0.01)^(2) .* 1.0e3
        efficiency = 100 * bias[indexPD] * (IV[indexPD] ./ heightDev) / (IncLightPowerDens)

        push!(JSCVec3, JSC)
        push!(FFVec3, fillfactor)
        push!(VOCVec3, open_circuit)
        push!(EfficiencyVec3, efficiency)
        push!(JMP3, IV[indexPD] ./ heightDev .* (0.01)^(2) .* 1.0e3)
        push!(VMP3, bias[indexPD])
    end

    ##################################################################
    vETL = 10
    vHTL = 10

    if printText
        println("Calculate for realistic cell (vETL = $(vETL), vHTL = $(vHTL))")
    end

    path = "vETL=$(vETL)vHTL=$(vHTL)/"

    ####################################################################################################

    Blues = get_cmap(:Blues)

    figure(figsize = figsizeIV)
    IVPL = readdlm(datadir("IV", "$path/IV-2D-forw-planar-generation-Maxwell-reco-all.dat"))
    IVNT1 = readdlm(datadir("IV", "$path/IV-2D-forw-nanotextured-ampl-$textampl1-generation-Maxwell-reco-all.dat"))
    IVNT2 = readdlm(datadir("IV", "$path/IV-2D-forw-nanotextured-ampl-$textampl2-generation-Maxwell-reco-all.dat"))

    PyPlot.plot(IVPL[:, 1], -IVPL[:, 2] .* (cm^2) .* 1.0e3 ./ heightDev, linewidth = 5, color = Blues(251), label = "planar")
    PyPlot.plot(IVNT1[:, 1], -IVNT1[:, 2] .* (cm^2) .* 1.0e3 ./ heightDev, linewidth = 5, color = Blues(201), label = "textured (300 nm)")
    PyPlot.plot(IVNT2[:, 1], -IVNT2[:, 2] .* (cm^2) .* 1.0e3 ./ heightDev, linewidth = 5, color = Blues(131), label = "textured (600 nm)")

    PyPlot.xlim(0.1, 1.27)
    PyPlot.xticks([0.0, 0.4, 0.8, 1.2])
    PyPlot.ylim(-1, 23)
    yticks([0, 5, 10, 15, 20])
    PyPlot.legend(fontsize = 12)
    PyPlot.xlabel("Appl. voltage [V]", fontsize = fontsize)
    PyPlot.ylabel("\$ J \$ [mA cm\$^{-2} \$]", fontsize = fontsize)
    PyPlot.tick_params(which = "both", labelsize = fontsize)
    subplots_adjust(left = 0.17, right = 0.97, top = 0.97, bottom = 0.15)

    if saveFig
        savefig(datadir("IV-planar-nanotexture-generation-Maxwell-vETL=$(vETL)vHTL=$(vHTL).pdf"))
    end

    ####################################################################################################

    JSCVec4 = zeros(0); FFVec4 = zeros(0)
    VOCVec4 = zeros(0); EfficiencyVec4 = zeros(0)
    JMP4 = zeros(0); VMP4 = zeros(0)

    IVPL = readdlm(datadir("IV", "$path/IV-2D-forw-planar-generation-Maxwell-reco-all.dat"))

    IV = -IVPL[:, 2]
    bias = IVPL[:, 1]

    powerDensity = bias .* (IV)           # power density function
    MaxPD, indexPD = findmax(powerDensity)

    open_circuit = compute_open_circuit_voltage(bias, IV)
    IncLightPowerDens = 1000.0 * W / m^2
    fillfactor = 100 * (bias[indexPD] * IV[indexPD]) / (IV[1] * open_circuit)
    JSC = IV[1] ./ heightDev .* (0.01)^(2) .* 1.0e3
    efficiency = 100 * bias[indexPD] * (IV[indexPD] ./ heightDev) / (IncLightPowerDens)

    push!(JSCVec4, JSC)
    push!(FFVec4, fillfactor)
    push!(VOCVec4, open_circuit)
    push!(EfficiencyVec4, efficiency)
    push!(JMP4, IV[indexPD] ./ heightDev .* (0.01)^(2) .* 1.0e3)
    push!(VMP4, bias[indexPD])

    ##############
    for textampl in Textampl
        if printText
            println("Texture height = ", textampl, " m")
        end

        IVNT = readdlm(datadir("IV", "$path/IV-2D-forw-nanotextured-ampl-$textampl-generation-Maxwell-reco-all.dat"))

        IV = -IVNT[:, 2]
        bias = IVNT[:, 1]

        powerDensity = bias .* (IV)           # power density function
        MaxPD, indexPD = findmax(powerDensity)

        open_circuit = compute_open_circuit_voltage(bias, IV)
        fillfactor = 100 * (bias[indexPD] * IV[indexPD]) / (IV[1] * open_circuit)
        JSC = IV[1] ./ heightDev .* (0.01)^(2) .* 1.0e3
        efficiency = 100 * bias[indexPD] * (IV[indexPD] ./ heightDev) / (IncLightPowerDens)

        push!(JSCVec4, JSC)
        push!(FFVec4, fillfactor)
        push!(VOCVec4, open_circuit)
        push!(EfficiencyVec4, efficiency)
        push!(JMP4, IV[indexPD] ./ heightDev .* (0.01)^(2) .* 1.0e3)
        push!(VMP4, bias[indexPD])
    end

    ##################################################################
    ##################################################################
    ## efficiency plot
    ##################################################################

    ampl2 = collect(0.5:0.5:7.5) .* 1.0e-7; ampl = vcat(0.0, ampl2)
    MKsize = 12
    ColBlue = [168 / 255, 168 / 255, 168 / 255]
    ColGold = [172 / 255, 235 / 255, 180 / 255]

    figure(figsize = figsize)
    PyPlot.plot(ampl ./ nm, EfficiencyVec1, "o", markersize = MKsize, color = ColBlue, label = "ETL = 2000, HTL = 500")
    PyPlot.plot(ampl ./ nm, EfficiencyVec2, "o", markersize = MKsize, color = Oranges(121), label = "HTL = 10")
    PyPlot.plot(ampl ./ nm, EfficiencyVec3, "o", markersize = MKsize, color = ColGold, label = "ETL = 10")
    PyPlot.plot(ampl ./ nm, EfficiencyVec4, "o", markersize = MKsize, color = "dodgerblue", label = "HTL = ETL = 10")

    PyPlot.xticks([0, 300, 600])
    PyPlot.xlim(-20.0, 770)
    PyPlot.ylim(12.5, 23.5)
    PyPlot.yticks([15, 17.5, 20, 22.5])
    PyPlot.xlabel("Texture height [nm]", fontsize = fontsize)
    PyPlot.ylabel("PCE [%]", fontsize = fontsize)
    PyPlot.tick_params(which = "both", labelsize = fontsize)
    PyPlot.legend(fontsize = 12)
    subplots_adjust(left = 0.27, right = 0.97, top = 0.97, bottom = 0.15)

    Max1 = maximum(EfficiencyVec1); ArgMax1 = argmax(EfficiencyVec1)
    Max2 = maximum(EfficiencyVec2); ArgMax2 = argmax(EfficiencyVec2)
    Max3 = maximum(EfficiencyVec3); ArgMax3 = argmax(EfficiencyVec3)
    Max4 = maximum(EfficiencyVec4); ArgMax4 = argmax(EfficiencyVec4)

    if printText
        println(".................................................")
        println(" ")
        println("Planar (standard):       PCE is $(round(EfficiencyVec1[1], digits = 3)) %.")
        println("Planar (HTL = 10):       PCE is $(round(EfficiencyVec2[1], digits = 3)) %.")
        println("Planar (ETL = 10):       PCE is $(round(EfficiencyVec3[1], digits = 3)) %.")
        println("Planar (HTL = ETL = 10): PCE is $(round(EfficiencyVec4[1], digits = 3)) %.")
        println(" ")
        println("Realistic (standard):       Maximum PCE is $(round(Max1, digits = 3)) %, reached at $(round(ampl[ArgMax1] ./ nm, digits = 1)) nm.")
        println("Realistic (HTL = 10):       Maximum PCE is $(round(Max2, digits = 3)) %, reached at $(round(ampl[ArgMax2] ./ nm, digits = 1)) nm.")
        println("Realistic (ETL = 10):       Maximum PCE is $(round(Max3, digits = 3)) %, reached at $(round(ampl[ArgMax3] ./ nm, digits = 1)) nm.")
        println("Realistic (HTL = ETL = 10): Maximum PCE is $(round(Max4, digits = 3)) %, reached at $(round(ampl[ArgMax4] ./ nm, digits = 1)) nm.")
    end

    if saveFig
        savefig(datadir("ampl-efficiency-generation-Maxwell.pdf"))
    end

    #################################################################
    #################################################################

    JSCOptical = readdlm(datadir("Max-PhotoCurrent-Optical.dat"))

    figure(figsize = figsize)
    PyPlot.plot(ampl ./ nm, JSCOptical, "o", markersize = MKsize, color = [255 / 255, 108 / 255, 108 / 255], label = "optical maximum")
    PyPlot.plot(ampl ./ nm, JSCVec1, "o", markersize = MKsize, color = ColBlue, label = "ETL = 2000, HTL = 500")
    PyPlot.plot(ampl ./ nm, JSCVec2, "o", markersize = 1.4 * MKsize, color = Oranges(121), label = "HTL = 10")
    PyPlot.plot(ampl ./ nm, JSCVec3, "o", markersize = MKsize, color = ColGold, label = "ETL = 10")
    PyPlot.plot(ampl ./ nm, JSCVec4, "o", markersize = MKsize, color = "dodgerblue", label = "HTL = ETL = 10")

    PyPlot.xticks([0, 300, 600])
    PyPlot.xlim(-20.0, 770)
    #PyPlot.ylim(21, 24)
    PyPlot.xlabel("Texture height [nm]", fontsize = fontsize)
    PyPlot.ylabel("\$ J_{\\mathrm{SC}}\$ [mA cm\$^{-2} \$]", fontsize = fontsize)
    PyPlot.tick_params(which = "both", labelsize = fontsize)
    PyPlot.legend(fontsize = 12)
    subplots_adjust(left = 0.27, right = 0.97, top = 0.97, bottom = 0.15)

    if saveFig
        savefig(datadir("ampl-JSC-generation-Maxwell.pdf"))
    end

    if printText
        println(".................................................")
        println(" ")
        println("Planar (standard):       JSC is $(round(JSCVec1[1], digits = 3)) mA/cm^2.")
        println("Planar (HTL = 10):       JSC is $(round(JSCVec2[1], digits = 3)) mA/cm^2.")
        println("Planar (ETL = 10):       JSC is $(round(JSCVec3[1], digits = 3)) mA/cm^2.")
        println("Planar (HTL = ETL = 10): JSC is $(round(JSCVec4[1], digits = 3)) mA/cm^2.")
        println(" ")
        println("Realistic (standard):       Maximum JSC is $(round(JSCVec1[ArgMax1], digits = 3)) mA/cm^2.")
        println("Realistic (HTL = 10):       Maximum JSC is $(round(JSCVec2[ArgMax2], digits = 3)) mA/cm^2.")
        println("Realistic (ETL = 10):       Maximum JSC is $(round(JSCVec3[ArgMax3], digits = 3)) mA/cm^2.")
        println("Realistic (HTL = ETL = 10): Maximum JSC is $(round(JSCVec4[ArgMax4], digits = 3)) mA/cm^2.")
    end

    ##################################################################
    figure(figsize = figsize)
    PyPlot.plot(ampl ./ nm, FFVec1, "o", markersize = MKsize, color = ColBlue, label = "ETL = 2000, HTL = 500")
    PyPlot.plot(ampl ./ nm, FFVec2, "o", markersize = MKsize, color = Oranges(121), label = "HTL = 10")
    PyPlot.plot(ampl ./ nm, FFVec3, "o", markersize = MKsize, color = ColGold, label = "ETL = 10")
    PyPlot.plot(ampl ./ nm, FFVec4, "o", markersize = MKsize, color = "dodgerblue", label = "HTL = ETL = 10")

    PyPlot.xticks([0, 300, 600])
    PyPlot.xlim(-20.0, 770)
    PyPlot.ylim(65, 85)
    PyPlot.yticks([70, 80])
    PyPlot.xlabel("Texture height [nm]", fontsize = fontsize)
    PyPlot.ylabel("FF [%]", fontsize = fontsize)
    PyPlot.tick_params(which = "both", labelsize = fontsize)
    PyPlot.legend(fontsize = 12)
    subplots_adjust(left = 0.27, right = 0.97, top = 0.97, bottom = 0.15)

    if saveFig
        savefig(datadir("ampl-FF-generation-Maxwell.pdf"))
    end

    if printText
        println(".................................................")
        println(" ")
        println("Planar (standard):       FF is $(round(FFVec1[1], digits = 3))  %.")
        println("Planar (HTL = 10):       FF is $(round(FFVec2[1], digits = 3))  %.")
        println("Planar (ETL = 10):       FF is $(round(FFVec3[1], digits = 3)) %.")
        println("Planar (HTL = ETL = 10): FF is $(round(FFVec4[1], digits = 3))  %.")
        println(" ")
        println("Realistic (standard):       Maximum FF is $(round(FFVec1[ArgMax1], digits = 3))  %.")
        println("Realistic (HTL = 10):       Maximum FF is $(round(FFVec2[ArgMax2], digits = 3))  %.")
        println("Realistic (ETL = 10):       Maximum FF is $(round(FFVec3[ArgMax3], digits = 3))  %.")
        println("Realistic (HTL = ETL = 10): Maximum FF is $(round(FFVec4[ArgMax4], digits = 3)) %.")
    end

    #################################

    figure(figsize = figsize)
    PyPlot.plot(ampl ./ nm, VOCVec1, "o", markersize = MKsize, color = ColBlue, label = "ETL = 2000, HTL = 500")
    PyPlot.plot(ampl ./ nm, VOCVec2, "o", markersize = MKsize, color = Oranges(121), label = "HTL = 10")
    PyPlot.plot(ampl ./ nm, VOCVec3, "o", markersize = MKsize, color = ColGold, label = "ETL = 10")
    PyPlot.plot(ampl ./ nm, VOCVec4, "o", markersize = MKsize, color = "dodgerblue", label = "HTL = ETL = 10")

    PyPlot.xticks([0, 300, 600])
    PyPlot.xlim(-20.0, 770)
    PyPlot.ylim(1.0, 1.27)
    PyPlot.yticks([1.0, 1.1, 1.2])
    PyPlot.xlabel("Texture height [nm]", fontsize = fontsize)
    PyPlot.ylabel("\$ V_{\\mathrm{OC}} \$ [V]", fontsize = fontsize)
    PyPlot.tick_params(which = "both", labelsize = fontsize)
    PyPlot.legend(fontsize = 12)
    subplots_adjust(left = 0.27, right = 0.97, top = 0.97, bottom = 0.15)

    if saveFig
        savefig(datadir("ampl-VOC-generation-Maxwell.pdf"))
    end

    if printText
        println(".................................................")
        println(" ")
        println("Planar (standard):       VOC is $(round(VOCVec1[1], digits = 3))  V.")
        println("Planar (HTL = 10):       VOC is $(round(VOCVec2[1], digits = 3))  V.")
        println("Planar (ETL = 10):       VOC is $(round(VOCVec3[1], digits = 3)) V.")
        println("Planar (HTL = ETL = 10): VOC is $(round(VOCVec4[1], digits = 3))  V.")
        println(" ")
        println("Realistic (standard):       Maximum VOC is $(round(VOCVec1[ArgMax1], digits = 3))  V.")
        println("Realistic (HTL = 10):       Maximum VOC is $(round(VOCVec2[ArgMax2], digits = 3))  V.")
        println("Realistic (ETL = 10):       Maximum VOC is $(round(VOCVec3[ArgMax3], digits = 3))  V.")
        println("Realistic (HTL = ETL = 10): Maximum VOC is $(round(VOCVec4[ArgMax4], digits = 3)) V.")
    end

    #################################
    VOCDiff1 = zeros(0)
    VOCDiff2 = zeros(0)
    VOCDiff3 = zeros(0)
    VOCDiff4 = zeros(0)

    for ii in 2:length(VOCVec1)

        VOCdiff1 = 1.0e3 * (VOCVec1[ii] - VOCVec1[1])
        VOCdiff2 = 1.0e3 * (VOCVec2[ii] - VOCVec2[1])
        VOCdiff3 = 1.0e3 * (VOCVec3[ii] - VOCVec3[1])
        VOCdiff4 = 1.0e3 * (VOCVec4[ii] - VOCVec4[1])

        push!(VOCDiff1, VOCdiff1)
        push!(VOCDiff2, VOCdiff2)
        push!(VOCDiff3, VOCdiff3)
        push!(VOCDiff4, VOCdiff4)

    end

    figure(figsize = figsize)
    PyPlot.plot(ampl[2:end] ./ nm, VOCDiff1, "o", markersize = MKsize, color = ColBlue, label = "ETL = 2000, HTL = 500")
    PyPlot.plot(ampl[2:end] ./ nm, VOCDiff2, "o", markersize = MKsize, color = Oranges(121), label = "HTL = 10")
    PyPlot.plot(ampl[2:end] ./ nm, VOCDiff3, "o", markersize = MKsize, color = ColGold, label = "ETL = 10")
    PyPlot.plot(ampl[2:end] ./ nm, VOCDiff4, "o", markersize = MKsize, color = "dodgerblue", label = "HTL = ETL = 10")

    PyPlot.xticks([0, 300, 600])
    PyPlot.xlim(-20.0, 770)
    PyPlot.ylim(-95, 25)
    PyPlot.yticks([-80, -40, -20, 0, 20])
    PyPlot.xlabel("Texture height [nm]", fontsize = fontsize)
    PyPlot.ylabel("\$ V_{\\mathrm{OC, NT}} - V_{\\mathrm{OC, PL}} \$ [mV]", fontsize = fontsize)
    PyPlot.tick_params(which = "both", labelsize = fontsize)
    PyPlot.legend(fontsize = 12)
    subplots_adjust(left = 0.27, right = 0.97, top = 0.97, bottom = 0.15)

    if saveFig
        savefig(datadir("ampl-VOC-difference-generation-Maxwell.pdf"))
    end

    return nothing
end

end
