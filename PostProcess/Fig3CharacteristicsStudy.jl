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

function main(;printText = true, saveFig = false,
            scanrate =  "1000p0", # "0p001", #  "10p0", # "1000p0", #
            generation = "Maxwell", # "uniform"
            parameter_file = scriptsdir("params_single_junction.jl"),
            )

    include(parameter_file)

    PyPlot.rc("font", family="sans-serif", size=14)
    PyPlot.rc("mathtext", fontset="dejavusans")
    PyPlot.close("all")

    path = "parameter-$paramsname/$scanrate"

    JSCVec = zeros(0); FFVec         = zeros(0)
    VOCVec = zeros(0); EfficiencyVec = zeros(0)
    JMP    = zeros(0); VMP           = zeros(0)

    @info "read " * datadir("IV", "$path/IV-2D-forw-planar-generation-$generation-reco-all.dat")
    IVPL = readdlm(datadir("IV", "$path/IV-2D-forw-planar-generation-$generation-reco-all.dat"), Float64)

    IV                = -IVPL[:, 2]
    bias              =  IVPL[:, 1]

    powerDensity      = bias .* (IV)           # power density function
    MaxPD, indexPD    = findmax(powerDensity)

    open_circuit      = compute_open_circuit_voltage(bias, IV)

    IncLightPowerDens = 1000.0 * W/m^2

    fillfactor        = 100 * (bias[indexPD] * IV[indexPD]) / (IV[1] * open_circuit)

    JSC               = IV[1] ./heightDev.*(0.01)^(2).*1.0e3

    efficiency        = 100 * bias[indexPD] * (IV[indexPD]./heightDev) / (IncLightPowerDens)

    push!(JSCVec, JSC);          push!(FFVec, fillfactor)
    push!(VOCVec, open_circuit); push!(EfficiencyVec, efficiency)
    push!(JMP, IV[indexPD]./heightDev.*(0.01)^(2).*1.0e3); push!(VMP, bias[indexPD])

    #############################################################################################

    Textampl = ["5p0e-8", "1p0e-7", "1p5e-7", "2p0e-7", "2p5e-7", "3p0e-7", "3p5e-7", "4p0e-7",
                "4p5e-7", "5p0e-7", "5p5e-7", "6p0e-7", "6p5e-7", "7p0e-7", "7p5e-7"]

    if printText
        println("Calculate for realistic cell")
    end
    for textampl in Textampl
        if printText
            println("Texture height = ", textampl, " m")
        end

        @info "read " * datadir("IV", "$path/IV-2D-forw-nanotextured-ampl-$textampl-generation-$generation-reco-all.dat")
        IVNT = readdlm(datadir("IV", "$path/IV-2D-forw-nanotextured-ampl-$textampl-generation-$generation-reco-all.dat"), Float64)

        IV                = -IVNT[:, 2]
        bias              =  IVNT[:, 1]

        powerDensity      = bias .* (IV)           # power density function
        MaxPD, indexPD    = findmax(powerDensity)

        open_circuit      = compute_open_circuit_voltage(bias, IV)

        IncLightPowerDens = 1000.0 * W/m^2

        fillfactor        = 100 * (bias[indexPD] * IV[indexPD]) / (IV[1] * open_circuit)

        JSC               = IV[1] ./heightDev.*(0.01)^(2).*1.0e3

        efficiency        = 100 * bias[indexPD] * (IV[indexPD]./heightDev) / (IncLightPowerDens)

        push!(JSCVec, JSC);          push!(FFVec, fillfactor)
        push!(VOCVec, open_circuit); push!(EfficiencyVec, efficiency)
        push!(JMP, IV[indexPD]./heightDev.*(0.01)^(2).*1.0e3);     push!(VMP, bias[indexPD])
    end

    ################################################################
    JSCVecRad = zeros(0); FFVecRad         = zeros(0)
    VOCVecRad = zeros(0); EfficiencyVecRad = zeros(0)
    JMPRad    = zeros(0); VMPRad           = zeros(0)

    @info "read " * datadir("IV", "$path/IV-2D-forw-planar-generation-$generation-reco-radiative.dat")
    IVPL = readdlm(datadir("IV", "$path/IV-2D-forw-planar-generation-$generation-reco-radiative.dat"), Float64)

    IV                = -IVPL[:, 2]
    bias              =  IVPL[:, 1]

    powerDensity      = bias .* (IV)           # power density function
    MaxPD, indexPD    = findmax(powerDensity)

    open_circuit      = compute_open_circuit_voltage(bias, IV)

    IncLightPowerDens = 1000.0 * W/m^2

    fillfactor        = 100 * (bias[indexPD] * IV[indexPD]) / (IV[1] * open_circuit)

    JSC               = IV[1] ./heightDev.*(0.01)^(2).*1.0e3

    efficiency        = 100 * bias[indexPD] * (IV[indexPD]./heightDev) / (IncLightPowerDens)

    push!(JSCVecRad, JSC);          push!(FFVecRad, fillfactor)
    push!(VOCVecRad, open_circuit); push!(EfficiencyVecRad, efficiency)
    push!(JMPRad, IV[indexPD]./heightDev.*(0.01)^(2).*1.0e3); push!(VMPRad, bias[indexPD])


    if printText
        println(" ")
        println("Calculate for ideal cell")
    end

    for textampl in Textampl
        if printText
            println("Texture height = ", textampl, " m")
        end

        IVNT = readdlm(datadir("IV", "$path/IV-2D-forw-nanotextured-ampl-$textampl-generation-$generation-reco-radiative.dat"))

        IV                = -IVNT[:, 2]
        bias              =  IVNT[:, 1]

        powerDensity      = bias .* (IV)           # power density function
        MaxPD, indexPD    = findmax(powerDensity)

        open_circuit      = compute_open_circuit_voltage(bias, IV)

        IncLightPowerDens = 1000.0 * W/m^2

        fillfactor        = 100 * (bias[indexPD] * IV[indexPD]) / (IV[1] * open_circuit)

        JSC               = IV[1] ./heightDev.*(0.01)^(2).*1.0e3

        efficiency        = 100 * bias[indexPD] * (IV[indexPD]./heightDev) / (IncLightPowerDens)

        push!(JSCVecRad, JSC);          push!(FFVecRad, fillfactor)
        push!(VOCVecRad, open_circuit); push!(EfficiencyVecRad, efficiency)
        push!(JMPRad, IV[indexPD]./heightDev.*(0.01)^(2).*1.0e3);  push!(VMPRad, bias[indexPD])
    end

    ampl2    = collect(0.5:0.5:7.5) .* 1.0e-7; ampl = vcat(0.0, ampl2)
    MKsize  = 12
    ColBlue = [168/255, 168/255, 168/255]
    ColGold = [172/255, 235/255, 180/255]

    ##################################################################
    ##################################################################
    ## efficiency plot
    ##################################################################

    PyPlot.plot(ampl./nm, EfficiencyVecRad, "o", markersize = MKsize, color = ColBlue, label = "ideal")
    PyPlot.plot(ampl./nm, EfficiencyVec,    "o", markersize = MKsize, color = ColGold, label = "realistic")

    PyPlot.xticks([0, 300, 600])
    PyPlot.xlim(-20.0, 770)
    PyPlot.ylim(20, 29.5)
    PyPlot.xlabel("Texture height [nm]", fontsize=18)
    PyPlot.ylabel("PCE [%]", fontsize=18)
    PyPlot.tick_params(which ="both", labelsize=18)
    PyPlot.legend(fontsize = 18)
    PyPlot.tight_layout()

    Max1Rad = maximum(EfficiencyVecRad); Max2Rad = argmax(EfficiencyVecRad)
    Max1    = maximum(EfficiencyVec);    Max2    = argmax(EfficiencyVec)

    if printText
        println(".................................................")
        println(" ")
        println("Planar (realistic): PCE is $(EfficiencyVec[1]) %.")
        println("Planar (ideal):     PCE is $(EfficiencyVecRad[1]) %.")

        println(" ")

        println("Realistic: Maximum PCE is $(Max1) %, reached at $(ampl[Max2]./nm) nm.")
        println("Ideal:     Maximum PCE is $(Max1Rad) %, reached at $(ampl[Max2Rad]./nm) nm.")
    end

    if saveFig
        savefig(datadir("ampl-efficiency-params-$paramsname-generation-$generation-scanrate-$scanrate.pdf"))
    end
    #####################

    EffDiff = zeros(0); EffDiffRad = zeros(0)

    for ii = 2:length(EfficiencyVec)
        Effdiff    = EfficiencyVec[ii]    - EfficiencyVec[1]
        EffdiffRad = EfficiencyVecRad[ii] - EfficiencyVecRad[1]

        push!(EffDiff, Effdiff); push!(EffDiffRad, EffdiffRad)
    end

    figure()
    PyPlot.plot(ampl[2:end]./nm, EffDiffRad, "o", markersize = MKsize, color = ColBlue, label = "ideal")
    PyPlot.plot(ampl[2:end]./nm, EffDiff,    "o", markersize = MKsize, color = ColGold, label = "realistic")

    PyPlot.xticks([0, 300, 600])
    PyPlot.xlim(-20.0, 770)
    PyPlot.xlabel("Texture height [nm]", fontsize=18)
    PyPlot.ylabel("\$ \\eta_{\\mathrm{OC, NT}} - \\eta_{\\mathrm{OC, PL}} \$ [%]", fontsize=18)
    PyPlot.tick_params(which ="both", labelsize=18)
    PyPlot.legend(fontsize = 18)
    PyPlot.tight_layout()

    if saveFig
        savefig(datadir("ampl-efficiency-difference-params-$paramsname-generation-$generation-scanrate-$scanrate.pdf"))
    end

    ampl2 = collect(0.5:0.5:7.5) .* 1.0e-7; ampl = vcat(0.0, ampl2)

    col = ampl./1.0e-6

    textampl1 = "3p0e-7"
    textampl2 = "5p0e-7"
    textampl3 = "7p0e-7"

    Blues  = get_cmap(:Greys)
    Greens = get_cmap(:Greens)

    figure()
    IVPL  = readdlm(datadir("IV", "$path/IV-2D-forw-planar-generation-$generation-reco-radiative.dat"))
    IVNT1 = readdlm(datadir("IV", "$path/IV-2D-forw-nanotextured-ampl-$textampl1-generation-$generation-reco-radiative.dat"))
    IVNT2 = readdlm(datadir("IV", "$path/IV-2D-forw-nanotextured-ampl-$textampl2-generation-$generation-reco-radiative.dat"))
    IVNT3 = readdlm(datadir("IV", "$path/IV-2D-forw-nanotextured-ampl-$textampl3-generation-$generation-reco-radiative.dat"))

    PyPlot.plot(IVPL[:, 1],   -IVPL[:, 2].*(cm^2).*1.0e3./heightDev, linewidth = 5, color = Blues(251),  alpha = 0.8)
    PyPlot.plot(IVNT1[:, 1], -IVNT1[:, 2].*(cm^2).*1.0e3./heightDev, linewidth = 5, color = Blues(211),  alpha = 0.8)
    PyPlot.plot(IVNT2[:, 1], -IVNT2[:, 2].*(cm^2).*1.0e3./heightDev, linewidth = 5, color = Blues(161),  alpha = 0.8)
    PyPlot.plot(IVNT3[:, 1], -IVNT3[:, 2].*(cm^2).*1.0e3./heightDev, linewidth = 5, color = Blues(101),  alpha = 0.8)

    IVPL  = readdlm(datadir("IV", "$path/IV-2D-forw-planar-generation-$generation-reco-all.dat"))
    IVNT1 = readdlm(datadir("IV", "$path/IV-2D-forw-nanotextured-ampl-$textampl1-generation-$generation-reco-all.dat"))
    IVNT2 = readdlm(datadir("IV", "$path/IV-2D-forw-nanotextured-ampl-$textampl2-generation-$generation-reco-all.dat"))
    IVNT3 = readdlm(datadir("IV", "$path/IV-2D-forw-nanotextured-ampl-$textampl3-generation-$generation-reco-all.dat"))

    PyPlot.plot(IVPL[:, 1],   -IVPL[:, 2].*(cm^2).*1.0e3./heightDev, linewidth = 5, color = Greens(251), label = "planar")
    PyPlot.plot(IVNT1[:, 1], -IVNT1[:, 2].*(cm^2).*1.0e3./heightDev, linewidth = 5, color = Greens(211), label = "textured (300 nm)")
    PyPlot.plot(IVNT2[:, 1], -IVNT2[:, 2].*(cm^2).*1.0e3./heightDev, linewidth = 5, color = Greens(161), label = "textured (500 nm)")
    PyPlot.plot(IVNT3[:, 1], -IVNT3[:, 2].*(cm^2).*1.0e3./heightDev, linewidth = 5, color = Greens(101), label = "textured (700 nm)")

    PyPlot.xlim(0.1, 1.40)
    PyPlot.ylim(-2, 25)
    yticks([0, 5, 15, 25])
    PyPlot.legend(fontsize = 18)
    PyPlot.xlabel("Appl. voltage [V]", fontsize=18)
    PyPlot.ylabel("\$ J \$ [mA cm\$^{-2} \$]", fontsize=18)
    PyPlot.tick_params(which ="both", labelsize=18)
    PyPlot.tight_layout()
    PyPlot.xticks([0.6, 1.0, 1.4])

    if saveFig
        savefig(datadir("IV-planar-nanotexture-params-$paramsname-generation-$generation-scanrate-$scanrate.pdf"))
    end

    figure()
    IVPL  = readdlm(datadir("IV", "$path/IV-2D-forw-planar-generation-$generation-reco-radiative.dat"))
    IVNT1 = readdlm(datadir("IV", "$path/IV-2D-forw-nanotextured-ampl-$textampl1-generation-$generation-reco-radiative.dat"))
    IVNT2 = readdlm(datadir("IV", "$path/IV-2D-forw-nanotextured-ampl-$textampl2-generation-$generation-reco-radiative.dat"))
    IVNT3 = readdlm(datadir("IV", "$path/IV-2D-forw-nanotextured-ampl-$textampl3-generation-$generation-reco-radiative.dat"))

    PyPlot.plot(IVPL[:, 1],  -IVPL[:, 2].*(cm^2).*1.0e3./heightDev,  linewidth = 5, color = Blues(251), alpha = 0.8)
    PyPlot.plot(IVNT1[:, 1], -IVNT1[:, 2].*(cm^2).*1.0e3./heightDev, linewidth = 5, color = Blues(211), alpha = 0.8)
    PyPlot.plot(IVNT2[:, 1], -IVNT2[:, 2].*(cm^2).*1.0e3./heightDev, linewidth = 5, color = Blues(161), alpha = 0.8)
    PyPlot.plot(IVNT3[:, 1], -IVNT3[:, 2].*(cm^2).*1.0e3./heightDev, linewidth = 5, color = Blues(101), alpha = 0.8)

    IVPL  = readdlm(datadir("IV", "$path/IV-2D-forw-planar-generation-$generation-reco-all.dat"))
    IVNT1 = readdlm(datadir("IV", "$path/IV-2D-forw-nanotextured-ampl-$textampl1-generation-$generation-reco-all.dat"))
    IVNT2 = readdlm(datadir("IV", "$path/IV-2D-forw-nanotextured-ampl-$textampl2-generation-$generation-reco-all.dat"))
    IVNT3 = readdlm(datadir("IV", "$path/IV-2D-forw-nanotextured-ampl-$textampl3-generation-$generation-reco-all.dat"))

    PyPlot.plot(IVPL[:, 1],  -IVPL[:, 2].*(cm^2).*1.0e3./heightDev,  linewidth = 5, color = Greens(251), label = "planar")
    PyPlot.plot(IVNT1[:, 1], -IVNT1[:, 2].*(cm^2).*1.0e3./heightDev, linewidth = 5, color = Greens(211), label = "textured (300 nm)")
    PyPlot.plot(IVNT2[:, 1], -IVNT2[:, 2].*(cm^2).*1.0e3./heightDev, linewidth = 5, color = Greens(161), label = "textured (500 nm)")
    PyPlot.plot(IVNT3[:, 1], -IVNT3[:, 2].*(cm^2).*1.0e3./heightDev, linewidth = 5, color = Greens(101), label = "textured (700 nm)")

    PyPlot.xlim(1.05, 1.40)
    if generation == "uniform"
        PyPlot.ylim(-2, 10)
        yticks([0, 5, 10])
    else
        PyPlot.ylim(-2, 16)
        yticks([0, 5, 10, 15])
    end
    PyPlot.legend(fontsize = 18)
    PyPlot.xlabel("Appl. voltage [V]", fontsize=18)
    PyPlot.ylabel("\$ J \$ [mA cm\$^{-2} \$]", fontsize=18)
    PyPlot.tick_params(which ="both", labelsize=18)
    PyPlot.tight_layout()
    PyPlot.xticks([1.1, 1.3])

    if saveFig
        savefig(datadir("IV-planar-nanotexture-zoom-params-$paramsname-generation-$generation-scanrate-$scanrate.pdf"))
    end

    # figure()
    # A = rand(20,20)
    # figure()
    # PyPlot.pcolor(A, cmap="Greys")
    # colorbar(ticks = [0.0, 0.25, 0.5, 0.75, 1.0])

    # if saveFig
    #     savefig(datadir("colorbar-grey.pdf"))
    # end

    #################################################################
    #################################################################

    JSCOptical = readdlm(datadir("Max-PhotoCurrent-Optical.dat"))

    figure()
    if generation == "Maxwell"
        PyPlot.plot(ampl./nm, JSCOptical[1:end-1], "o", markersize = MKsize, color = [255/255, 118/255, 118/255], label = "optical maximum")
    end
    PyPlot.plot(ampl./nm, JSCVecRad,           "o", markersize = MKsize, color = ColBlue,   label = "ideal")
    PyPlot.plot(ampl./nm, JSCVec,              "o", markersize = MKsize, color = ColGold,   label = "realistic")

    PyPlot.xticks([0, 300, 600])
    PyPlot.xlim(-20.0, 770)
    PyPlot.ylim(21, 24)
    PyPlot.xlabel("Texture height [nm]", fontsize=18)
    PyPlot.ylabel("\$ J_{\\mathrm{SC}}\$ [mA cm\$^{-2} \$]", fontsize=18)
    PyPlot.tick_params(which ="both", labelsize=18)
    PyPlot.legend(fontsize = 18)
    PyPlot.tight_layout()

    if saveFig
        savefig(datadir("ampl-JSC-params-$paramsname-generation-$generation-scanrate-$scanrate.pdf"))
    end

    if printText
        println(".................................................")
        println(" ")
        println("Planar (realistic): JSC is $(JSCVec[1]) mA/cm^2.")
        println("Planar (ideal):     JSC is $(JSCVecRad[1]) mA/cm^2.")

        println(" ")

        println("Realistic: JSC is $(JSCVec[Max2]) mA/cm^2.")
        println("Ideal:     JSC is $(JSCVecRad[Max2Rad]) mA/cm^2.")
    end

    ##################################################################
    figure()
    PyPlot.plot(ampl./nm, FFVecRad, "o", markersize = MKsize, color = ColBlue, label = "ideal")
    PyPlot.plot(ampl./nm, FFVec,    "o", markersize = MKsize, color = ColGold, label = "realistic")

    PyPlot.xticks([0, 300, 600])
    PyPlot.xlim(-20.0, 770)
    PyPlot.xlabel("Texture height [nm]", fontsize=18)
    PyPlot.ylabel("FF [%]", fontsize=18)
    PyPlot.tick_params(which ="both", labelsize=18)
    PyPlot.legend(fontsize = 18)
    PyPlot.tight_layout()

    if saveFig
        savefig(datadir("ampl-FF-params-$paramsname-generation-$generation-scanrate-$scanrate.pdf"))
    end

    if printText
        println(".................................................")
        println(" ")
        println("Planar (realistic): FF is $(FFVec[1]) %.")
        println("Planar (ideal):     FF is $(FFVecRad[1]) %.")

        println(" ")

        println("Realistic: FF is $(FFVec[Max2]) %.")
        println("Ideal:     FF is $(FFVecRad[Max2Rad]) %.")
    end

    #################################

    figure()
    PyPlot.plot(ampl./nm, VOCVecRad, "o", markersize = MKsize, color = ColBlue, label = "ideal")
    PyPlot.plot(ampl./nm, VOCVec,    "o", markersize = MKsize, color = ColGold, label = "realistic")

    PyPlot.xticks([0, 300, 600])
    PyPlot.xlim(-20.0, 770)
    PyPlot.ylim(1.16, 1.37)
    PyPlot.xlabel("Texture height [nm]", fontsize=18)
    PyPlot.ylabel("\$ V_{\\mathrm{OC}} \$ [V]", fontsize=18)
    PyPlot.tick_params(which ="both", labelsize=18)
    PyPlot.legend(fontsize = 18)
    PyPlot.tight_layout()

    if saveFig
        savefig(datadir("ampl-VOC-params-$paramsname-generation-$generation-scanrate-$scanrate.pdf"))
    end

    if printText
        println(".................................................")
        println(" ")
        println("Planar (realistic): VOC is $(VOCVec[1]) V.")
        println("Planar (ideal):     VOC is $(VOCVecRad[1]) V.")

        println(" ")

        println("Realistic: VOC is $(VOCVec[Max2]) V.")
        println("Ideal:     VOC is $(VOCVecRad[Max2Rad]) V.")
    end

    #################################
    VOCDiff = zeros(0); VOCDiffRad = zeros(0)

    for ii = 2:length(VOCVec)
        VOCdiff    = VOCVec[ii]    - VOCVec[1]
        VOCdiffRad = VOCVecRad[ii] - VOCVecRad[1]

        push!(VOCDiff, VOCdiff); push!(VOCDiffRad, VOCdiffRad)
    end

    figure()
    PyPlot.plot(ampl[2:end]./nm, VOCDiffRad.*1.0e3, "o", markersize = MKsize, color = ColBlue, label = "ideal")
    PyPlot.plot(ampl[2:end]./nm, VOCDiff.*1.0e3,    "o", markersize = MKsize, color = ColGold, label = "realistic")

    PyPlot.xticks([0, 300, 600])
    PyPlot.xlim(-20.0, 770)
    PyPlot.ylim(-1, 21)
    PyPlot.xlabel("Texture height [nm]", fontsize=18)
    PyPlot.ylabel("\$ V_{\\mathrm{OC, NT}} - V_{\\mathrm{OC, PL}} \$ [mV]", fontsize=18)
    PyPlot.tick_params(which ="both", labelsize=18)
    PyPlot.legend(fontsize = 18)
    PyPlot.tight_layout()

    if saveFig
        savefig(datadir("ampl-VOC-difference-params-$paramsname-generation-$generation-scanrate-$scanrate.pdf"))
    end

    # #################################
    # JMPVMPDiff = zeros(0); JSCVOCDiff = zeros(0)

    # for ii = 1:length(JMP)
    #     JMPVMPdiff = JMP[ii] * VMP[ii] #- JMP[1] * VMP[1]
    #     JSCVOCdiff = JSCVec[ii] * VOCVec[ii] #- JSCVec[1] * VOCVec[1]

    #     push!(JMPVMPDiff, JMPVMPdiff); push!(JSCVOCDiff, JSCVOCdiff)

    # end

    # figure()
    # PyPlot.plot(ampl./nm, 1.0e-1.*JSCVOCDiff, "o", markersize = MKsize, color = "midnightblue", label = "\$ J_{\\mathrm{SC}} V_{\\mathrm{OC}} \$")
    # PyPlot.plot(ampl./nm, 1.0e-1.*JMPVMPDiff, "o", markersize = MKsize, color = "dodgerblue", label = "\$ J_{\\mathrm{MP}} V_{\\mathrm{MP}} \$")

    # PyPlot.xticks([0, 300, 600])
    # PyPlot.xlim(-20.0, 770)
    # #PyPlot.ylim(18, 21)
    # PyPlot.xlabel("Texture height [nm]", fontsize=18)
    # PyPlot.ylabel("Power [W]", fontsize=18)
    # PyPlot.tick_params(which ="both", labelsize=18)
    # PyPlot.legend(fontsize = 18)
    # PyPlot.tight_layout()

    # if saveFig
    #     savefig(datadir("ampl-JMP-VMP-JSC-VOC-params-$paramsname-generation-$generation-scanrate-$scanrate.pdf"))
    # end


    return nothing
end

end