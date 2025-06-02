#=

Code for visualizing Recombination currents

=#

module Fig4RecombinationCurrents

using PyPlot
using DelimitedFiles
using ChargeTransport
using PyCall
using TexturedPerovskiteSolarCells

# for convenience
datadir = TexturedPerovskiteSolarCells.datadir
scriptsdir = TexturedPerovskiteSolarCells.scriptsdir

# https://stackoverflow.com/questions/29443369/how-to-make-a-custom-colormap-using-pyplot-not-matplotlib-proper
@pyimport matplotlib.colors as matcolors

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

function main(;printText = true, saveFig   = false,
              scanrate   = "1000p0",
              generation = "Maxwell", # "uniform"
              parameter_file = scriptsdir("params_single_junction.jl"),
              )

    include(parameter_file)

    PyPlot.rc("font", family="sans-serif", size=14)
    PyPlot.rc("mathtext", fontset="dejavusans")
    PyPlot.close("all")

    path = "parameter-$paramsname/$scanrate"

    Textampl = ["5p0e-8", "1p0e-7", "1p5e-7", "2p0e-7", "2p5e-7", "3p0e-7", "3p5e-7", "4p0e-7",
    "4p5e-7", "5p0e-7", "5p5e-7", "6p0e-7", "6p5e-7", "7p0e-7", "7p5e-7"]

    ampl2    = collect(0.5:0.5:7.5) .* 1.0e-7; ampl = vcat(0.0, ampl2)

    col = ampl./1.0e-6

    #############################################################
    if printText
        println("Plot JSRH")
    end

    IVPL = readdlm(datadir("IV", "$path/JSRH-2D-forw-planar-generation-$generation-reco-all.dat"))

    semilogy(IVPL[:, 1],  abs.(IVPL[:, 2].*(cm^2).*1.0e3./heightDev),   linewidth = 5, color = parula_map(col[1]))

    ii = 0
    for textampl in Textampl[1:2:end]
        ii = ii + 2

        if printText
            println("amplitude = ", textampl)
        end

        IVNT = readdlm(datadir("IV", "$path/JSRH-2D-forw-nanotextured-ampl-$textampl-generation-$generation-reco-all.dat"))

        semilogy(IVNT[:, 1],  abs.(IVNT[:, 2].*(cm^2).*1.0e3./heightDev),   linewidth = 5, color = parula_map(col[ii]))
    end

    xlim(0.0, 1.2)
    ylim(2.0e-7, 1.9e1)
    PyPlot.xlabel("Appl. voltage [V]", fontsize=18)
    PyPlot.ylabel("\$ J \$ [mA cm\$^{-2} \$]", fontsize=18)
    PyPlot.tick_params(which ="both", labelsize=18)
    PyPlot.title("JSRH")
    PyPlot.tight_layout()

    if saveFig
        savefig(datadir("J-reco-SRH-generation-$generation-scanrate-$scanrate.pdf"))
    end

    figure()
    IVPL = readdlm(datadir("IV", "$path/JSRH-2D-forw-planar-generation-$generation-reco-all.dat"))
    semilogy(IVPL[:, 1],  abs.(IVPL[:, 2].*(cm^2).*1.0e3./heightDev),   linewidth = 5, color = parula_map(col[1]))

    textampl = "3p0e-7"; ii = 7
    IVNT     = readdlm(datadir("IV", "$path/JSRH-2D-forw-nanotextured-ampl-$textampl-generation-$generation-reco-all.dat"))
    semilogy(IVNT[:, 1],  abs.(IVNT[:, 2].*(cm^2).*1.0e3./heightDev),   linewidth = 5, color = parula_map(col[ii]))

    textampl = "5p0e-7"; ii = 11
    IVNT     = readdlm(datadir("IV", "$path/JSRH-2D-forw-nanotextured-ampl-$textampl-generation-$generation-reco-all.dat"))
    semilogy(IVNT[:, 1],  abs.(IVNT[:, 2].*(cm^2).*1.0e3./heightDev),   linewidth = 5, color = parula_map(col[ii]))

    textampl = "7p0e-7"; ii = 15
    IVNT     = readdlm(datadir("IV", "$path/JSRH-2D-forw-nanotextured-ampl-$textampl-generation-$generation-reco-all.dat"))
    semilogy(IVNT[:, 1],  abs.(IVNT[:, 2].*(cm^2).*1.0e3./heightDev),   linewidth = 5, color = parula_map(col[ii]))

    xlim(0.97, 1.20)
    xticks([1.0, 1.1, 1.2])
    ylim(5.0e-1, 1.9e1)
    PyPlot.xlabel("Appl. voltage [V]", fontsize=18)
    PyPlot.ylabel("\$ J\$ [mA cm\$^{-2} \$]", fontsize=18)
    PyPlot.tick_params(which ="both", labelsize=18)
    PyPlot.title("JSRH")
    PyPlot.tight_layout()

    if saveFig
        savefig(datadir("J-reco-SRH-zoom-generation-$generation-scanrate-$scanrate.pdf"))
    end

    ############################################################
    if printText
        println("  ")
        println("Plot JRad")
    end

    IVPL = readdlm(datadir("IV", "$path/JRad-2D-forw-planar-generation-$generation-reco-all.dat"))

    figure()
    semilogy(IVPL[:, 1],  abs.(IVPL[:, 2].*(cm^2).*1.0e3./heightDev),   linewidth = 5, color = parula_map(col[1]))

    ii = 0
    for textampl in Textampl[1:2:end]
        ii = ii + 2

        if printText
            println("amplitude = ", textampl)
        end

        IVNT = readdlm(datadir("IV", "$path/JRad-2D-forw-nanotextured-ampl-$textampl-generation-$generation-reco-all.dat"))

        semilogy(IVNT[:, 1],  abs.(IVNT[:, 2].*(cm^2).*1.0e3./heightDev),   linewidth = 5, color = parula_map(col[ii]))
    end

    xlim(0.0, 1.2)
    ylim(2.0e-7, 1.9e1)
    PyPlot.xlabel("Appl. voltage [V]", fontsize=18)
    PyPlot.ylabel("\$ J\$ [mA cm\$^{-2} \$]", fontsize=18)
    PyPlot.tick_params(which ="both", labelsize=18)
    PyPlot.title("Jrad")
    PyPlot.tight_layout()

    if saveFig
        savefig(datadir("J-reco-rad-generation-$generation-scanrate-$scanrate.pdf"))
    end

    ############################################################
    if printText
        println("  ")
        println("Plot JSR, ETL")
    end

    ## left boundary
    IVPL = readdlm(datadir("IV", "$path/JSRnL-2D-forw-planar-generation-$generation-reco-all.dat"))

    figure()
    semilogy(IVPL[:, 1],  abs.(IVPL[:, 2].*(cm^2).*1.0e3./heightDev),   linewidth = 5, color = parula_map(col[1]))

    ii = 0
    for textampl in Textampl[1:2:end]
        ii = ii + 2

        if printText
            println("amplitude = ", textampl)
        end

        IVNT = readdlm(datadir("IV", "$path/JSRnL-2D-forw-nanotextured-ampl-$textampl-generation-$generation-reco-all.dat"))

        semilogy(IVNT[:, 1],  abs.(IVNT[:, 2].*(cm^2).*1.0e3./heightDev),   linewidth = 5, color = parula_map(col[ii]))
    end

    xlim(0.0, 1.2)
    ylim(2.0e-7, 1.9e1)
    PyPlot.xlabel("Appl. voltage [V]", fontsize=18)
    PyPlot.ylabel("\$ J \$ [mA cm\$^{-2} \$]", fontsize=18)
    PyPlot.tick_params(which ="both", labelsize=18)
    PyPlot.title("JSR ETL")
    PyPlot.tight_layout()

    if saveFig
        savefig(datadir("J-reco-SR-ETL-generation-$generation-params-$paramsname-scanrate-$scanrate.pdf"))
    end

    ############################################################
    if printText
        println("  ")
        println("Plot JSR, HTL")
    end

    ## right boundary
    IVPL = readdlm(datadir("IV", "$path/JSRnR-2D-forw-planar-generation-$generation-reco-all.dat"))

    figure()
    semilogy(IVPL[:, 1],  abs.(IVPL[:, 2].*(cm^2).*1.0e3./heightDev),   linewidth = 5, color = parula_map(col[1]))

    ii = 0
    for textampl in Textampl[1:2:end]
        ii = ii + 2

        if printText
            println("amplitude = ", textampl)
        end

        IVNT = readdlm(datadir("IV", "$path/JSRnR-2D-forw-nanotextured-ampl-$textampl-generation-$generation-reco-all.dat"))

        semilogy(IVNT[:, 1],  abs.(IVNT[:, 2].*(cm^2).*1.0e3./heightDev),   linewidth = 5, color = parula_map(col[ii]))
    end

    xlim(0.0, 1.2)
    ylim(2.0e-7, 1.9e1)
    PyPlot.xlabel("Appl. voltage [V]", fontsize=18)
    PyPlot.ylabel("\$ J \$ [mA cm\$^{-2} \$]", fontsize=18)
    PyPlot.tick_params(which ="both", labelsize=18)
    PyPlot.title("JSR HTL")
    PyPlot.tight_layout()

    if saveFig
        savefig(datadir("J-reco-SR-HTL-generation-$generation-params-$paramsname-scanrate-$scanrate.pdf"))
    end

    return nothing
end

end