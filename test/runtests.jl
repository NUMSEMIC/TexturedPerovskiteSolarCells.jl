using Aqua
using Test
using TexturedPerovskiteSolarCells

@testset "Aqua - Quality Test" begin
    Aqua.test_all(
        TexturedPerovskiteSolarCells;
        stale_deps = (ignore = [:DelimitedFiles, :PyPlot, :PyCall],),
    )
end


@testset "SingleJunction" begin

    include(joinpath("..", "scripts", "SingleJunction.jl"))

    @test SingleJunction.test(gridDim = 1, demo_run = true) == true
    @test SingleJunction.test(gridDim = 2, typeGrid = "nanotextured", amplitude = 2.0e-7, demo_run = true) == true

end


@testset "Figures" begin

    if !Sys.iswindows()
        for Fig in [
                :Fig2Photogeneration,
                :Fig3CharacteristicsStudy,
                :Fig4RecombinationCurrents,
                :Fig5Density1D,
                :Fig5ElectricField,
                :Fig5RatioElectronsAndHoles,
                :FigS5VacancyDensity,
                :FigS6ElectronAndHoleDensity,
            ]

            @info "plot $Fig"

            include(joinpath("..", "PostProcess", string(Fig) * ".jl"))

            @eval begin
                @test $(Fig).main(printText = false, saveFig = true) === nothing
            end
        end
    end
end
