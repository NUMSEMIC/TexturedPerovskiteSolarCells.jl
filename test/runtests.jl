using Aqua
using Test
using TexturedPerovskiteSolarCells


@testset "Aqua - Quality Test" begin
    Aqua.test_all(
        TexturedPerovskiteSolarCells;
        stale_deps = (ignore = [:DelimitedFiles, :PyPlot],),
    )
end


@testset "SingleJunction" begin

    include(joinpath("..", "scripts", "SingleJunction.jl"))

    @test SingleJunction.test(gridDim =1, demo_run = true) == true
    @test SingleJunction.test(gridDim = 2, typeGrid = "nanotextured", amplitude = 2.0e-7, demo_run = true) == true

end

