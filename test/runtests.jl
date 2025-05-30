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

    @test SingleJunction.test(demo_run = true) == true

end

