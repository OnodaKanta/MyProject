# test/runtests.jl

using RenormalizationSim
using Test

@testset "RenormalizationSim.jl" begin
    
    @testset "calculate_sigma" begin
        
        @test calculate_sigma(1.0, 0.0, 100.0) == 0.0

        sigma = calculate_sigma(64.0, 20.0, 100.0)
        @test sigma isa Float64

        @test sigma > 0.0
    end

   
    @testset "run_full_simulation" begin
        println("\nテスト: run_full_simulationがエラーなく実行されるか確認...")
       
        @test begin
            run_full_simulation()
            true 
        end
    end
end