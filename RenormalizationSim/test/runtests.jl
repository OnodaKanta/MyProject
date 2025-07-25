# test/runtests.jl

using RenormalizationSim 
using Test              

@testset "RenormalizationSim.jl" begin
    # --- calculate_sigma 関数のテスト ---
    @testset "calculate_sigma" begin
        # 結合定数が0なら、自己エネルギーは0になるはず
        @test calculate_sigma(1.0, 1.0, 0.0, 100.0) == 0.0

        # 簡単なケースで型が正しいかテスト
        sigma = calculate_sigma(64.0, 64.0, 20.0, 100.0)
        @test sigma isa Float64

        # 結果が正の値であるべきことをテスト
        @test sigma > 0.0
    end

    # --- シミュレーション実行のテスト ---
    @testset "run_full_simulation" begin
        println("\nテスト: run_full_simulationがエラーなく実行されるか確認...")
        @test begin
            run_full_simulation()
            true 
        end
    end
end
