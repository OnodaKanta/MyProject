# test/runtests.jl

using RenormalizationSim
using Test

@testset "RenormalizationSim.jl" begin
    # --- calculate_sigma 関数のテスト ---
    @testset "calculate_sigma" begin
        # 修正点: 呼び出しを3引数に変更 (m_squared_phys, lambda_coupling, cutoff_Lambda)
        # 結合定数が0なら、自己エネルギーは0になるはず
        @test calculate_sigma(1.0, 0.0, 100.0) == 0.0

        # 修正点: 呼び出しを3引数に変更
        # 簡単なケースで型が正しいかテスト
        sigma = calculate_sigma(64.0, 20.0, 100.0)
        @test sigma isa Float64

        # 結果が正の値であるべきことをテスト
        @test sigma > 0.0
    end

    # --- シミュレーション実行のテスト ---
    # このテストは元々成功していたので、変更は不要です。
    @testset "run_full_simulation" begin
        println("\nテスト: run_full_simulationがエラーなく実行されるか確認...")
        # run_full_simulation() がエラーを投げずに正常終了すればテスト成功とみなす
        @test begin
            run_full_simulation()
            true # エラーがなければ true を返す
        end
    end
end
