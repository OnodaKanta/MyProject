module RenormalizationSim

using Plots
using QuadGK
using Printf

export calculate_sigma, run_full_simulation

"""
スカラー場の1ループ自己エネルギー Σ(p^2) を計算します。

"""
function calculate_sigma(p_squared_euclidean::Float64, m_squared_phys::Float64, lambda_coupling::Float64, cutoff_Lambda::Float64)
    prefactor = lambda_coupling / (16.0 * π^2)
    integrand(k) = k^3 / (k^2 + m_squared_phys)
    integral_val, _ = quadgk(integrand, 0, cutoff_Lambda)
    return prefactor * integral_val
end

"""
シミュレーション全体を実行し、結果をプロットします。

"""
function run_full_simulation()
    m_phys = 7.0
    m_squared_phys = m_phys^2
    lambda_coupling = 20.0
    epsilon = 1e-6

    # --- グラフ1: 繰り込みフローを計算 ---
    println("Graph 1: Calculating renormalization flow...")
    cutoffs_for_running = 1.0:1.0:1000.0
    bare_mass_sq_values = []
    sigma_values = []

    for Lambda in cutoffs_for_running
        sigma_at_mass_pole = calculate_sigma(m_squared_phys, m_squared_phys, lambda_coupling, Lambda)
        bare_mass_sq = m_squared_phys - sigma_at_mass_pole
        push!(sigma_values, sigma_at_mass_pole)
        push!(bare_mass_sq_values, bare_mass_sq)
    end
    println("Calculation complete.")

    # --- グラフ2: 伝播関数の安定性を計算 ---
    println("\nGraph 2: Calculating propagator stability...")
    p_squared_minkowski_range = 0.0:0.02:100.0
    p_squared_euclidean_range = -p_squared_minkowski_range
    low_cutoff = 1.0
    high_cutoff = 10000.0
    
    sigma_low = calculate_sigma(m_squared_phys, m_squared_phys, lambda_coupling, low_cutoff)
    m0_sq_low = m_squared_phys - sigma_low
    println("Bare mass for low cutoff (Λ=$(low_cutoff)): m₀² = $(round(m0_sq_low, digits=4))")
    
    sigma_high = calculate_sigma(m_squared_phys, m_squared_phys, lambda_coupling, high_cutoff)
    m0_sq_high = m_squared_phys - sigma_high
    println("Bare mass for high cutoff (Λ=$(high_cutoff)): m₀² = $(round(m0_sq_high, digits=4))")

    propagator_free = []
    propagator_low_cutoff = []
    propagator_high_cutoff = []

    for (i, p_sq_minkowski) in enumerate(p_squared_minkowski_range)
        p_sq_euclidean = p_squared_euclidean_range[i]

        # 1. 自由な伝播関数
        D_free_inv = p_sq_minkowski - m_squared_phys + im * epsilon
        push!(propagator_free, abs2(1.0 / D_free_inv))

        # 2. 低いカットオフを持つ裸でない伝播関数
        sigma_p_low = calculate_sigma(p_sq_euclidean, m_squared_phys, lambda_coupling, low_cutoff)
        D_low_inv = p_sq_minkowski - m0_sq_low - sigma_p_low + im * epsilon
        push!(propagator_low_cutoff, abs2(1.0 / D_low_inv))

        # 3. 高いカットオフを持つ裸でない伝播関数
        sigma_p_high = calculate_sigma(p_sq_euclidean, m_squared_phys, lambda_coupling, high_cutoff)
        D_high_inv = p_sq_minkowski - m0_sq_high - sigma_p_high + im * epsilon
        push!(propagator_high_cutoff, abs2(1.0 / D_high_inv))
    end
    println("Calculation complete.")

    # --- 結果をプロット ---
    plot1 = plot(
        cutoffs_for_running,
        [sigma_values, bare_mass_sq_values],
        title="Graph 1: Renormalization Flow",
        xlabel="Momentum Cutoff Λ (GeV)",
        ylabel="Mass Squared (GeV²)",
        label=["Self-Energy Correction Σ(Λ)" "Bare Mass Squared m₀²(Λ)"],
        linewidth=2,
        legend=:bottomleft
    )
    hline!(plot1, [m_squared_phys], linestyle=:dash, color=:black, label="Physical Mass Squared m_phys²")

    plot2 = plot(
        p_squared_minkowski_range,
        [propagator_free, propagator_low_cutoff, propagator_high_cutoff],
        title="Graph 2: Propagator Stability",
        xlabel="Momentum Squared p² (GeV²)",
        ylabel="Propagator Absolute Value Squared |D'(p)|²",
        label=["Free Propagator" "Dressed Propagator (Low Λ)" "Dressed Propagator (High Λ)"],
        linewidth=[2 2 2],
        linestyle=[:dot :solid :dash],
        legend=:topleft,
        yaxis=:log
    )
    vline!(plot2, [m_squared_phys], linestyle=:dash, color=:black, label="Physical Mass m_phys²")

    final_plot = plot(plot1, plot2, layout=(2, 1), size=(800, 800))
    savefig(final_plot, "renormalization_plots.png")
end

end
