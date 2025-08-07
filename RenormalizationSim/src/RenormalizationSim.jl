module RenormalizationSim

using Plots
using Measures 
gr()
using QuadGK
using Printf

export calculate_sigma, run_full_simulation

function calculate_sigma(m_squared_phys::Float64, lambda_coupling::Float64, cutoff_Lambda::Float64)
    prefactor = lambda_coupling / (16.0 * π^2)
    integrand(k) = k^3 / (k^2 + m_squared_phys)
    integral_val, _ = quadgk(integrand, 0, cutoff_Lambda)
    return prefactor * integral_val
end

function calculate_sigma_derivative(m_squared_phys::Float64, lambda_coupling::Float64, cutoff_Lambda::Float64)
   
    prefactor = -lambda_coupling / (16.0 * π^2)
    integrand(k) = k^3 / (6.0 * (k^2 + m_squared_phys)^2)
    integral_val, _ = quadgk(integrand, 0, cutoff_Lambda)
    return prefactor * integral_val
end

function run_full_simulation()
  
    m_phys = 7.0               
    m_squared_phys = m_phys^2   
    lambda_coupling = 20.0      
    epsilon = 1e-6              

    # --- グラフ1 & 2: 繰り込みフローを計算 ---
    println("Graph 1 & 2: Calculating renormalization flow for m₀²(Λ) and Z(Λ)...")
    cutoffs_for_running = 1.0:1.0:1000.0
    bare_mass_sq_values = []
    sigma_values = []
    z_factor_values = []

    for Lambda in cutoffs_for_running
        
        sigma_at_mass_pole = calculate_sigma(m_squared_phys, lambda_coupling, Lambda)
        bare_mass_sq = m_squared_phys - sigma_at_mass_pole
        push!(sigma_values, sigma_at_mass_pole)
        push!(bare_mass_sq_values, bare_mass_sq)

        sigma_deriv = calculate_sigma_derivative(m_squared_phys, lambda_coupling, Lambda)
        z_factor = 1.0 / (1.0 - sigma_deriv)
        push!(z_factor_values, z_factor)
    end
    println("Calculation complete.")

    # --- グラフ3: 伝播関数の安定性を計算 ---
    println("\nGraph 3: Calculating propagator stability...")
    p_squared_minkowski_range = 0.0:0.02:100.0
    low_cutoff = 10.0
    high_cutoff = 1000.0
    
    # 低いカットオフでのパラメータ
    sigma_low = calculate_sigma(m_squared_phys, lambda_coupling, low_cutoff)
    m0_sq_low = m_squared_phys - sigma_low
    sigma_deriv_low = calculate_sigma_derivative(m_squared_phys, lambda_coupling, low_cutoff)
    z_low = 1.0 / (1.0 - sigma_deriv_low)
    @printf("Low cutoff (Λ=%.1f): m₀²=%.4f, Z=%.4f\n", low_cutoff, m0_sq_low, z_low)
    
    # 高いカットオフでのパラメータ
    sigma_high = calculate_sigma(m_squared_phys, lambda_coupling, high_cutoff)
    m0_sq_high = m_squared_phys - sigma_high
    sigma_deriv_high = calculate_sigma_derivative(m_squared_phys, lambda_coupling, high_cutoff)
    z_high = 1.0 / (1.0 - sigma_deriv_high)
    @printf("High cutoff (Λ=%.1f): m₀²=%.4f, Z=%.4f\n", high_cutoff, m0_sq_high, z_high)

    propagator_free = []
    propagator_low_cutoff = []
    propagator_high_cutoff = []

    for p_sq_minkowski in p_squared_minkowski_range
        # 1. 自由な伝播関数 (Z=1)
        D_free_inv = p_sq_minkowski - m_squared_phys + im * epsilon
        push!(propagator_free, abs2(1.0 / D_free_inv))

        # 2. 低いカットオフを持つ繰り込み済み伝播関数
        sigma_p_low = calculate_sigma(m_squared_phys, lambda_coupling, low_cutoff) 
        D_low_inv = p_sq_minkowski - m0_sq_low - sigma_p_low + im * epsilon
        push!(propagator_low_cutoff, abs2(z_low / D_low_inv)) 

        # 3. 高いカットオフを持つ繰り込み済み伝播関数
        sigma_p_high = calculate_sigma(m_squared_phys, lambda_coupling, high_cutoff) 
        D_high_inv = p_sq_minkowski - m0_sq_high - sigma_p_high + im * epsilon
        push!(propagator_high_cutoff, abs2(z_high / D_high_inv)) 
    end
    println("Calculation complete.")

    # --- グラフ4: 3D伝播関数サーフェスを計算 ---
    println("\nGraph 4: Calculating 3D propagator surface...")
    p_squared_range_3d = 0.0:2.0:100.0
    cutoff_range_3d = 10.0:10.0:500.0
    propagator_surface = zeros(length(cutoff_range_3d), length(p_squared_range_3d))

    for (i, Lambda) in enumerate(cutoff_range_3d)
        
        sigma = calculate_sigma(m_squared_phys, lambda_coupling, Lambda)
        m0_sq = m_squared_phys - sigma
        sigma_deriv = calculate_sigma_derivative(m_squared_phys, lambda_coupling, Lambda)
        z = 1.0 / (1.0 - sigma_deriv)
        
        for (j, p_sq) in enumerate(p_squared_range_3d)
            
            sigma_p = sigma 
            
            D_inv = p_sq - m0_sq - sigma_p + im * epsilon
            propagator_surface[i, j] = abs2(z / D_inv)
        end
    end
    log_propagator_surface = log10.(propagator_surface)

    finite_values = filter(isfinite, log_propagator_surface)

    if !isempty(finite_values)
        min_val = minimum(finite_values)
        replace!(log_propagator_surface, NaN => min_val - 1)
        replace!(log_propagator_surface, -Inf => min_val - 1)
    else
        replace!(log_propagator_surface, NaN => 0)
        replace!(log_propagator_surface, -Inf => 0)
        replace!(log_propagator_surface, Inf => 0)
        println("Warning: Could not find any finite values for the 3D plot surface. The plot may be empty.")
    end

    println("Calculation complete.")

    
    plot1 = plot(
        cutoffs_for_running,
        [sigma_values, bare_mass_sq_values],
        title="Graph 1: Renormalization Flow of Mass",
        xlabel="Momentum Cutoff Λ (GeV)",
        ylabel="Mass Squared (GeV²)",
        label=["Self-Energy Σ(Λ)" "Bare Mass m₀²(Λ)"],
        linewidth=2,
        legend=:bottomleft
    )
    hline!(plot1, [m_squared_phys], linestyle=:dash, color=:black, label="Physical Mass m_phys²")

    plot2 = plot(
        cutoffs_for_running,
        z_factor_values,
        title="Graph 2: Renormalization Flow of Z-factor",
        xlabel="Momentum Cutoff Λ (GeV)",
        ylabel="Field Strength Renormalization Z(Λ)",
        label="Z(Λ)",
        linewidth=2,
        legend=:topright 
    )
    hline!(plot2, [1.0], linestyle=:dash, color=:black, label="Z = 1 (Free Theory)")

    plot3 = plot(
        p_squared_minkowski_range,
        [propagator_free, propagator_low_cutoff, propagator_high_cutoff],
        title="Graph 3: Dressed Propagator",
        xlabel="Momentum Squared p² (GeV²)",
        ylabel="Propagator |D'(p)|²",
        label=["Free Propagator" "Dressed Propagator (Low Λ)" "Dressed Propagator (High Λ)"],
        linewidth=[2 2 2],
        linestyle=[:dot :solid :dash],
        legend=:topleft,
        yaxis=:log
    )
    vline!(plot3, [m_squared_phys], linestyle=:dash, color=:black, label="Physical Mass Pole")

    plot4 = surface(
        p_squared_range_3d,
        cutoff_range_3d,
        log_propagator_surface,
        title="Graph 4: Propagator vs p² and Cutoff Λ",
        xlabel="Momentum Squared p² (GeV²)",
        ylabel="Cutoff Λ (GeV)",
        zlabel="log₁₀|D'(p,Λ)|²", 
        camera=(30, 45), 
        colorbar_title="log₁₀|D'|²", 
        xrotation = -30, 
        yrotation = 45  
    )

    final_plot = plot(
        plot1, plot2, plot3, plot4, 
        layout=(4, 1), 
        size=(900, 1800),
        left_margin=20mm, 
        bottom_margin=20mm 
    )
    savefig(final_plot, "renormalization_plots_advanced.png")
    println("\nPlots saved to 'renormalization_plots_advanced.png'")
end

end
