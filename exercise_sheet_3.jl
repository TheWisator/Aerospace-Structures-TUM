using Plots, LinearAlgebra, Unitful, DataFrames
#plotly()

## Material Definition of Ti6Al4V
E_ti = 114_000
ν_ti = 0.33
R_p02_ti = 830
ρ_ti = 4_430

## Material Definition of Al 7075
E_al = 70_000
R_p02_al = 460
ρ_al = 2_800

## Material Definition of stainless steel
E_ss = 210_000
R_p02_ss = 960
ρ_ss = 7_800

## Material Definition of CFRP
E_cf = 90_000
R_p02_cf = 800
ρ_cf = 1_500

## Definition Case 1: Orbiter / B747
L_1 = 1796.6
k_c1 = (E_ti * A_c1 / L_1)
Fx_c1 = 90_300
Fy_c1 = 154_800
ϕ_c1 = 52.8 * π / 180

## Definition Case 2: Orbiter / B747
L_2 = 1796.6
k_c2 = (E_ti * A_c2 / L_2)
Fx_c2 = 404_800
Fy_c2 = 529_300
ϕ_c2 = 52.8 * π / 180



# define Master Stiffness Matrix
K(p1, p2, k1, k2) = ustrip.([   k1 .* cos.(p1).^2 -k1 .* sin.(p1) .* cos.(p1) -k1 .* cos.(p1).^2 k1 .* sin.(p1) .* cos.(p1) 0 0;
                -k1 .* sin.(p1) .* cos.(p1) k1 .* sin.(p1).^2 k1 .* sin.(p1) .* cos.(p1) -k1 .* sin.(p1).^2 0 0;
                -k1 .* cos.(p1).^2 k1 .* sin.(p1) .* cos.(p1) k1 .* cos.(p1).^2 .+ k2 .* cos.(p2).^2 -k1 .* sin.(p1) .* cos.(p1) + k2 .* sin.(p2) .* cos.(p2) -k2 .* cos.(p2).^2 -k2.* sin.(p2) .* cos.(p2);
                k1 .* sin.(p1) .* cos.(p1) -k1 .* sin.(p1).^2 -k1 .* sin.(p1) .* cos.(p1) .+ k2 .* sin.(p2) .* cos.(p2) k1 .* sin.(p1).^2 .+ k2 .* sin.(p2).^2 -k2 .* sin.(p2) .* cos.(p2) -k2 .* sin.(p2).^2;
                0 0 -k2 .* cos.(p2).^2 -k2 .* sin.(p2) .* cos.(p2) k2 .* cos.(p2).^2 k2 .* sin.(p2) .* cos.(p2);
                0 0 -k2 .* sin.(p2) .* cos.(p2) -k2 .* sin.(p2).^2 k2 .* sin.(p2) .* cos.(p2) k2 .* sin.(p2).^2]) .* unit(k_c1)


# define External Forces vector
F(fx, fy) = [   0;
                0;
                fx;
                fy;
                0;
                0]

# define general Displacement vector
u(u3, u4)   =  [0;
                0;
                u3;
                u4;
                0;
                0]


function solve(ϕ1, ϕ2, Fx, Fy, k1, k2, A1, A2)
    # 6x6 stiffness matrix according to sheet
    K_mat = K(ϕ1, ϕ2, k1, k2)
    # 6x1 displacement vector
    F_vec = F(Fx, Fy)
    # solve stiffness equation for two free displacements
    u_free = inv(K_mat[3:4, 3:4]) * F_vec[3:4]
    # calculate reaction Forces
    F_r = K_mat * u(u_free[1], u_free[2])
    # rotation matrix
    T(α) = [    cos(α) -sin(α);
                sin(α) cos(α)]
    
    # stresses in element 1
    u1 = T((-ϕ1)) * u_free
    σ1 = u1[1] * k1 / A1

    # stresses in element 2
    u2 = T((ϕ2)) * u_free
    σ2 = u2[1] * k2 / A2

    σ_r = [σ1, σ2]

    return (u_free, F_r, σ_r)
end

function iterate_for_stress(ϕ1, ϕ2, Fx, Fy, E, ρ, RF, Rp02)
    A_iter = 100
    k_iter = E * A_iter / L_1
    σ_req = Rp02 / RF
    solved_iter = solve(ϕ1, ϕ2, Fx, Fy, k_iter, k_iter, A_iter, A_iter)
    σ_cache = maximum(solved_iter[3])
    residual = σ_req - σ_cache



    while abs(residual) > 0.001
        A_iter = A_iter * (1 - residual / σ_req)
        k_iter = E * A_iter / L_1
        solved_iter = solve(ϕ1, ϕ2, Fx, Fy, k_iter, k_iter, A_iter, A_iter)
        σ_cache = maximum(solved_iter[3])
        residual = σ_req - σ_cache
    end
    
    # results of optimal area and diameter assuming a circular rod
    A_opt = A_iter
    D_opt = sqrt(4A_opt/π)
    m_opt = 2 * (D_opt/1000)^2 * π/4 * L_1/1000 * ρ
    return (A_opt, D_opt, m_opt)
end

function iterate_for_deformation(ϕ1, ϕ2, Fx, Fy, E, ρ, u_max)
    A_iter = 100
    k_iter = E * A_iter / L_1
    u_req = u_max
    solved_iter = solve(ϕ1, ϕ2, Fx, Fy, k_iter, k_iter, A_iter, A_iter)
    u_cache = sqrt(solved_iter[1][1]^2 + solved_iter[1][2]^2)
    residual = u_req - u_cache



    while abs(residual) > 0.001
        A_iter = A_iter * (1 - residual / u_req)
        k_iter = E * A_iter / L_1
        solved_iter = solve(ϕ1, ϕ2, Fx, Fy, k_iter, k_iter, A_iter, A_iter)
        u_cache = sqrt(solved_iter[1][1]^2 + solved_iter[1][2]^2)
        residual = u_req - u_cache
    end
    
    # results of optimal area and diameter assuming a circular rod
    A_opt = A_iter
    D_opt = sqrt(4A_opt/π)
    m_opt = 2 * (D_opt/1000)^2 * π/4 * L_1/1000 * ρ
    return (A_opt, D_opt, m_opt)
end


function iterate(ϕ1, ϕ2, Fx, Fy, E, ρ, RF_stress, Rp02, u_max)
    result_stress = iterate_for_stress(ϕ1, ϕ2, Fx, Fy, E, ρ, RF_stress, Rp02)
    result_deformation = iterate_for_deformation(ϕ1, ϕ2, Fx, Fy, E, ρ, u_max)
    if result_deformation[3] > result_stress[3]
        println("Deformation is dimensioning")
        return result_deformation
    else
        println("Stress is dimensioning")
        return result_stress
    end
end

## Design parameters
RF_min = 3
u_max = 3

## Case 1
solved = solve(ϕ_c1,ϕ_c1, Fx_c1, Fy_c1, k_c1, k_c1, A_c1, A_c1)

println("CASE 1: B747")
println("Solved Case 1 (Area [mm²], Diameter [mm], Mass [kg]):")
println(iterate(ϕ_c1,ϕ_c1, Fx_c1, Fy_c1, E_ti, ρ_ti, RF_min, R_p02_ti, u_max))
println("----")

## Case 2
solved = solve(ϕ_c1,ϕ_c1, Fx_c2, Fy_c2, k_c2, k_c2, A_c2, A_c2)

println("CASE 2: Tank")
println("Solved Area Case 2 (Area [mm²], Diameter [mm], Mass [kg]):")
println(iterate(ϕ_c2, ϕ_c2, Fx_c2, Fy_c2, E_ti, ρ_ti, RF_min, R_p02_ti, u_max))
println("----")

println("Because Case 2 is dimensioning, the following calculations will only be performed for this case")
println("")
print("Ti6Al4V: ")
ret_ti = (iterate(ϕ_c2, ϕ_c2, Fx_c2, Fy_c2, E_ti, ρ_ti, RF_min, R_p02_ti, u_max))

print("Alu 7075: ")
ret_al = (iterate(ϕ_c2, ϕ_c2, Fx_c2, Fy_c2, E_al, ρ_al, RF_min, R_p02_al, u_max))

print("Stainless Steel: ")
ret_ss = (iterate(ϕ_c2, ϕ_c2, Fx_c2, Fy_c2, E_ss, ρ_ss, RF_min, R_p02_ss, u_max))

print("CFRP: ")
ret_cf = (iterate(ϕ_c2, ϕ_c2, Fx_c2, Fy_c2, E_cf, ρ_cf, RF_min, R_p02_cf, u_max))

names = ["Ti6Al4V"; "AL-7075"; "Stainless Steel"; "CFRP"]
areas = [ret_ti[1]; ret_al[1]; ret_ss[1]; ret_cf[1]]
diams = [ret_ti[2]; ret_al[2]; ret_ss[2]; ret_cf[2]]
masses = [ret_ti[3]; ret_al[3]; ret_ss[3]; ret_cf[3]]

result = DataFrame()
result.material = names
result."Area [mm]" = areas
result."Diameter [mm]" = diams
result."Mass [kg]" = masses
display(result)

# A = result."Area [mm]"[1]
# k = E_ti * A / L_2
# println(solve(ϕ_c2, ϕ_c2, Fx_c2, Fy_c2, k, k, A, A))
