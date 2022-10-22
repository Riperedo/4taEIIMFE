"""
calcula las raices de la ecuación de van der Waals
"""
function roots_vdW(P, P_0)
    N = length(P)
    index = []
    cond = P .> P_0
    for i in 2:N
        if (cond[i-1] == 1) & (cond[i] == 0) append!(index, i) end
    end
    if length(index) == 1 append!(index, length(cond)) end
    for i in index[1]:index[2]
        if (cond[i-1] == 0) & (cond[i] == 1) append!(index, i) end
    end
    return sort(index)
end

"""
encuentra el máximo y mínimo de la ecuación de van der Waals
"""
function max_vdW(P)
    cond = diff(P) .< 0
    index = []
    N = length(P)
    for i in 2:N-1
        if (cond[i-1] == 1) & (cond[i] == 0) 
            append!(index, i) 
            break
        end
    end
    for i in (index[1]+1):N-1
       if (cond[i-1] == 0) & (cond[i] == 1) 
            append!(index, i) 
            break
        end 
    end
    return index
end

"""
Criterio de Maxwell para áreas iguales
"""
function areas_iguales(V, P)
    function integral(x, f)
        dx = diff(x)
        return sum(f[1:end-1].*dx)
    end
    index_m = max_vdW(P)
    V_min = V[index_m[1]]
    V_max = V[index_m[2]]
    P_min = P[index_m[1]]
    P_max = P[index_m[2]]
    function Δv(P_0)
        index_r = roots_vdW(P, P_0)
        V1 = integral(V[index_r[1]:index_r[2]], (P.-P_0)[index_r[1]:index_r[2]])
        V2 = integral(V[index_r[2]:index_r[3]], (P.-P_0)[index_r[2]:index_r[3]])
        return abs(V1 + V2)
    end
    p = collect(P_min:0.001:P_max)
    dV = Δv.(p)
    minimo = findmin(dV)[2]
    return p[minimo]
end

"""
P_vdW(V, T, a, b)
`V::Real` densidad reducida
`T::Real` Temperatura reducia
Presión de la ecuación de van der Waals en escala reducuda:
* P_0 = a/27b^2
* V_0 = 3nb
* T_0 = 8a/27Rb
"""
P_vdW(V::Real, T::Real) = 8*T/(3*V-1) - 3/(V^2)


#V = collect(0.4:0.01:4)
#P = P_vdW.(V, 0.85)
#areas_iguales(V,P)
#index = roots_vdW(P, 0.5)
#println(index)
#max_vdW(P)