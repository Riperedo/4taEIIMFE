include("Dynamics.jl")


"""
HS
"""
StructureFactor_HS(ϕ::Real, k::Real) = SF_PY(ϕ, k)

"""
HSVW
"""
function StructureFactor_HSVW(ϕ::Real, k::Real)
	σᵛʷ = (1 - ϕ/16)^(1/3)
	ϕᵛʷ = ϕ - (ϕ^2)/16
	return StructureFactor_HS(ϕᵛʷ, k*σᵛʷ)
end
D₀ = 1.0

"""
WCA
"""
function StructureFactor_WCA(ϕ::Real, T:: Real, k::Real)
	λ = blip(T)
	λ³ = λ^3
	return StructureFactor_HSVW(ϕ*λ³, k*λ)
end

"""
Función de dispersión intermedia
F(k,τ)
"""
function ISF(S::Array, k::Array, τ::Real; D₀ = 1.0 ::Real)
	return S.*exp.(-D₀*τ*(k.^2)./S)
end

"""
Auto-función de dispersión intermedia
F(k,τ)
"""
function sISF(k::Array, τ::Real; D₀ = 1.0 ::Real)
	return exp.(-D₀*τ*(k.^2))
end

### SCGLE

"""
Regresa un objeto que maneja una colección de parámetros necesarios para evaluar al SCGLE.
"""
struct Input_SCGLE
	L :: Liquid
	G :: Grid
	S :: Array
end
# Constructores
vector_de_onda(I::Input_SCGLE) = I.G.x
estructura(I::Input_SCGLE) = I.S
estructura_DHS(I::Input_SCGLE) = I.S[1], I.S[2], I.S[3]
fraccion_de_volumen(I::Input_SCGLE) = I.L.ϕ[1]
temperatura(I::Input_SCGLE) = I.L.T
"""
Regresa un objeto que maneja una colección de parámetros necesarios para evaluar al sistema de esferas duras dentro de la SCGLE.
"""
function Input_HS(kₘᵢₙ::Real, kₘₐₓ::Real, N::Integer, ϕ::Real; VW = false::Bool)
	G = Grid(kₘᵢₙ, kₘₐₓ, N)
	L = Liquid()
	L.setDistribution([ϕ], [1.0], phi = true)
	S = SF(L, G, HS, VerletWeis = VW)
	return Input_SCGLE(L, G, S)
end

"""
Regresa un objeto que maneja una colección de parámetros necesarios para evaluar al sistema de esferas suaves dentro de la SCGLE.
"""
function Input_WCA(kₘᵢₙ::Real, kₘₐₓ::Real, N::Integer, ϕ::Real, Temp::Real)
	G = Grid(kₘᵢₙ, kₘₐₓ, N)
	L = Liquid()
	L.setDistribution([ϕ], [1.0], phi = true)
	L.Soft()
	L.T = Temp
	S = SF(L, G, HS)
	return Input_SCGLE(L, G, S)
end

"""
regresa la longitud de localización
"""
function longitud_de_localizacion(I::Input_SCGLE; flag = true::Bool)
	iteraciones, gammas, sistema = Asymptotic(I.L, I.G, I.S, flag = flag)
	return gammas[Integer(maximum(iteraciones))]
end

"""
Regresa la evolución dinámica de un sistema coloidal
"""
function SCGLE(I::Input_SCGLE, k_max; dt = 1e-10::Real, nT = 6::Integer, decimaciones=50::Integer, flag = false)
	if flag print("Calculando ...") end
	τ, Fs, F, Δζ, Δη, D, W = dynamics(I.L, I.G, I.S, k_max, dt, nT, decimaciones)
	if flag println(" ¡Listo!") end
	return τ, Fs, F, Δζ, Δη
end


"""
Regresa la estructura de esferas duras dipolares
"""
function Input_DHS(kₘᵢₙ::Real, kₘₐₓ::Real, N::Integer, ϕ::Real, Temp::Real)
	G = Grid(kₘᵢₙ, kₘₐₓ, N)
	L = Liquid()
	L.setDistribution([ϕ], [1.0], phi = true)
	L.T = Temp
	S⁰⁰, S¹⁰, S¹¹ = SF(L, G, dd, VerletWeis = true)
	return Input_SCGLE(L, G, [S⁰⁰, S¹⁰, S¹¹])
end

print("Hola")