include("utils.jl")
include("Grid.jl")
include("Liquid.jl")
include("Potentials.jl")

#############################################
# 	Structure Factor for a monodisperese 	#
#		dipole hard spheres colloid				#
#		Wertheim 1971							#
#############################################

#################################
#	Some derivates funtions		#
#################################

# Some special functions
js(ν, x) = 1.0#sphericalbesselj(ν, x)
Si(x) = sinint(x)
Ci(x) = cosint(x)


# Percus-Yevick solution's coeficients
α₁(ϕ) =  -(1.0 + 2.0*ϕ)^2/(1.0 - ϕ)^4
α₂(ϕ) = 6.0*ϕ*(1.0 + 0.5*ϕ)^2/(1.0 - ϕ)^4
α₃(ϕ) = 0.5*ϕ*α₁(ϕ)

# Some usefull integrals
# ∫(0,1) dx x² j₀(kx)
I₁(k) = k == 0.0 ? 1/3 : (sin(k) - k*cos(k))/k^3
# ∫(0,1) dx x² xj₀(kx)
I₂(k) = k == 0.0 ? 1/4 : (-(k^2 - 2.0)*cos(k) + 2.0*k*sin(k) - 2.0)/k^4
# ∫(0,1) dx x² x³j₀(kx)
I₃(k) = k == 0.0 ? 1/6 : (4.0*k*(k^2 - 6.0)*sin(k) - (k^4 - 12.0*k^2 + 24.0)*cos(k) + 24.0)/k^6
# ∫(0,1) dx x² xj₂(kx)
I₄(k) = k == 0.0 ? 0.0 : ((k^2 - 8.0)*cos(k) - 5.0*k*sin(k) + 8.0)/k^4
# ∫(0,1) dx x² x³j₂(kx)
I₅(k) = k == 0.0 ? 0.0 : ((48.0*k - 7.0*k^3)*sin(k) + (k^4 - 24.0*k^2 + 48.0)*cos(k) - 48.0)/k^6
# ∫(1,∞) dx x² j₂(kx)/x³
I₆(k) = k < 0.5 ? 1/3 - k^2/30 + k^4/840 : js(1,k)/k

# FT of Wertheim's direct correlation functions
c(ϕ, k) = α₁(ϕ)*I₁(k) + α₂(ϕ)*I₂(k) + α₃(ϕ)*I₃(k)
cΔ(ϕ, κ, k) = 2*κ*((α₁(2*κ*ϕ) - α₁(-κ*ϕ))*I₁(k) + (α₂(2*κ*ϕ) - α₂(-κ*ϕ))*I₂(k) + (α₃(2*κ*ϕ) - α₃(-κ*ϕ))*I₃(k))
cD(ϕ, κ, λ, k) = -0.25*κ*(2*α₂(2*κ*ϕ) + α₂(-κ*ϕ))*I₄(k) - 0.5*κ*(2*α₃(2.0*κ*ϕ) + α₃(-κ*ϕ))*I₅(k) - λ*I₆(k)

# Static structure factor
function SF_PY(ϕ :: Real, k :: Real)
	ρ = ϕ2ρ(ϕ)
	return 1.0/(1.0 - 4.0*π*ρ*c(ϕ, k))
end

function SF_Wertheim(ϕ :: Real, κ :: Real, λ :: Real, k :: Real)
	ρ = ϕ2ρ(ϕ)
	c¹¹₀ = 1.0 - 4.0*π*ρ*(cΔ(ϕ, κ, k) + 2*cD(ϕ, κ, λ, k))/3
	c¹¹₁ = 1.0 - 4.0*π*ρ*(cΔ(ϕ, κ, k) - cD(ϕ, κ, λ, k))/3
	return 1/c¹¹₀, 1/c¹¹₁
end

function SF(L ::  Liquid, G :: Grid, P :: Potential; VerletWeis = false)
	ϕ = L.ϕ[1]
	σ = L.σ[1]
	if P.spherical
		if VerletWeis & (P == HS)
			σ = (1 - L.ϕ[1]/16)^(1/3)
			ϕ = L.ϕ[1] - (L.ϕ[1]^2)/16
		end
		if L.soft
			T = L.T
			λ = T != 0.0 ? blip(1/T) : 1.0
			λ3 = λ^3
			σ = λ*((1 - λ3*L.ϕ[1]/16)^(1/3))
			ϕ = λ3*L.ϕ[1]*(1 - λ3*L.ϕ[1]/16)
		end
		S = zeros(G.n)
		ρ = L.ρ_total
		for i in 1:G.n
			k = G.x[i]
			Sᵛʷ = SF_PY(ϕ, k*σ)
			U = P.FT(k, L)[1,1]
			Sₛₛ⁻¹ = 1/Sᵛʷ + ρ*U # (1 - ρ(CHS - βU))
			S[i] = 1/Sₛₛ⁻¹
		end
		return S
	else
		if VerletWeis
			ϕ = L.ϕ[1] - (L.ϕ[1]^2)/16
		end
		λ = 1/L.T
		y = 8*λ*ϕ/3
		κ = ξ(y)/ϕ
		S00 = zeros(G.n)
		S10 = zeros(G.n)
		S11 = zeros(G.n)
		for i in 1:G.n
			k = G.x[i]
			S00[i] = SF_PY(ϕ, k)
			S10[i], S11[i] = SF_Wertheim(ϕ, κ, 1/L.T, k)
		end
		return S00, S10, S11
	end
end

#########################################
#	Structure Factor Blip function 	#
#########################################

function blip(ϵ :: Float64; ν = 6, flag = false)
	λ = 0.0
	dx = 1/1000
	for i in 1:1000
		x = i*dx
		λ -= x*x*exp(-ϵ*(1/x^(2*ν)-2/x^ν+1))
	end
	λ *= 3*dx
	λ += 1.0
	λ = λ^(1/3)
	if flag println("blip function computed ", λ) end
	return λ
end

#####################
#	κ parameter 	#
#####################

function ξ(y₀ :: Real)
	δξ = 0.0001
	ξf = 0.5
	q(ϕ) = -α₁(ϕ)
	function condicion(ξ)
		y(ξ) = (q(2*ξ)-q(-ξ))/3
		return y(ξ) < y₀
	end
	ξ = biseccion(condicion, 0.0, ξf, δξ) # Aquiles inicia en cero, la tortuga está en ξf
	return ξ
end


"""
Function: Auxiliary function that computes for the FT of the perturbation potential -β u_p(k)
System: Square Well
Inputs:
	-T* = Dimensionless temperature / Domain (0:infty)
	-lambda range of the well in terms of the diameter sigma (1:infty)
	-k: wave vector magnitude / type double / Domain [0:infty)
Outputs:
	-β u_p(k): FT of the perturbation well
"""
function sw_m_beta_uk(T :: Real, λ :: Real, k :: Real)
	k2 = k * k
	if k > 0.0750
		λk = λ * k 
		sink = sin(k)
		cosk = cos(k)
		sinλk = sin(λk)
		cosλk = cos(λk)
		c_aux = ((cosk - λ * cosλk) / k) + ((sinλk - sink) / k2) 
		c_aux = c_aux / k 	
	else 
		λ3 = λ ^ 3
		λ5 = λ3 * λ * λ
		c_aux = (1.0/3.0) * (λ3 - 1.0) - (1.0/30.0) * (λ5 - 1.0) * k2
	end
	return 4.0 * π * c_aux / T
end

"""
Function: FT Direct correlation function
System: Hard Sphere + Square Well
Approximation: Verlet-Weiss + Sharma-Sharma
Inputs:
	-phi: Volume fraction / type double / Domain (0:1)
	-T* = Dimensionless temperature / Domain (0:infty)
	-lambda range of the well in terms of the diameter sigma / Domain (1:infty)
	-k: wave vector magnitude / type double / Domain [0:infty)
Outputs:
	-c: FT direct correlation / type double / Range NA
"""
function ck_hssw_vwsh(ϕ :: Real, T :: Real, λ :: Real, k :: Real)
	chs = ck_hs_vw(ϕ, k)
	return chs + sw_m_beta_uk(T, λ, k)
end



"""
Function: Inverse of static structure factor
System: Hard Sphere
Approximation: Percus Yevick
Inputs:
	-phi: Volume fraction / type double / Domain (0:1)
	-k: wave vector magnitude / type double / Domain [0:infty)
Outputs:
	-is: Inverse of static structure factor / type double / Range (0:infty)
"""
function is_hssw_vwsh(ϕ :: Real, T :: Real, λ :: Real, k :: Real)
	ϕ_vw = ϕ*(1.0 - (ϕ / 16.0))
	chs = ck_hs_vw(ϕ, k)
	is = (ϕ_vw * chs) + (ϕ * sw_m_beta_uk(T, λ, k))
	return 1.0 - 6.0 * is / π
end

"""
Function: Static structure factor
System: Hard Sphere
Approximation: Percus Yevick
Inputs:
  -phi: Volume fraction / type double / Domain (0:1)
  -k: wave vector magnitude / type double / Domain [0:infty)
Outputs:
  -s: Static structure factor / type double / Range [0:infty)
"""
function sk_hssw_vwsh(ϕ :: Real, T :: Real, λ :: Real, k :: Real)
	return 1.0 / is_hssw_vwsh(ϕ, T, λ, k)
end

"""
Function: FT Direct correlation function
System: Hard Sphere
Approximation: Percus Yevick
Inputs:
	-ϕ: Volume fraction / type double / Domain (0:1)
	-k: wave vector magnitude / type double / Domain [0:infty)
Outputs:
	-c: FT direct correlation / type double / Range NA
"""
function ck_hs_py(ϕ :: Real, k :: Real)
	α1 = - (1.0 + 2.0 * ϕ)^2
	α2 = 6.0 * ϕ * ((1.0 + (ϕ / 2.0))^2)
	α3 = - 0.5 * ϕ * ((1.0 + (2.0 * ϕ))^2)
	aux_den = (1.0 - ϕ)
	aux_den = (aux_den^4)
	α1 = α1 / aux_den
	α2 = α2 / aux_den
	α3 = α3 / aux_den
	k² = k * k
	if k > 0.0750
		sink = sin(k)
		cosk = cos(k)
		k³ = k² * k
		k⁴ = k³ * k
		k⁶ = k⁴ * k²
		c = (α1 * (sink - k * cosk) / (k³)) + (α2 * (((2.0 * k) * sink) + ((- (k²) + 2.0) * cosk) - 2.0) / (k⁴)) + (α3 * (((4.0 * (k³) - 24.0 * k) * sink) + ((- (k⁴) + 12.0 * (k²) - 24.0) * cosk) + 24.0) / (k⁶))
	else
		c = α1 * ((1.0 / 3.0) - (k² / 30.0)) + α2 * (0.250 -	(k² / 36.0)) + α3 * ((1.0 / 6.0) - (k² / 48.0))
	end
	c *= 4.0 * π
	return c
end

"""
Function: Inverse of static structure factor
System: Hard Sphere
Approximation: Percus Yevick
Inputs:
	-ϕ: Volume fraction / type double / Domain (0:1)
	-k: wave vector magnitude / type double / Domain [0:infty)
Outputs:
	-is: Inverse of static structure factor / type double / Range (0:infty)
"""
function is_hs_py(ϕ :: Real, k :: Real)
	return 1 - 6.0 * ϕ * ck_hs_py(ϕ, k) / π
end

"""
Function: Static structure factor
System: Hard Sphere
Approximation: Percus Yevick
Inputs:
	-ϕ: Volume fraction / type double / Domain (0:1)
	-k: wave vector magnitude / type double / Domain [0:infty)
Outputs:
	-s: Static structure factor / type double / Range [0:infty)
"""
function sk_hs_py(ϕ :: Real, k :: Real)
	return 1.0 / is_hs_py(ϕ, k)
end

"""
Function: Inverse of static structure factor
System: Hard Sphere
Approximation: Percus Yevick + Verlet Weis
Ref: doi=10.1103/PhysRevA.5.939
Inputs:
	-ϕ: Volume fraction / type double / Domain (0:1)
	-k: wave vector magnitude / type double / Domain [0:infty)
Outputs:
	-is: Inverse of static structure factor / type double / Range (0:infty)
"""
function ck_hs_vw(ϕ :: Real, k :: Real)
	ϕ_vw = ϕ*(1.0 - (ϕ / 16.0)) # Density correction from Verlet-Weiss
	k_vw = k * ((ϕ_vw / ϕ)^(1.0 / 3.0)) # Wave vector correction from Verlet-Weiss
	return ck_hs_py(ϕ_vw,k_vw)
end

"""
Function: Inverse of static structure factor
System: Hard Sphere
Approximation: Percus Yevick + Verlet Weis
Ref: doi=10.1103/PhysRevA.5.939
Inputs:
  -phi: Volume fraction / type double / Domain (0:1)
  -k: wave vector magnitude / type double / Domain [0:infty)
Outputs:
  -is: Inverse of static structure factor / type double / Range (0:infty)
"""
function is_hs_vw(ϕ :: Real, k :: Real)
	ϕ_vw = ϕ*(1.0 - (ϕ / 16.0)) # Density correction from Verlet-Weiss
	k_vw = k * ((ϕ_vw / ϕ)^(1.0 / 3.0)) # Wave vector correction from Verlet-Weiss
	return is_hs_py(ϕ_vw, k_vw)
end

"""
Function: Static structure factor
System: Hard Sphere
Approximation: Percus Yevick + Verlet Weis
Ref: doi=10.1103/PhysRevE.87.052306
Inputs:
  -phi: Volume fraction / type double / Domain (0:1)
  -k: wave vector magnitude / type double / Domain [0:infty)
Outputs:
  -s: Static structure factor / type double / Range [0:infty)
"""
function sk_hs_vw(ϕ :: Real, k :: Real)
	return 1.0 / is_hs_vw(ϕ, k)
end

"""
Blip function
Ref: doi=10.1103/PhysRevE.87.052306
λ³(T, ν) = 1 - 3∫₀¹dx x²exp(-(1/T)(1/x^(2ν)-2/x^ν + 1))
"""
function blip(T :: Real; ν = 6)
	if T == 0 return 1.0
	else
		λ³ = 0.0
		dx = 1/1000
		for i in 1:1000
			x = i*dx
			λ -= x*x*exp(-(1/T)*(1/x^(2*ν)-2/x^ν+1))
		end
		λ³ *= 3*dx
		λ³ += 1.0
		return λ³^(1/3), λ³
	end
end


"""
Function: Inverse of static structure factor
System: Hard Sphere
Approximation: Percus Yevick + Verlet Weis + Blip function
Ref: doi=10.1103/PhysRevE.87.052306
Inputs:
	-ϕ: Volume fraction / type double / Domain (0:1)
	-T: Temperature / type double / Domain (0: infty)
	-k: wave vector magnitude / type double / Domain [0:infty)
Outputs:
	-is: Inverse of static structure factor / type double / Range (0:infty)
"""
function ck_hs_vw_blip(ϕ :: Real, T :: Real, k :: Real; ν = 6)
	λ, λ³ = blip(T, ν)
	return ck_hs_vw(λ³*ϕ, λ*k)
end

"""
Function: Inverse of static structure factor
System: Hard Sphere
Approximation: Percus Yevick + Verlet Weis + Blip function
Ref: doi=10.1103/PhysRevE.87.052306
Inputs:
  -phi: Volume fraction / type double / Domain (0:1)
  -k: wave vector magnitude / type double / Domain [0:infty)
Outputs:
  -is: Inverse of static structure factor / type double / Range (0:infty)
"""
function is_hs_vw_blip(ϕ :: Real, T :: Real, k :: Real; ν = 6)
	λ, λ³ = blip(T, ν = ν)
	return is_hs_vw(λ³*ϕ, λ*k)
end

"""
Function: Static structure factor
System: Hard Sphere
Approximation: Percus Yevick + Verlet Weis + Blip function
Ref: doi=10.1103/PhysRevE.87.052306
Inputs:
  -phi: Volume fraction / type double / Domain (0:1)
  -k: wave vector magnitude / type double / Domain [0:infty)
Outputs:
  -s: Static structure factor / type double / Range [0:infty)
"""
function sk_hs_vw_blip(ϕ :: Real, T :: Real, k :: Real; ν = 6)
	return 1.0 / is_hs_vw_blip(ϕ, T, k, ν = ν)
end

