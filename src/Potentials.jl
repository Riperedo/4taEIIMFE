abstract type object end

#########################################
# 		Potential definition			#
#########################################

mutable struct Potential <: object
	name :: String
	U :: Function
	info :: String
	isHS :: Bool
	spherical :: Bool
	FT :: Function
	# methods
	setName :: Function
	setPotential :: Function
	setInfo :: Function
	getInfo :: Function
	setHS :: Function
	setFT :: Function
	NonSpherical :: Function
	function FT_HS(k :: Real, L :: Liquid) end
	#constructor
	function Potential()
		this = new()
		this.name = ""
		this.setName =	
			function (Name :: String)
				this.name = Name
			end
		this.U = function() end
		this.setPotential = 
			function (u :: Function)
				this.U = u
			end
		this.info = ""
		this.setInfo = 
			function (prop :: String)
				this.info = prop
			end
		this.getInfo = 
			function ()
				if this.info == ""
					println("No info.")
				else
					println(this.info)
				end
			end
		this.isHS = false
		this.setHS = 
			function ()
				this.isHS = true
			end
		this.spherical = true
		this.NonSpherical =
			function ()
				this.spherical = false
			end
		this.FT = FT_HS
		this.setFT =
			function(F :: Function)
				this.FT = F
			end
		return this
	end
end

#########################################
#		Particular potentials 			#
#########################################


#########################################################
#	Fourier Transform of some potentials of the form	#
#			 | 0 	r < σ 								#
#	βU(r) = -|											#
#			 | v(r)	r > σ 								#
#########################################################

#=		Hard Sphere potential 		=#
function FT_HS(k :: Real, L :: Liquid) return 0.0 end
HS = Potential()
HS.setName("HardSphere")
function U_HS(r :: Real, L :: Liquid)
	u = zeros(L.species, L.species)
	for α in 1:L.species
		for β in 1:L.species
			σ = 0.58*(L.σ[α]+L.σ[β])
			if r < σ
				u[α,β] = INFTY
			end
		end
	end
	return u
end
HS.setPotential(U_HS)
HS.setInfo("This is the Hard sphere potential field which returns zero if the argument r is bigger than it's diameter σ or infinity in other case. This function have one argument r and one default value σ = 1.0")
HS.setHS()
HS.setFT(FT_HS)

#=		WCA potential 		=#

function FT_WCA(k :: Float64, L :: Liquid)
	u = zeros(L.species, L.species)
	z = L.z
	for α in 1:L.species
		for β in 1:L.species
			σ = 0.5*(L.σ[α]+L.σ[β])
			ϵ = 0.5*(L.ϵ[α]+L.ϵ[β])
			u₀ = ϵ*4.0*π*σ^3
			q = k*σ
			I₁(x) = -(sinint(q*x)*(q*x)^10+((q*x)^8-6.0*(q*x)^6+120.0*(q*x)^4-5040.0*(q*x)^2+362880.0)*sin(q*x) + (q*x)*((q*x)^8-2.0*(q*x)^6+24.0*(q*x)^4-720.0*(q*x)^2+40320.0)*cos(q*x))/(3628800*q*x^10)
			I₂(x) = (sinint(k*x)*(k*x)^4+((k*x)^2-6.0)*sin(k*x)+k*x*((k*x)^2-2.0)*cos(k*x))/(24.0*k*x^4)
			I₃(x) = (sin(q*x)-q*x*cos(q*x))/(q^3)
			b = 2.0^(1.0/6.0)
			u[α,β] = u₀*4.0*(I₁(b)-I₁(1.0)-I₂(b)+I₂(1.0))+u*(I₃(b)-I₃(1.0))
		end
	end
	return u
end

WCA = Potential()
WCA.setName("WCA")
function U_WCA(r :: Real, L :: Liquid)
	u = zeros(L.species, L.species)
	for α in 1:L.species
		for β in 1:L.species
			σ = 0.5*(L.σ[α]+L.σ[β])
			ϵ = 0.5*(L.ϵ[α]+L.ϵ[β])
			x6 = power(r/σ, 6)
			x12 = x6*x6
			xc = power(2.0, 1.0/6.0)
			if r < X_MIN/σ
				u[α,β] = INFTY
			else
				if r/σ < xc
					u[α,β] = 4.0*ϵ*(1.0/x12 - 1.0/x6) + ϵ
				else
					u[α,β] = 0.0
				end
			end
		end
	end
	return u	
end
WCA.setPotential(U_WCA)
WCA.setInfo("This potential field is called Weeks-Chandler-Anderson(WCA). Is composed by the repulsive part of Lennard-Jones (12-6) potential shifted up by it's characteristic energy unit ϵ. After it's minimum value at r = 2^{1/6}σ, where σ is its characteristic scale length, this field have value equals to zero. This field requires one argument r and two default arguments: ϵ = 1.0 and σ =1.0.")
WCA.setFT(FT_WCA)

#=		Soft Sphere potential 		=#
SS = Potential()
SS.setName("SoftSphere")
function U_SS(r :: Real, L :: Liquid; n = 12 :: Integer)
	u = zeros(L.species, L.species)
	for α in 1:L.species
		for β in 1:L.species
			σ = 0.5*(L.σ[α]+L.σ[β])
			ϵ = 0.5*(L.ϵ[α]+L.ϵ[β])
			x = r/σ
			xⁿ = power(x, n)
			if r/σ < X_MIN
				u[α,β] = INFTY
			else
				u[α,β] = ϵ/xⁿ
			end
		end
	end
	return u
end
SS.setPotential(U_SS)
SS.setInfo("The soft sphere potential field return simply the expression ϵσⁿ/rⁿ where ϵ is the energy unit and σ is the length unit. This field requires one argument r and have three default values: n = 12, ϵ = 1.0, σ = 1.0")

#=		Lennard-Jones potential 		=#

function FT_LJ(k :: Real, L :: Liquid)
	u = zeros(L.species, L.species)
	z = L.z
	for α in 1:L.species
		for β in 1:L.species
			σ = 0.5*(L.σ[α]+L.σ[β])
			ϵ = 0.5*(L.ϵ[α]+L.ϵ[β])
			u₀ = ϵ*4.0*π*σ^3
			q = k*σ
			I₁(x) = -(sinint(q*x)*(q*x)^10+((q*x)^8-6.0*(q*x)^6+120.0*(q*x)^4-5040.0*(q*x)^2+362880.0)*sin(q*x) + (q*x)*((q*x)^8-2.0*(q*x)^6+24.0*(q*x)^4-720.0*(q*x)^2+40320.0)*cos(q*x))/(3628800*q*x^10)
			I₂(x) = (sinint(q*x)*(q*x)^4+((q*x)^2-6.0)*sin(q*x)+q*x*((q*x)^2-2.0)*cos(q*x))/(24.0*q*x^4)
			u[α,β] = u₀*4.0*(-π*k^9/(7257600.0)-I₁(1.0)-π*k^3/48.0+I₂(1.0))
		end
	end
	return u
end

LJ = Potential()
LJ.setName("Lennard-Jones")
function U_LJ(r :: Real, L :: Liquid)
	u = zeros(L.species, L.species)
	for α in 1:L.species
		for β in 1:L.species
			σ = 0.5*(L.σ[α]+L.σ[β])
			ϵ = 0.5*(L.ϵ[α]+L.ϵ[β])	
			x = r/σ
			xᵐ = power(x, m)
			xⁿ = power(x, n)
			if r/σ < X_MIN
				u[α,β] = INFTY
			else
				u[α,β] = 4.0*ϵ*(1.0/xᵐ - 1.0/xⁿ)
			end
		end
	end
	return u
end
LJ.setPotential(U_LJ)
LJ.setInfo("The Lennard-Jones Potential is composed by the sum of a repulsive part ϵ/xᵐ and an attractive part - ϵ/xⁿ where x = r*σ. This potential have one argument r and four default arguments: m = 12, n = 6, ϵ = 1.0 and σ = 1.0.")
LJ.setFT(FT_LJ)

#= Square Well potential 	=#

function FT_Square_Well(q :: Float64, L :: Liquid)
	u = zeros(L.species, L.species)
	for α in 1:L.species
		for β in 1:L.species
			σ₁ = 0.5*(L.σ[α]+L.σ[β])
			σ₂ = 0.5*(L.σ₂[α]+L.σ₂[β])
			ϵ = 0.5*(L.ϵ[α]+L.ϵ[β])
			k = q*σ
			if k < 0.5
				u[α,β] = 4.0*π*ϵ*(σ₂^3-σ₁^3)/3.0 - 4.0*π*k*k*(σ₂^5-σ₁^5)/30.0
			else
				u[α,β] = 4.0*π*ϵ*(sin(σ₂*k) - σ₂*k*cos(σ₂*k) - sin(σ₁*k) + σ₁*k*cos(σ₁*k))/(k^3)
			end
		end
	end
	return u
end

SW = Potential()
SW.setName("Square Well")
function U_SW(r :: Real, L :: Liquid)
	u = zeros(L.species, L.species)
	for α in 1:L.species
		for β in 1:L.species
			σ₀ = 0.5*(L.σ[α]+L.σ[β])
			σ₁ = 0.5*(L.σ₂[α]+L.σ₂[β])
			ϵ = 0.5*(L.ϵ[α]+L.ϵ[β])	
			if r < σ₀
				u[α,β] = INFTY
			elseif r < σ₁
				u[α,β] = -ϵ
			else
				u[α,β] = 0.0
			end
		end
	end
	return u
end
SW.setPotential(U_SW)
SW.setFT(FT_Square_Well)

#= Square Shoulder potential 	=#

FT_Square_Shoulder(q :: Float64, L :: Liquid) = -FT_Square_Well(q, L)

SSh = Potential()
SSh.setName("Square Shoulder")
function U_SSh(r :: Real, L :: Liquid)
	u = zeros(L.species, L.species)
	for α in 1:L.species
		for β in 1:L.species
			σ₀ = 0.5*(L.σ[α]+L.σ[β])
			σ₁ = 0.5*(L.σ₂[α]+L.σ₂[β])
			ϵ = 0.5*(L.ϵ[α]+L.ϵ[β])	
			if r < σ₀
				u[α,β] = INFTY
			elseif r < σ₁
				u[α,β] = ϵ
			else
				u[α,β] = 0.0
			end
		end
	end
	return u
end
SSh.setPotential(U_SSh)
SSh.setFT(FT_Square_Shoulder)

#=		Yukawa potential 		=#

function FT_Yukawa(q :: Float64, L :: Liquid)
	u = zeros(L.species, L.species)
	z = L.z
	for α in 1:L.species
		for β in 1:L.species
			σ = 0.5*(L.σ[α]+L.σ[β])
			ϵ = 0.5*(L.ϵ[α]+L.ϵ[β])
			x = q*σ
			u₀ = ϵ*4.0*π*σ^3
			if x < 0.1
				factor = (z + 1)/z^2 - (x^2*(z^3 + 3*z^2 + 6*z + 6))/(6*z^4) + (x^4*(z^5 + 5*z^4 + 20*z^3 + 60*z^2 + 120*z + 120))/(120*z^6)
				u[α,β] = -u₀*factor
			else
				u[α,β] = -u₀*(x*cos(x)+z*sin(x))/(x*(x*x+z*z))
			end
		end
	end
	return u
end

Yk = Potential()
Yk.setName("Yukawa")
function U_Yk(r :: Real, L :: Liquid)
	u = zeros(L.species, L.species)
	for α in 1:L.species
		for β in 1:L.species
			σ = 0.5*(L.σ[α]+L.σ[β])
			ϵ = 0.5*(L.ϵ[α]+L.ϵ[β])
			z = L.z
			x = r/σ
			#xᵐ = power(x, m)
			#xⁿ = power(x, n)
			if r/σ < 1.0
				u[α,β] = INFTY
			else
				u[α,β] = -ϵ*exp(-z*(x-1))/x
			end
		end
	end
	return u
end
Yk.setPotential(U_Yk)
#Yk.setInfo("The Lennard-Jones Potential is composed by the sum of a repulsive part ϵ/xᵐ and an attractive part - ϵ/xⁿ where x = r*σ. This potential have one argument r and four default arguments: m = 12, n = 6, ϵ = 1.0 and σ = 1.0.")
Yk.setFT(FT_Yukawa)

#		Dipole-Dipole potential (only de radial component)		#
function FT_dd(q :: Float64, L :: Liquid) 
	j1(k) = sin(k)/k^2 - cos(k)/k
	# ∫(1,∞) dx x² j₂(kx)/x³
	I₆(k) = k < 0.5 ? 1/3 - k^2/30 + k^4/840 : j1(k)/k
	u = zeros(L.species, L.species)
	for α in 1:L.species
		for β in 1:L.species
			σ = 0.5*(L.σ[α]+L.σ[β])
			ϵ = 0.5*(L.ϵ[α]+L.ϵ[β])
			λ = L.μ[α]*L.μ[β]/L.T
			u[α,β] = -4*π*λ*I₆(q*σ)*σ^3
		end
	end
	return 0.0
end
dd = Potential()
dd.setName("Dipole-Dipole")
dd.NonSpherical()
function U_dd(r :: Real, L :: Liquid)
	u = 0.0
	for α in 1:L.species
		for β in 1:L.species
			σ = 0.5*(L.σ[α]+L.σ[β])
			ϵ = 0.5*(L.ϵ[α]+L.ϵ[β])
			λ = L.μ[α]*L.μ[β]/L.T
			x = r/σ
			if x < 0.1
				u[α,β] = INFTY
			else
				u[α,β] = -λ/x^3
			end
		end
	end
	return u
end
dd.setPotential(U_dd)
dd.setFT(FT_dd)
dd.setInfo("This potential only includes de radial component -λ/r^3.")
#=
using DelimitedFiles
if !isdir("data\\potentials") #if there are not directory "blimp"
	mkdir("data\\potentials") #then make it
end

potentials = [HS, SS, WCA, LJ]
r = collect(X_MIN:0.01:4.0)
for p in potentials
	#p.setσ(2)
	p.setϵ(1.5)
	u = map(R -> p.U(R), r)
	save_data("data\\potentials\\"*p.name*".dat", [r u])
end
=#




#println(WCA.name)
#println(typeof(WCA))
#WCA.getInfo()
#println(WCA.U(2.0))
#println(HS.σ)
#println(HS.U(HS.σ-0.1))
#HS.setσ(1.5)
#println(HS.σ)
#println(HS.U(HS.σ-0.1))
