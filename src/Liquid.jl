#####################
#	Liquid Library	#
#####################
abstract type object end
power(x, n) = x^n

# Perfil de distribución de densidades
mutable struct Liquid <: object
	#atributos
	ρ_total :: Real
	ϕ_total :: Real
	ρ :: Array
	σ :: Array
	σ₂ :: Array
	μ :: Array
	ϕ :: Array
	ϵ :: Array
	T :: Real
	z :: Real
	soft :: Bool
	fase :: Any
	species :: Integer
	monodisperse :: Bool
	#add is_HS :: Bool
	#m'etodos
	setρ :: Function
	setϕ :: Function
	setDistribution :: Function
	setT :: Function
	setz :: Function
	setϵ :: Function
	setσ₂ :: Function
	setμ :: Function
	saveinfo :: Function
	setfase :: Function
	Soft :: Function
	function Liquid()
		this = new()
		this.ρ_total = 0.1
		this.ρ = collect(this.ρ_total:this.ρ_total) # por defecto es una sola densidad
		this.species = length(this.ρ)
		this.σ = collect(1.0:1.0)
		this.σ₂ = collect(1.0:1.0)
		this.μ = collect(0.0:0.0)
		this.ϵ = collect(1.0:1.0)
		this.ϕ = [π*this.σ[1]*this.ρ[1]/6.0]
		this.T = 1.0
		this.z = 0.0
		this.soft = false
		this.monodisperse = true
		this.setϕ = function (phi = false)
			if phi
				ρ = zeros(this.species)
				for i in 1:this.species
					ρ[i] = this.ϕ[i]*6.0/(π*power(this.σ[i],3))
				end
				this.ρ = ρ
				this.ρ_total = sum(this.ρ)
			else
				ϕ = zeros(this.species)
				for i in 1:this.species
					ϕ[i] = π*power(this.σ[i],3)*this.ρ[i]/6.0
				end
				this.ϕ = ϕ
				this.ϕ_total = sum(this.ϕ)
				@assert this.ϕ_total <= 1.0 "Inconsistencia entre la relación de radios y la distribucion de densidad. Obtienes una fracción de volumen ϕ = "*string(this.ϕ)
				#if this.ϕ >= 0.582 
				#	println("Warning: la fracción de volumen es alta ϕ = "*string(this.ϕ))
				#end
			end
		end
		this.setρ = function (rho :: Real)
			this.ρ = rho
		end
		this.setDistribution = function (rho :: Array, sigmas :: Array; phi = false)
			@assert length(rho) == length(sigmas) "La distribución de σs debe ser de la misma longitud que distribution de densidades. length(σ) = "*string(length(σ))*", length(ρ) = "*length(this.ρ)
			this.σ = sigmas
			this.σ₂ = sigmas
			if phi
				this.ϕ = rho
				this.ϕ_total = sum(rho)
			else
				this.ρ = rho
				this.ρ_total = sum(rho)
			end
			this.μ = this.μ[1]*ones(size(sigmas))
			this.ϵ = this.ϵ[1]*ones(size(sigmas))
			# actualiza ϕ
			this.setϕ(phi)
			#actualiza ρ
			this.species = length(this.ρ)
			if sigmas[1]*ones(size(sigmas)) == sigmas
				this.monodisperse = true
			else
				this.monodisperse = false
			end
		end
		this.setT = function (temperatura :: Real)
			this.T = temperatura
		end
		this.Soft = function ()
			this.soft = !this.soft
		end
		this.saveinfo = function(texto = "")
			info = "Densidad de partículas: "*string(this.ρ_total)*"\n"
			info *= "Fracción de volumen: "*string(sum(this.ϕ))*"\n"
			info *= "Temperatura : "*string(this.T)*"\n"
			info *= "Apantallamiento : "*string(this.z)*"\n"
			info *= "Número de especies: "*string(this.species)*"\n"
			info *= "Distribución de densidades: "*string(this.ρ)*"\n"
			info *= "Distribución de fracciones de volumen: "*string(this.ϕ)*"\n"
			info *= "Distribución de diámetros: "*string(this.σ)*"\n"
			info *= "Distribución de diámetros auxiliares: "*string(this.σ₂)*"\n"
			info *= "Distribución de energias auxiliares: "*string(this.ϵ)*"\n"
			info *= texto
			open("Liquid_info.dat", "w") do f
				write(f, info)
			end
		end
		this.setϵ = function(epsilon :: Array{Float64})
			@assert length(epsilon) == this.species
			this.ϵ = epsilon
		end
		this.setσ₂ = function(sigma :: Array{Float64})
			@assert length(sigma) == this.species
			this.σ₂ = sigma
		end
		this.setμ = function(mu :: Array{Float64})
			@assert length(mu) == this.species
			this.μ = mu
		end
		this.setz = function(zeta :: Float64)
			this.z = zeta
		end
		this.fase = Dump
		this.setfase = function(system)
			this.fase = system
		end
		return this
	end
end

#################################
#	algunas distribuciones 		#
#################################

function uniform(x :: Real, N :: Integer)
	@assert N > 0 "N no es mayor que cero"
	Δx = x/N
	dist = zeros(N)
	for i in 1:N
		dist[i] = Δx
	end
	return dist*x/sum(dist)
end

function gaussian(xmin :: Real, xmax :: Real, N :: Integer, μ :: Real, σ :: Real)
	@assert xmin > 0.0 "xmin es negativo"
	@assert xmax > 0.0 "xmax es negativo"
	@assert xmax > xmin "xmax es menor que xmin"
	@assert xmax > μ > xmin "la media no está en el rango xmin = "*string(xmin)*" y xmax = "*string(xmax)
	@assert N > 0 "N no es mayor que cero"
	Δx = (xmax - xmin)/(N-1)
	x = collect(xmin:Δx:xmax)
	dist = zeros(N)
	for i in 1:N
		dist[i] = exp(-(x[i] - μ)^2/(2.0*σ^2))*Δρ
	end
	return dist*μ/sum(dist)
end

function monodisperso(x :: Real, N :: Integer)
	dist = zeros(N)
	dist[1] = x
	return dist
end

function porcentual(x :: Real, A :: Array)
	dist = A
	return x*dist/sum(A)
end

#############
#	fases 	#
#############

struct Zero
	function Zero()
		this = new()
		return this
	end
end

struct Crystal
	function Crystal()
		this = new()
		return this
	end
end

struct Glass
	function Glass()
		this = new()
		return this
	end
end

struct Mixed
	function Mixed()
		this = new()
		return this
	end
end

struct Fluid
	function Fluid()
		this = new()
		return this
	end
end

struct Dump
	function Dump()
		this = new()
		return this
	end
end

#############################
#	fases in dipolar system	#
#############################

struct TotalArrest
	function TotalArrest()
		this = new()
		return this
	end
end

struct TraslationalArrest
	function TraslationalArrest()
		this = new()
		return this
	end
end

struct RotationalArrest
	function RotationalArrest()
		this = new()
		return this
	end
end

struct Basura
	function Basura()
		this = new()
		return this
	end
end
