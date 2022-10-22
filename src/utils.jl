#########################################
#		 Some useful Definitions		#
#########################################
using DelimitedFiles
using LinearAlgebra
#using SpecialFunctions #https://juliamath.github.io/SpecialFunctions.jl/latest/

#################################################
# We define a abstract type "objet". It's 		#
# not necessary but it works as mental anchor	#
#################################################
abstract type object end

const INFTY = prevfloat(Inf) #9999999.9999
const X_MIN = 0.65
power(x, y) = x^(y)
ρ2ϕ(ρ) = π*ρ/6.0
ϕ2ρ(ϕ) = 6.0*ϕ/π

function bool2text(boolean)
	txt = ""
	if boolean
		txt = "True"
	else
		txt = "False"
	end
end

function num2text(x)
	txt = ""
	texto = string(x)
	if length(texto) > 8
		txt = texto[1:8]
	else
		txt = texto
	end
	if '.' in txt
		txt = replace(txt, "." => "p")
	end
	return txt
end

#save_data("datos.dat", [k, S])
function save_data(nombre, formato; flag = true)
	@assert typeof(nombre) == typeof("hola") "El primer argumento debe ser texto"
	# el formato debe estar entre parentesis cuadrados [x, y, z, ...]
	open(nombre, "w") do io
		writedlm(io, formato)
	end
	if flag	println("Data saved as ", nombre) end
end

function read_data(name)
	a = readdlm(name, Float64)
	return a[:, 1], a[:, 2]
end

function biseccion(condicion :: Function, A :: Real, T :: Real, tolerancia :: Real; flag = false)
	Aquiles = A
	Tortuga = T
	δ = 0.5
	paso(t, a, Δ) = Δ*(t - a)
	while Aquiles != Tortuga
		δA = paso(Tortuga, Aquiles, δ)
		if condicion(Aquiles+δA) # Si la condicion se cumple se acepta el paso
			Aquiles += δA
		else
			δ *= 0.5
		end

		if abs(δA) < tolerancia
			break
		end
		if flag println(Aquiles, " ", δA, " ", Aquiles+δA) end
	end
	return Aquiles
end

function sphericalbesselj(nu, x::T) where {T}
    besselj_nuhalf_x = besselj(nu + one(nu)/2, x)
    if abs(x) ≤ sqrt(eps(real(zero(besselj_nuhalf_x))))
        nu == 0 ? one(besselj_nuhalf_x) : zero(besselj_nuhalf_x)
    else
        √((float(T))(π)/2x) * besselj_nuhalf_x
    end
end

