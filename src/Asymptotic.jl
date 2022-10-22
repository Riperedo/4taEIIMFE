include("StructureFactor.jl")

#############################################################
#		calcula y determina comportamiento asintótico		#
#############################################################

function γᴵ(L :: Liquid, G :: Grid, γ :: Real, F₀ :: Array{Float64})
	kc = 1.305*2*π*L.σ[1]
	gammaI = 0.0
	integrando = zeros(G.n)
	for i in 1:G.n
		k = G.x[i]
		λ = 1.0/(1.0 + (k/kc)^2)
		S = F₀[i]
		fs = λ/(λ + (k*k)*γ)
		f = λ*S/(λ*S + (k*k)*γ)
		integrando[i] = (k^4)*(S-1)*f*fs*(S-1)/S
	end
	
	I = integral(G, integrando, rule = "trapezoidal")
	gammaI = max(I/(6.0*π*π*L.ρ_total), 1e-12)
	return gammaI
end

function selector(L :: Liquid, γ :: Real, M :: Real, Error_Relativo :: Float64, infinito :: Float64)
	Converge = abs(M) < Error_Relativo 
	Diverge = γ > infinito

	if Diverge
		return Fluid()
	elseif Converge
		return Glass()
	else
		return Dump()
	end
end

function Asymptotic(L :: Liquid, G :: Grid, S :: Array{Float64, 1}; flag= false)
	# iteraciones
	It_MAX = 1000
	decimo =div(It_MAX, 50)
	if flag
		println("|-------------------------|------------------------| <- 100%")
		print("|")
	end

	#outputs
	sistema = Dump()
	iteraciones = zeros(It_MAX+1)
	gammas = zeros(It_MAX+1)

	# seed
	γ = 0.00001
	gammas[1] = γ

	#ciclo principal
	it = 1
	while true
		if it % decimo == 0 && flag
		    print("#")
		end
		#sistema= Dump()
		γ_new = γᴵ(L, G, γ, S)
		if γ_new > 0.0
			γ_new = 1/γ_new
			convergencia = (γ - γ_new)/γ_new
		    γ = γ_new
		else
			println(γ_new)
			sistema = Zero()
			break
		end
		
		iteraciones[it+1] = it
	
		gammas[it+1] = γ
		
		sistema = selector(L, γ, convergencia, 0.0001, 1e10)
		if typeof(sistema) == Fluid
			if flag print("Fluid") end
			break
		elseif typeof(sistema) == Glass
			if flag print("Glass") end
			break
		end
		it += 1
		if it > It_MAX break end
	end
	if it < It_MAX && flag
		qwerty = div(It_MAX-it, decimo)
		for dot in 1:qwerty-4 print(" ") end
	end
	if flag println("| ¡Listo!") end
	return iteraciones, gammas, sistema
end

function Asymptotic_structure(L, G, S)
	iteraciones, gammas, sistema = Asymptotic(L, G, S)
	γ = gammas[Int(maximum(iteraciones))]
	kc = 1.305*2*π*L.σ[1]
	F∞ = zeros(G.n)
	Fs∞ = zeros(G.n)
	integrando = zeros(G.n)
	for i in 1:G.n
		k = G.x[i]
		λ = 1/(1+(k/kc)^2)
		F∞[i] = λ*S[i]*S[i]/(λ*S[i] + (k^2)*γ)
		Fs∞[i] = λ/(λ + (k^2)*γ)
		integrando[i] += (k^4)*((S[i]-1)^2)*F∞[i]*Fs∞[i]/(S[i]^2)
	end
	Δζ∞ = integral(G, integrando, rule = "trapezoidal")/(6.0*π*π*L.ρ[1])
	return Δζ∞, Fs∞, F∞
end
