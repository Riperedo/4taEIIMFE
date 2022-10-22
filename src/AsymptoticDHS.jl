include("StructureFactor.jl")

#############################################################
#		calcula y determina comportamiento asintótico		#
#############################################################

mutable struct TR <: object
	T :: Any
	R :: Any
	function TR(A, B)
		this = new()
		this.T = A
		this.R = B
		return this
	end
end

function γᴵ(L :: Liquid, G :: Grid, γ :: TR, F₀₀ :: Array{Float64}, F₁₀ :: Array{Float64}, F₁₁ :: Array{Float64})
	kc = 1.305*2*π
	integrandoT = zeros(G.n)
	integrandoR = zeros(G.n)
	for i in 1:G.n
		k = G.x[i]
		λ = TR(1.0/(1.0 + (k/kc)^2), 1.0)
		S₀₀ = F₀₀[i]
		S₁₀ = F₁₀[i]
		S₁₁ = F₁₁[i]
		fˢ₀₀ = λ.T*λ.R/(λ.T*λ.R + (k^2)*γ.T*λ.R + 0*(0+1)*γ.R*λ.T)
		fˢ₁₀ = λ.T*λ.R/(λ.T*λ.R + (k^2)*γ.T*λ.R + 1*(1+1)*γ.R*λ.T)
		fˢ₁₁ = λ.T*λ.R/(λ.T*λ.R + (k^2)*γ.T*λ.R + 1*(1+1)*γ.R*λ.T)
		f₀₀ = S₀₀*λ.T*λ.R/(S₀₀*λ.T*λ.R + (k^2)*γ.T*λ.R + 0*(0+1)*γ.R*λ.T)
		f₁₀ = S₁₀*λ.T*λ.R/(S₁₀*λ.T*λ.R + (k^2)*γ.T*λ.R + 1*(1+1)*γ.R*λ.T)
		f₁₁ = S₁₁*λ.T*λ.R/(S₁₁*λ.T*λ.R + (k^2)*γ.T*λ.R + 1*(1+1)*γ.R*λ.T)
		#  à la Gory-Ernesto
		integrandoT[i] = (k^4)*(((1-1/S₀₀)^2)*S₀₀*fˢ₀₀*f₀₀ + 3*((1-1/S₁₀)^2)*S₁₀*fˢ₁₀*f₁₀)
		integrandoR[i] = 12*(k^2)*((S₁₀-1)^2)*fˢ₁₁*f₁₁/S₁₁
	end
	
	IT = integral(G, integrandoT, rule = "trapezoidal")
	IR = integral(G, integrandoR, rule = "trapezoidal")
	gammaTI = max(IT/(6.0*π*π*L.ρ_total), 1e-12)
	gammaRI = max(IR/(16.0*π*π*L.ρ_total), 1e-12)
	#print("T = ", 1/gammaTI)
	#println(", R = ", 1/gammaRI)
	return TR(gammaTI, gammaRI)
end

function selector(L :: Liquid, γ :: TR, M :: TR, Error_Relativo :: Float64, infinito :: Float64)
	converge = TR(false, false)
	diverge = TR(false, false)
	converge.T = M.T < Error_Relativo
	converge.R = M.R < Error_Relativo
	diverge.T = γ.T > infinito
	diverge.R = γ.R > infinito
	# si algo divenge entonces no converge
	if diverge.T converge.T = false end
	if diverge.R converge.R = false end
	#print("diverge ", diverge)
	#println(", converge ", converge)
	if converge.T && converge.R # si la parte T y R convergen
		return TotalArrest()
	elseif converge.T && diverge.R # si la parte T converge y R diverge
		return TraslationalArrest()
	elseif diverge.T && converge.R # si la parte T diverge y R converge
		return RotationalArrest()
	elseif diverge.T && diverge.R # si la parte T diverge y R diverge
		return Fluid()
	else # ninguna de las anteriores
		return Dump()
	end
end

function inversa(M :: TR)
	inversa = true
	if M.T == 0.0 inversa = false end
	if M.R == 0.0 inversa = false end
	return inversa
end

function Asymptotic(L :: Liquid, G :: Grid, S00, S10, S11; flag= false)
	# iteraciones
	It_MAX = 1000
	decimo = div(It_MAX, 50)
	if flag
		println("|-------------------------|------------------------| <- 100%")
		print("|")
	end

	#outputs
	sistema = Dump()
	iteraciones = zeros(It_MAX+1)
	gammas = zeros(2, It_MAX+1)

	# seed
	γ = TR(1e-6, 1e-6)
	convergencia = TR(0, 0)
	gammas[1, 1] = γ.T
	gammas[2, 1] = γ.R

	#ciclo principal
	it = 1
	while true
		if it % decimo == 0 && flag
			print("#")
		end
		γ_new = γᴵ(L, G, γ, S00, S10, S11)
		if inversa(γ_new)
			γ_new = TR(1/γ_new.T, 1/γ_new.R)
			convergencia = TR(((γ.T - γ_new.T)/γ_new.T)^2, ((γ.R - γ_new.R)/γ_new.R)^2)
			γ = γ_new
		else
			println("γ no tiene inversa")
			println(γ_new)
			break
		end
		iteraciones[it+1] = it
	
		gammas[1, it+1] = γ.T
		gammas[2, it+1] = γ.R
		
		sistema = selector(L, γ, convergencia, 0.0001, 1e10)
		if typeof(sistema) == Fluid
			if flag print("Fluid") end
			break
		elseif typeof(sistema) == TotalArrest
			if flag print("Total Arrest") end
			break
		elseif typeof(sistema) == TraslationalArrest
			if flag print("Traslational Arrest") end
			break
		elseif typeof(sistema) == RotationalArrest
			if flag print("Rotational Arrest") end
			break
		elseif typeof(sistema) == Basura
			if flag print("Basura") end
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

function Asymptotic_structure(L, G, F₀₀, F₁₀, F₁₁)
	iteraciones, gammas, sistema = Asymptotic(L, G, F₀₀, F₁₀, F₁₁)
	γ = TR(gammas[1, Int(maximum(iteraciones))], gammas[2, Int(maximum(iteraciones))])
	kc = 1.305*2*π
	fˢ₀₀ = zeros(G.n)
	fˢ₁₀ = zeros(G.n)
	fˢ₁₁ = zeros(G.n)
	F₀₀ = zeros(G.n)
	F₁₀ = zeros(G.n)
	F₁₁ = zeros(G.n)
	integrandoT = zeros(G.n)
	integrandoR = zeros(G.n)
	for i in 1:G.n
		k = G.x[i]
		λ = TR(1/(1+(k/kc)^2), 1.0)
		S₀₀ = F₀₀[i]
		S₁₀ = F₁₀[i]
		S₁₁ = F₁₁[i]
		fˢ₀₀[i] = λ.T*λ.R/(λ.T*λ.R + (k^2)*γ.T*λ.R + 0*(0+1)*γ.R*λ.T)
		fˢ₁₀[i] = λ.T*λ.R/(λ.T*λ.R + (k^2)*γ.T*λ.R + 1*(1+1)*γ.R*λ.T)
		fˢ₁₁[i] = λ.T*λ.R/(λ.T*λ.R + (k^2)*γ.T*λ.R + 1*(1+1)*γ.R*λ.T)
		F₀₀[i] = S₀₀*S₀₀*λ.T*λ.R/(S₀₀*λ.T*λ.R + (k^2)*γ.T*λ.R + 0*(0+1)*γ.R*λ.T)
		F₁₀[i] = S₁₀*S₁₀*λ.T*λ.R/(S₁₀*λ.T*λ.R + (k^2)*γ.T*λ.R + 1*(1+1)*γ.R*λ.T)
		F₁₁[i] = S₁₁*S₁₁*λ.T*λ.R/(S₁₁*λ.T*λ.R + (k^2)*γ.T*λ.R + 1*(1+1)*γ.R*λ.T)
		integrandoT[i] = (k^4)*(((1-1/S₀₀)^2)*fˢ₀₀[i]*F₀₀[i] + 3*((1-1/S₁₀)^2)*fˢ₁₀[i]*F₁₀[i])
		integrandoR[i] = 12*(k^2)*((S₁₀-1)^2)*fˢ₁₁[i]*F₁₁[i]/(S₁₁^2)
	end
	ΔζT = integral(G, integrandoT, rule = "trapezoidal")/(6.0*π*π*L.ρ[1])
	ΔζR = integral(G, integrandoR, rule = "trapezoidal")/(16.0*π*π*L.ρ[1])
	return TR(ΔζT, ΔζR), fˢ₀₀, fˢ₁₀, fˢ₁₁, F₀₀, F₁₀, F₁₁
end
