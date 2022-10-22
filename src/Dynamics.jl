include("Asymptotic.jl")

#########################
#		Dynamics		#
#########################

function Δη_(G, F, S)
	k = G.x
	dk = diff(k)
	dS = diff(S)
	dSdk = dS./dk
	integrando = (k[2:end].^4).*(dSdk.^2).*(F[2:end].^2)./(S[2:end].^4)
	return sum(integrando.*dk)/(60*π*π)
end

function memoria(L, G, S, F, Fs)
	ONE = ones(G.n)
	k = G.x
	integrando = (k.^4).*((S-ONE).^2).*F.*Fs./(S.^2)
	dk = diff(k)
	#return integral(G, integrando; rule = "trapezoidal")/(6*π*π*L.ρ[1])
	return sum(integrando[2:end].*dk)/(6*π*π*L.ρ[1])
end

function decimation(G, S, T, Δζ, D, W, Fs, F, T_save, Δζ_save, D_save, W_save, Fs_save, F_save, Δη_save, index)
	N = length(T)
	for n in Int(N//2):N
		if n%2 == 1
			append!(T_save, T[n]/2)
			append!(Δζ_save, Δζ[n])
			append!(Fs_save, Fs[index, n])
			append!(F_save, F[index, n])
			append!(Δη_save, Δη_(G, F[:,n], S))
			append!(D_save, D[n])
			append!(W_save, W[n])
		end
	end

	for n in 1:Int(N//2)
		T[n] = T[2*n]
		Δζ[n] = Δζ[2*n]
		D[n] = D[2*n]
		W[n] = W[2*n]
		Fs[:,n] = Fs[:,2*n]
		F[:,n] = F[:,2*n]
	end
	for n in Int(N//2):N
		T[n] = 0.0
		Δζ[n] = 0.0
		D[n] = 0.0
		W[n] = 0.0
		Fs[:,n] .= 0.0
		F[:,n] .= 0.0
	end
	return T, Δζ, D, W, Fs, F, T_save, Δζ_save, D_save, W_save, Fs_save, F_save, Δη_save
end

function step_free(L, G, S, F, Fs, Δζ, T, D, W, dt, n)
	ONE = ones(G.n)
	k = G.x
	α = ONE + dt*(k.*k)./S
	F[:,n+1] = F[:,n]./α
	α = ONE + dt*(k.*k)
	Fs[:,n+1] = Fs[:,n]./α
	Δζ[n+1] = memoria(L, G, S, F[:, n+1], Fs[:, n+1])
	T[n+1] = dt*n
	D[n+1] = 1 - Δζ[n+1]*T[n+1]
	W[n+1] = D[n+1]*dt + W[n]
	return T, F, Fs, Δζ, D, W
end

function Σᵢ(Δζ, F, n)
	Total = zeros(length(F[:,1]))
	for i in 2:(n-1)
		Total = Total + (Δζ[n+1-i]-Δζ[n-i])*F[:, i]
	end
	return Total
end

function ΣD(Δζ, D, n)
	Total = 0.0
	for i in 1:n-1
		#Total = Total + (Δζ[n+1-i]-Δζ[n-i])*D[i] + (D[n+1-i]-D[n-i])*Δζ[i]
		Total = Total + D[i]*(Δζ[n-i+1] - Δζ[n-i])
	end
	return Total
end

function ΣΔr2(Δζ, W, n)
	msdtsum = 0.0
	for i in 2:n-1
		msdtsum = msdtsum + (Δζ[i] * (W[n+1-i] - W[n-i]) )
	end
	return msdtsum
end


function F_dump(L, G, S, F, Δζ, dt, n)
	"""
F_dump(k, n) = α(k)λ(k){Δζₙ₋₁F₁(k) - Σᵢ₌₂ⁿ⁻¹[Δζₙ₊₁₋ᵢ-Δζₙ₋ᵢ]Fᵢ(k)} + Δt⁻¹α(k)Fₙ₋₁(k)
	"""
	ONE = ones(G.n)
	k = G.x
	kc = 2*π*1.305
	λ = ONE./(ONE + (k.*k)/(kc^2))
	α = ONE./(ONE/dt + (k.*k)./S + λ*Δζ[1])
	return α.*(λ.*(Δζ[n-1]*F[:,1] - Σᵢ(Δζ, F, n)) + F[:,n-1]/dt)
end

function F_it(L, G, S, F, Δζ, dt, n)
	"""
F_it(k, n) = α(k)λ(k)Δζₙ[S(k)-F₁(k)]
	"""
	ONE = ones(G.n)
	k = G.x
	kc = 2*π*1.305
	λ = ONE./(ONE + (k.*k)/(kc^2))
	α = ONE./(ONE/dt + (k.*k)./S + λ*Δζ[1])
	return α.*λ.*(S[:]-F[:,1])*Δζ[n]
end

function step(L, G, S, F, Fs, Δζ, T, D, W, dt, n)
	"""
Fₙ(k) = α(k)λ(k){Δζₙ₋₁F₁(k) - Σᵢ₌₂ⁿ⁻¹[Δζₙ₊₁₋ᵢ-Δζₙ₋ᵢ]Fᵢ(k)} + Δt⁻¹α(k)Fₙ₋₁(k)
		+ α(k)λ(k)Δζₙ[S(k)-F₁(k)]
	 = F_dump(k, n) + F_it(k, n)
where
F_dump(k, n) = α(k)λ(k){Δζₙ₋₁F₁(k) - Σᵢ₌₂ⁿ⁻¹[Δζₙ₊₁₋ᵢ-Δζₙ₋ᵢ]Fᵢ(k)} + Δt⁻¹α(k)Fₙ₋₁(k)
F_it(k, n) = α(k)λ(k)Δζₙ[S(k)-F₁(k)]
and
α(k) = [Δτ⁻¹I + k²DS⁻¹(k) + λ(k)Δζ₁(k)]⁻¹
	"""
	s = ones(G.n)
	Dump = F_dump(L, G, S, F, Δζ, dt, n)
	dump = F_dump(L, G, s, Fs, Δζ, dt, n)
	Δζ[n+1] = Δζ[n]
	F[:,n+1] = Dump + F_it(L, G, S, F, Δζ, dt, n)
	Fs[:,n+1] = dump + F_it(L, G, s, Fs, Δζ, dt, n)
	Δζ[n+1] = memoria(L, G, S, F[:, n+1], Fs[:, n+1])
	while true
	#for t in 1:10
		F[:,n+1] = Dump + F_it(L, G, S, F, Δζ, dt, n)
		Fs[:,n+1] = dump + F_it(L, G, s, Fs, Δζ, dt, n)
		Δζ_new = memoria(L, G, S, F[:, n+1], Fs[:, n+1])
		if (Δζ[n+1]-Δζ_new)^2 < 1e-10 break end
		Δζ[n+1] = Δζ_new
	end
	T[n+1] = dt*n
	# Calculando difusion
	suma = ΣD(Δζ, D, n)
	#d =(1/(1/dt+Δζ[1]))*(-Δζ[Int((n+1 - (n+1)%2)//2)]*D[Int((n+1 - (n+1)%2)//2)] - (Δζ[n+1]-Δζ[n])*D[1] + (1/dt+Δζ[1])*D[n] - suma)
	D[n+1] = (-suma + D[n]/dt)/( (1/dt) + Δζ[1])
	# calculando MSD
	msdtsum = ΣΔr2(Δζ, W, n) - (Δζ[1] + (1/dt)) * W[n-1]
	#W[n+1] = D[n+1]*dt + W[n]
	W[n+1] = dt*(1 - msdtsum - (Δζ[n+1]*W[1]))/(1 + dt*Δζ[1])
	return T, F, Fs, Δζ, D, W
end


function dynamics(L, G, S, k_max, dt, nT, decimaciones)
	# grid temporal
	N = 2<<nT
	F = zeros(G.n, N)
	Fs = zeros(G.n, N)
	Δζ = zeros(N)
	T = zeros(N)
	D = zeros(N)
	W = zeros(N)

	# output
	index = G.find(k_max)
	T_save = []
	Δζ_save = []
	Fs_save = []
	F_save = []
	Δη_save = []
	D_save = []
	W_save = []
	#initial conditions
	ONE = ones(G.n)
	kc = 2*π*1.305
	F[:, 1] = S
	Fs[:, 1] = ONE
	Δζ[1] = memoria(L, G, S, F[:, 1], Fs[:, 1])
	# first steps
	# free diffusion
	for n in 1:N-1
		T, F, Fs, Δζ, D, W = step_free(L, G, S, F, Fs, Δζ, T, D, W, dt, n)
	end
	for d in 1:decimaciones
		# decimation
		T, Δζ, D, W, Fs, F, T_save, Δζ_save, D_save, W_save, Fs_save, F_save, Δη_save = decimation(G, S, T, Δζ, D, W, Fs, F, T_save, Δζ_save, D_save, W_save, Fs_save, F_save, Δη_save, index)
		dt *= 2
		for n in Int(N//2):N
			#T, F, Fs, Δζ, D, W = step_free(L, G, S, F, Fs, Δζ, T, D, W, dt, n-1)
			T, F, Fs, Δζ, D, W = step(L, G, S, F, Fs, Δζ, T, D, W, dt, n-1)
		end
	end

	# saving final steps
	for n in (Int(N//2) + 1):N
		append!(T_save, T[n]/2)
		append!(Δζ_save, Δζ[n])
		append!(D_save, D[n])
		append!(W_save, W[n])
		append!(Fs_save, Fs[index, n])
		append!(F_save, F[index, n])
		append!(Δη_save, Δη_(G, F[:,n], S))
	end

	return T_save, Fs_save, F_save, Δζ_save, Δη_save, D_save, W_save
end

