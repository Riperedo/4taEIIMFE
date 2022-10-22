include("Dynamics.jl")
include("Viscoelasticity.jl")


η∞(ϕ) = (1 + 1.5*ϕ*(1 + ϕ - 0.189*ϕ*ϕ))/(1 - ϕ*(1 + ϕ - 0.189*ϕ*ϕ))

function dη(τ, Δη)
	dη = 0.0
	n = length(τ)
	for i in 2:n
		dlog_τ = log(τ[i]) - log(τ[i-1])
		dη += 0.5*(Δη[i]*τ[i]+Δη[i-1]*τ[i-1])*dlog_τ
	end
	return dη
end

function mobility(t, Δζ)
	n = length(t)
	bI = 0.0
	for i in 2:n
		dlog_t = log(t[i]) - log(t[i-1])
		bI += 0.5*(Δζ[i]*t[i]+Δζ[i-1]*t[i-1])*dlog_t
	end
	bI += 1.0
	return bI
end


function NESF(Li ::  Liquid, Lf ::  Liquid, G :: Grid, P :: Potential, u; VerletWeis = false)
	Si = SF(Li, G, P, VerletWeis = VerletWeis)
	Sf = SF(Lf, G, P, VerletWeis = VerletWeis)
	Su = zeros(G.n)
	for i in 1:G.n
		k = G.x[i]
		Su[i] = Si[i]*exp(-2*(k*k*u)/Sf[i]) + Sf[i]*(1-exp(-2*(k*k*u)/Sf[i]))
	end
	return Su
end
#=
function MSD(τ, Δζ, Δη, ϕ)

	W = zeros(length(τ))
	δζ = zeros(length(τ))
	D = zeros(length(τ))
	δη = zeros(length(τ))
	for ii in 1:length(τ)-1
		dt = τ[ii+1] - τ[ii]
		δζ[ii+1] = δζ[ii] + dt*0.5*(Δζ[ii] + Δζ[ii+1])
		δη[ii+1] = δη[ii] + dt*0.5*(Δη[ii] + Δη[ii+1])
		#W[ii+1] = W[ii] + dt*0.5*(1/(1+δζ[ii+1])+1/(1+δζ[ii]))
		W[ii+1] = W[ii] + dt*0.5*(1/(η_∞(ϕ)+δη[ii+1])+1/(η_∞(ϕ)+δη[ii]))
		D[ii+1] = 1/(η_∞(ϕ)+δη[ii+1])
	end
	return δζ, δη, W, D
end

function Mason(τ, Δr²)
	N = length(τ)
	z = []
	G = []
	G´ = []
	G´´ = []
	for i in 3:N
		s = 1/τ[i]
		α = (log(Δr²[i-1])-log(Δr²[i]))/(log(τ[i-1])-log(τ[i]))
		g = 6/(Δr²[i]*gamma(1+α))
		g´ = g*cos(π*α/2)
		g´´ = g*sin(π*α/2)
		append!(z, s)
		append!(G, g)
		append!(G´, g´)
		append!(G´´, g´´)
	end
	return z, G, G´, G´´
end
=#
function instantaneous_process(phi_i, Ti, phi_f, Tf, G; soft = false, k₀ = 7.2, VerletWeis = false, u_max = 10.0, u_min = 1e-6, factor  =2.0, mute= false, save_path = "")
	#save path
	#main = "monodisperse_hard_sphere"
	main = "HSSW_lambda1p5"
	if soft
		#main = "monodisperse_soft_sphere"
		main = "HSSW_lambda1p5"
	end
	folder = save_path*main
	if !isdir(folder) #if there are not directory "plots"
		mkdir(folder) #then make it
	end

	ϕi = phi_i	# phi_i nos va a servir después para la estructura
	ϕf = phi_f	# phi_i nos va a servir después para la estructura
	phi_save = ""
	if ϕi == ϕf
		phi_save *= "phi"*num2text(ϕi)
	else
		phi_save *= "phii"*num2text(ϕi)*"phif"*num2text(ϕf)
	end

	folder *= "/"*phi_save
	if !isdir(folder) #if there are not directory "plots"
		mkdir(folder) #then make it
	end

	T_save = ""
	if Ti == Tf
		T_save = "T"*num2text(Ti)
	else
		T_save = "Ti"*num2text(Ti)*"Tf"*num2text(Tf)
	end

	folder *= "/"*T_save
	if !isdir(folder) #if there are not directory "plots"
		mkdir(folder) #then make it
	end

	println("Start")

	tw = 0.0
	u = 1e-6
	bI = 1.0

	#some arrays to save data
	tw_save = []
	u_save = []
	bi_save = []
	b_save = []
	tau_save = []
	eta_save = []
	Δeta_save = []
	#arrays to save Viscosity TFL initial data

	Δt = 1e-6
	j = 0
	j_save = []

	# Initial structure
	Li = Liquid()
	ϕi = [phi_i]
	σi = [1.0]
	Li.setDistribution(ϕi, σi, phi = true)
	Li.T = Ti
	if soft
		Li.Soft()
	end

	Si = SF(Li, G, HS, VerletWeis = VerletWeis)

	τ, Fs, F, Δζ, Δη, D, W = dynamics(Li, G, Si, k₀, 1e-10, 5, 60)

	save_data(folder*"/inputI.dat", [G.x Si])
	save_data(folder*"/fsI.dat", [τ Fs F Δζ D W Δη])

	z, G_r, G´, G´´ = Mason(τ, W)
	save_data(folder*"/ViscoelasticityI.dat", [z G_r G´ G´´ z])

	δζ, δη, W, D = MSD(τ, Δζ, Δη, phi_i)
	save_data(folder*"/prop_transpI.dat", [τ δζ δη W D])

	z, G_r, G´, G´´ = Mason(τ, W)
	save_data(folder*"/Viscoelasticity0I.dat", [z G_r G´ G´´ z])

	# Final structure
	Lf = Liquid()
	ϕf = [phi_f]
	σf = [1.0]
	Lf.setDistribution(ϕf, σf, phi = true)
	Lf.T = Tf
	if soft
		Lf.Soft()
	end

	Sf = SF(Lf, G, HS, VerletWeis = VerletWeis)

	#=
	τ, Fs, F, Δζ, Δη, D, W = dynamics(Lf, G, Sf, k₀, 1e-10, 5, 60)

	save_data(folder*"/inputF.dat", [G.x Si])
	save_data(folder*"/fsF.dat", [τ Fs F Δζ D W Δη])
	z, G_r, G´, G´´ = Mason(τ, W)
	save_data(folder*"/ViscoelasticityF.dat", [z G_r G´ G´´ z])

	δζ, δη, W, D = MSD(τ, Δζ, Δη, phi_f)
	save_data(folder*"/prop_transpF.dat", [τ δζ δη W D])

	z, G_r, G´, G´´ = Mason(τ, W)
	save_data(folder*"/Viscoelasticity0F.dat", [z G_r G´ G´´ z])
	=#

	while true
		S = NESF(Li, Lf, G, HS, u; VerletWeis = VerletWeis)

		τ, Fs, F, Δζ, Δη, D, W = dynamics(Lf, G, S, k₀, 1e-10, 5, 60)
		#τ, Fs, F, Δζ, Δη, D, W = dynamics(Lf, G, S, k₀, 1e-10, 5, 80)

		tau_dump = 0.0
		for ii in 1:length(τ)-1
			if Fs[ii] > 1/exp(1) tau_dump = τ[ii] end
		end

		save_data(folder*"/input"*num2text(j)*".dat", [G.x Si])
		save_data(folder*"/fs"*num2text(j)*".dat", [τ Fs F Δζ D W Δη])

		z, G_r, G´, G´´ = Mason(τ, W)
		save_data(folder*"/Viscoelasticity"*num2text(j)*".dat", [z G_r G´ G´´ z*tau_dump z*tw])

		δζ, δη, W, D = MSD(τ, Δζ, Δη, phi_f)
		save_data(folder*"/prop_transp"*num2text(j)*".dat", [τ δζ δη W D])

		z, G_r, G´, G´´ = Mason(τ, W)
		save_data(folder*"/Viscoelasticity0"*num2text(j)*".dat", [z G_r G´ G´´ z*tau_dump z*tw])

		bI = mobility(τ, Δζ)

		u += Δt/bI

		append!(tw_save, tw)
		append!(b_save, 1/bI)
		append!(bi_save, bI)
		append!(tau_save, tau_dump)
		append!(u_save, u)
		append!(j_save, j)
		append!(Δeta_save, dη(τ, Δη))
		append!(eta_save, η∞(phi_f) + Δeta_save[end])
		save_data(folder*"/instantaneous_process.dat", [tw_save u_save bi_save b_save tau_save j_save eta_save Δeta_save])

		if tw > 1e8 break end
		if bI > 1e7 break end
		#if tw > 1e11 break end
		#if bI > 1e15 break end
		if Fs[end] > 1/exp(1) break end #Corta Δη para que no genere ruido
		tw += Δt
		j += 1
		Δt *= 1.5
	end

	# Arrest structure
	factor = 1.005
	Sf= NESF(Li, Lf, G, HS, factor*u_save[end]; VerletWeis = VerletWeis)

	τ, Fs, F, Δζ, Δη, D, W = dynamics(Lf, G, Sf, k₀, 1e-10, 5, 60)

	save_data(folder*"/inputA.dat", [G.x Si])
	save_data(folder*"/fsA.dat", [τ Fs F Δζ D W Δη])
	z, G_r, G´, G´´ = Mason(τ, W)
	save_data(folder*"/ViscoelasticityA.dat", [z G_r G´ G´´ z])

	δζ, δη, W, D = MSD(τ, Δζ, Δη, phi_f)
	save_data(folder*"/prop_transpA.dat", [τ δζ δη W D])

	z, G_r, G´, G´´ = Mason(τ, W)
	save_data(folder*"/Viscoelasticity0A.dat", [z G_r G´ G´´ z])

	#save_data(folder*"/instantaneous_process.dat", [tw_save u_save bi_save b_save tau_save j_save eta_save Δeta_save])
end
