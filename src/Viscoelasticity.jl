η_∞(ϕ) = (1+1.5*ϕ*(1+ϕ-0.189*ϕ*ϕ))/(1-ϕ*(1+ϕ-0.189*ϕ*ϕ))

function MSD(τ, Δζ, Δη, ϕ)
	#= calculando MSD de la forma
	        <Δr²(τ)>    τ             τ         D₀
	W(τ) = ---------- = ∫ dt' D(t') = ∫ dt' ----------
	          2n        0             0     1 + δζ(t')
	donde   τ
	δη(τ) = ∫dt'Δη(t')/η₀
			0
	----------------------------
		      k_BT	         D₀3πη₀σ		   D₀
	D(τ) = --------- = ----------------- = ----------------
		    3πη(τ)σ 	3π[η∞ + δη(t)]σ     (η∞/η₀ + δη(t)/η₀)
	donde   t
	δη(t) = ∫dt'Δη(t')
			0
	=#
	W = zeros(length(τ))
	δζ = zeros(length(τ))
	D = zeros(length(τ))
	δη = zeros(length(τ))
	for ii in 1:length(τ)-1
		dt = τ[ii+1] - τ[ii]
		δζ[ii+1] = δζ[ii] + dt*0.5*(Δζ[ii] + Δζ[ii+1])
		δη[ii+1] = δη[ii] + dt*0.5*(Δη[ii] + Δη[ii+1])
		W[ii+1] = W[ii] + dt*0.5*(1/(1+δζ[ii+1])+1/(1+δζ[ii]))
		D[ii+1] = 1/(1+δζ[ii+1])
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
		#if τ[i] == τ[i+1] continue end
		#println(i, " ", log(Δr²[i-1])-log(Δr²[i]), " ", log(τ[i-1])-log(τ[i]), " ", α)
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

#=
δζ, δη, W, D = MSD(τ, Δζ, Δη, ϕ)

save_data("DAT\\prop_transp"*num2text(ϕ)*".dat", [τ δζ δη W D])

z, G_r, G´, G´´ = Mason(τ, W)

save_data("DAT\\Viscoelasticity"*num2text(ϕ)*".dat", [z G_r G´ G´´])
=#