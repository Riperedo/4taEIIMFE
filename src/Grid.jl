#####################
#	Grid Library 	#
#####################

mutable struct Grid <: object
	# atributos
	n :: Integer
	x₀ :: Real
	xₙ :: Real
	grid :: Array{Tuple{Float64, Float64}}
	uniform :: Bool
	x :: Array{Float64}
	w :: Array{Float64}
	# métodos
	setn :: Function
	setx₀ :: Function
	setxₙ :: Function
	setGrid :: Function
	addGrid :: Function
	find :: Function
	function Grid(x0 :: Real, xn :: Real, N :: Integer)
		@assert x0 < xn "El primer arguento debe ser menor que el segundo."
		this = new()
		this.n = N
		this.x₀ = x0
		this.xₙ = xn
		dx = (xn - x0)/(N-1.0)
		index = collect(1:N)
		grid = collect(x0:dx:xn)
		if length(grid) < N
		    append!(grid, xn)
		end
		weight = dx*ones(N)
		this.x = grid
		this.w = weight
		this.grid = map(i -> (grid[i], weight[i]), index)
		this.uniform = true
		# métodos públicos
		this.setx₀ = function (x :: Real)
			this.x₀ = x
			reload()
		end
		this.setxₙ = function (x :: Real)
			this.xₙ = x
			reload()
		end
		this.setn = function (x :: Integer)
			this.n = x
			reload()
		end
		this.setGrid = function (grid :: Array{Float64}, weight :: Array{Float64})
			@assert length(grid) == length(weight) "Arrays de diferente tamano"
			this.x₀ = grid[1]
			this.xₙ = grid[end]
			this.n = length(grid)
			index = collect(1:this.n)
			dx = (this.xₙ - this.x₀)/(this.n-1.0)
			this.x = grid
			this.w = weight
			this.grid = map(i -> (grid[i], weight[i]), index)
			checkUniform()
		end
		this.setGrid = function (grid :: Array{Float64})
			weight = zeros(length(grid))
			for i in 2:length(grid)
				weight[i] = grid[i+1]-grid[i]
			end
			weight[1] = weight[2]
			this.x₀ = grid[1]
			this.xₙ = grid[end]
			this.n = length(grid)
			index = collect(1:this.n)
			dx = (this.xₙ - this.x₀)/(this.n-1.0)
			this.x = grid
			this.w = weight
			this.grid = map(i -> (grid[i], weight[i]), index)
			checkUniform()
		end
		this.addGrid = function (x :: Grid)
			@assert x.x₀ < x.xₙ "Tu grid debe estar ordenado."
			x0 = x.x₀
			xn = x.xₙ
			tx0 = this.x₀
			txn = this.xₙ
			# Caso 1
			# x0   xn  tx0            txn
			# /----/   |--------------|
			if xn < tx0
				x.w[end] = this.x[1] - x.x[end]
				this.grid = concatenate(x.grid, this.grid)
			end
			# Caso 2
			# tx0           txn  x0   xn
			# |--------------|   /----/
			if x0 >= txn
				this.w[end] = x.x[1] - this.x[end]
				this.grid = concatenate(this.grid, x.grid)
			end
			# Caso 3
			# tx0   x0   xn  txn
			# |-----/----/---|
			if tx0 < x0 < txn && tx0 < xn < txn 
				grid = pop(x0, xn, this.grid)
				this.grid = concatenate(grid, x.grid)
				sort!(this.grid, by = x -> x[1])
			end
			# Caso 4
			# x0 tx0  xn  txn
			# /---|---/---------|
			if x0 < tx0 < xn && tx0 < xn < txn 
				grid = pop(x0, xn, this.grid)
				this.grid = concatenate(x.grid, grid)
			end
			# Caso 5
			# tx0       x0 txn  xn
			# |---------/---|---/ 
			if tx0 < x0 < txn && x0 < txn < xn 
				grid = pop(x0, xn, this.grid)
				this.grid = concatenate(grid, x.grid)
			end
			# Caso 6
			# x0   tx0  txn  xn
			# /-----|----|---/
			if x0 < tx0 < xn && x0 < txn < xn 
				this.setGrid(x.x, x.w)
			end
			this.x₀ = this.grid[1][1]
			this.xₙ = this.grid[end][1]
			this.n = length(this.grid)
			getArrays()
			checkUniform()
		end
		this.find = function (xᵢ)
			Xi = (this.x - xᵢ*ones(this.n)).^2
			x = minimum(Xi)
			index = 1
			for i in 1:this.n
				if Xi[i] == x
					index = i
				end
			end
			return index
		end
		# metodos privados
		function reload()
			dx = (this.xₙ - this.x₀)/(this.n-1.0)
			grid = collect(this.x₀:dx:this.xₙ)
			weight = dx*ones(this.n)
			this.setGrid(grid, weight)
		end
		function concatenate(A :: Array{Tuple{Float64, Float64}}, B :: Array{Tuple{Float64, Float64}})
			if size(A)[1] == 1 && size(B)[1] == 1 # arreglos horizontales
				return hcat(A, B)
			else # arreglo vertical
				return vcat(A, B)
			end
		end
		function pop(x0 :: Float64, xn :: Float64, A :: Array{Tuple{Float64, Float64}})
			index = collect(1:length(A))
			grid = map(j -> A[j], filter(i -> !(x0 <= A[i][1] <= xn), index))
			return grid
		end
		function checkUniform()
			dx = (this.xₙ - this.x₀)/(this.n - 1.0)
			uniform = collect(this.x₀:dx:this.xₙ)
			weight = dx*ones(this.n)
		    if this.x == uniform && this.w == weight 
				this.uniform = true
			else
				this.uniform = false
			end
		end
		function getArrays()
			index = collect(1:length(this.grid))
		    this.x = map(i -> this.grid[i][1], index)
		    this.w = map(i -> this.grid[i][2], index)
		end
		return this
	end
end

function T(n, x)
	k = collect(1:n)
	polinomio = 2.0^(n-1)
	for i in 1:n
		polinomio *= x - cospi((k[i] + 0.5)/n)
	end
	return polinomio
end

function ChebyshevGrids(a :: Real, b :: Real, N :: Integer)
	@assert a < b "el primer argumento debe ser menor que el primero"
	i = collect(1:N)
	x = map(i -> 0.5*(a+b) + 0.5*(b-a)*cospi((2*i-1)/(2*N)), i)
	y = [-1.0]
	append!(y, map(i -> cospi((2*i-1)/(2*N)), i))
	append!(y, 1.0)
	w = map(y -> π*0.5*(b-a)*sqrt(1.0 - y*y)/N, y)
	append!(x, a)
	append!(x, b)
	sort!(x)
	return x, w
end

function integral(grid :: Grid, F :: Array; rule = "rectangle" :: String)
	lista_de_metodos = ["rectangle", "trapezoidal"]#, "Simpson 1/3"]
	@assert rule in lista_de_metodos "Los métodos programados se encuentran en la siguiente lista: "*string(metodos)
	@assert length(F) == grid.n "La función a integrar debe tener la misma longitud que el grid."
	if rule == "rectangle"
		I = 0.0
		for i in 1:grid.n-1
		    I += F[i]*grid.w[i]
		end
		return I
	end
	if rule == "trapezoidal"
		I = 0.0
		for i in 1:grid.n-1
		    I += 0.5*(F[i] + F[i+1])*grid.w[i]
		end
		return I
	end
	if rule == "Simpson 1/3"
		# no est'a bien programado
		I = 0.0
		N = grid.n - 1
		if grid.uniform
			for i in 1:2:round(Int, N/2)
				a = grid.x[i]
				b = grid.x[i+2]
				I += (b-a)*(F[i] + 4*F[i+1] + F[i+2])/6.0
			end
			return I
		else
			for i in 1:2:round(Int, N/2)
				a = grid.x[i]
				m = grid.x[i+1]
				b = grid.x[i+2]
				Ia = (b^2 + b*a + a^2 - 1.5*(m+b)*(b+a) + 3.0*m*b)/(3.0*(m-a))
				Im = (b-a)*(b^2 + a*b + a^2 - 1.5*(b^2+2*a*b+a^2) + 3.0*a*b)/(3.0*(m-b)*(m-a))
				Ib = (b^2 + b*a + a^2 - 1.5*(a+m)*(b+a) + 3.0*a*m)/(3.0*(b-m))
				I += F[i]*Ia + F[i+1]*Im + F[i+2]*Ib
			end
			return I
		end
	end
end
