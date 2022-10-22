# Función de distribución par analítica
#######################################################################################################
#		g(r)
#######################################################################################################

"""
`g(η, r)`

# Argumentos
`η::Real` Fracción de volumen.

`r::Real` distacia respecto de separación entre partículas.

Función de distribución par para una fracción de volumen dado η.
"""
function g_HS(x, r)

	power(a, b) = a^b

	function f(x)
		return 3.0 + 3.0*x - x*x
	end

	function yp(x)
		F = f(x)
		uno = (2.0*x*F)*(sqrt(2.0*(x^4)/(F*F) + 1.0) + 1.0)
		return power(uno, 1.0/3.0)
	end

	function ym(x)
		F = f(x)
		uno = (2.0*x*F)*(sqrt(2.0*x*x*x*x/(F*F) + 1.0) - 1.0)
		return power(uno, 1.0/3.0)
	end

	function zd(x)
		return yp(x) - ym(x)
	end

	function zs(x)
		return yp(x) + ym(x)
	end

	function A(x)
		return (-2.0*x + zd(x))/(1.0 - x)
	end

	function B(x)
		return (-2.0*x - 0.5*zd(x))/(1.0 - x)
	end

	function C(x)
		return 0.5*sqrt(3.0)*zs(x)/(1.0 - x)
	end

	#######################################################################################################
	#		First Shell
	#######################################################################################################

	function a1(x)
		z = zd(x)
		numerador = -2.0*x*(1.0 - x - 3.0*x*x) + (1.0 - 3.0*x - 4.0*x*x)*z + (1.0 + 0.5*x)*z*z
		denominador = 3.0*(2.0*x*x + z*z)*(1.0 - x)*(1.0 - x)
		return numerador/denominador
	end

	function a2(x)
		z = zd(x)
		numerador = x*(2.0 + 4.0*x - 3.0*x*x) - (1.0 - 3.0*x - 4.0*x*x)*z + (2.0 + x)*z*z
		denominador = 3.0*(2.0*x*x + z*z)*(1.0 - x)*(1.0 - x)
		return numerador/denominador
	end

	function a3(x)
		z = zd(x)
		numerador = (1.0 - 3.0*x - 4.0*x*x)*(4.0*x*x + z*z) + x*(2.0 - 5.0*x*x)*z
		denominador = sqrt(3.0)*zs(x)*(2.0*x*x + z*z)*(1.0 - x)*(1.0 - x)
		return numerador/denominador
	end

	function H1(x, r)
		return a1(x)*exp(A(x)*(r - 1.0)) + a2(x)*exp(B(x)*(r - 1.0))*cos(C(x)*(r - 1.0)) + a3(x)*exp(B(x)*(r - 1.0))*sin(C(x)*(r - 1.0))
	end

	function g1(x, r)
		g = 0.0
		if r < 1.0
			return 0.0
		else
			return H1(x, r)/r
		end
	end

	#######################################################################################################
	#		Second Shell
	#######################################################################################################

	function b1(x)
		z = zd(x)
		x2 = x*x
		x3 = x*x2
		x4 = x2*x2
		numerador = (2.0*x2 - z*z)*(1.0 - 6.0*x - 3.0*x2 + 20.0*x3 + 15.0*x4) + z*x*(16.0 + 24.0*x - 21.0*x2 - 13.0*x3 + 21.0*x4)
		denominador = (2.0*x2 + z*z)*(2.0*x2 + z*z)*(2.0*x2 + z*z)*(1.0 - x)*(1.0 - x)
		return -4.0*x*numerador/(3.0*denominador)
	end

	function b2(x)
		return -b1(x)
	end

	function b3(x)
		z = zd(x)
		w = zs(x)
		x2 = x*x
		x3 = x*x2
		x4 = x2*x2
		x5 = x2*x3
		x6 = x3*x3
		numerador = (4.0*x*x + z*z)*(2.0 - 10.0*x - 24.0*x2 + 30.0*x3 + 79.0*x4 + 21.0*x5 - 17.0*x6) + 2.0*z*x*(16.0 + 40.0*x - x2 - 50.0*x3 + 11.0*x4 + 52.0*x5 + 13.0*x6)
		denominador = w*w*w*(2.0*x*x + z*z)*(2.0*x*x + z*z)*(2.0*x*x + z*z)*(1.0 - x)*(1.0 - x)
		return 8.0*sqrt(3.0)*x2*numerador/denominador
	end

	function b4(x)
		z = zd(x)
		x2 = x*x
		x3 = x*x2
		x4 = x2*x2
		numerador = 2.0*x*(10.0 + 28.0*x + 21.0*x2 - 13.0*x3 - 19.0*x4) + z*(2.0 - 12.0*x - 18.0*x2 + 28.0*x3 + 27.0*x4) + z*z*(4.0 - 6.0*x - 18.0*x2 - 7.0*x3)
		denominador = (2.0*x2 + z*z)*(2.0*x2 + z*z)*(1.0 - x)*(1.0 - x)*(1.0 - x)
		return -2.0*x*numerador/(3.0*denominador)
	end

	function b5(x)
		z = zd(x)
		w = zs(x)
		x2 = x*x
		x3 = x*x2
		x4 = x2*x2
		x5 = x2*x3
		x6 = x3*x3
		numerador = 4.0*(6.0 - 30.0*x - 82.0*x2 + 58.0*x3 + 222.0*x4 + 94.0*x5 - 25.0*x6) - z*(24.0 - 10.0*x - 164.0*x2 - 156.0*x3 + 22.0*x4 +41.0*x5) - z*z*(10.0 + 32.0*x + 15.0*x2 - 31.0*x3 - 26.0*x4)
		denominador = w*w*(2.0*x2 + z*z)*(2.0*x2 + z*z)*(1.0 - x)*(1.0 - x)*(1.0 - x)
		return -4.0*x2*numerador/(3.0*denominador)	
	end

	function b6(x)
		z = zd(x)
		w = zs(x)
		x2 = x*x
		x3 = x*x2
		x4 = x2*x2
		x5 = x2*x3
		numerador = 24.0 - 10.0*x - 164.0*x2 - 156.0*x3 + 22.0*x4 + 41.0*x5 - z*(10.0 + 32.0*x + 15.0*x2 - 31.0*x3 - 26*x4)
		denominador = w*(2.0*x2 + z*z)*(2.0*x2 + z*z)*(1.0 - x)*(1.0 - x)*(1.0 - x)
		return 4.0*x2*numerador/(sqrt(3.0)*denominador)	
	end

	function H2(x, r)
		return b1(x)*exp(A(x)*(r - 2.0)) + b2(x)*exp(B(x)*(r - 2.0))*cos(C(x)*(r - 2.0)) + b3(x)*exp(B(x)*(r - 2.0))*sin(C(x)*(r - 2.0)) + b4(x)*(r - 2.0)*exp(A(x)*(r - 2.0)) + b5(x)*(r - 2.0)*exp(B(x)*(r - 2.0))*cos(C(x)*(r - 2.0)) + b6(x)*(r - 2.0)*exp(B(x)*(r - 2.0))*sin(C(x)*(r - 2.0))
	end

	function g2(x, r)
		g = 0.0
		if r < 2.0
			return 0.0
		else
			return H2(x, r)/r
		end
	end

	#######################################################################################################
	#		Third Shell
	#######################################################################################################

	function c1(x)
		z = zd(x)
		w = zs(x)
		x2 = x*x
		x3 = x*x2
		x4 = x2*x2
		x5 = x2*x3
		x6 = x3*x3
		x7 = x3*x4
		x8 = x4*x4
		numerador = (2.0*x2 - z*z)*(4.0 - 56.0*x - 136.0*x2 + 130.0*x3 + 343.0*x4 - 268.0*x5 - 665.0*x6 - 170.0*x7 + 89.0*x8) + z*x*(92.0 + 116.0*x - 744.0*x2 - 1468.0*x3 + 332.0*x4 + 1902.0*x5 + 251.0*x6 - 925.0*x7 - 258*x8)
		denominador = w*w*power(2.0*x2 + z*z, 5)*(1.0 - x)*(1.0 - x)
		return -32.0*x3*numerador/denominador
	end

	function c2(x)
		return -c1(x)
	end

	function c3(x)
		z = zd(x)
		w = zs(x)
		x2 = x*x
		x3 = x*x2
		x4 = x2*x2
		x5 = x2*x3
		x6 = x3*x3
		x7 = x3*x4
		x8 = x4*x4
		x9 = x4*x5
		x10 = x5*x5
		numerador = (4.0*x2 + z*z)*(24.0 - 312.0*x - 1252.0*x2 - 40.0*x3 + 3854.0*x4 + 1658.0*x5 - 6616.0*x6 - 6376.0*x7 + 593.0*x8 + 1799.0*x9 + 107.0*x10) + 2.0*z*x*(276.0 + 624.0*x - 1960.0*x2 - 6976.0*x3 - 3208.0*x4 + 8690.0*x5 + 7499.0*x6 - 4996.0*x7 - 6541.0*x8 - 610.0*x9 + 641.0*x10)
		denominador = w*w*w*w*w*power(2.0*x2 + z*z, 5)*(1.0 - x)*(1.0 - x)
		return 64.0*sqrt(3.0)*x4*numerador/denominador
	end

	function c4(x)
		z = zd(x)
		w = zs(x)
		x2 = x*x
		x3 = x*x2
		x4 = x2*x2
		x5 = x2*x3
		x6 = x3*x3
		x7 = x3*x4
		numerador = 2.0*x*(36.0 - 92.0*x - 360.0*x2 + 156.0*x3 + 982.0*x4 + 747.0*x5 + 195.0*x6 +37.0*x7) - 3.0*z*x*(4.0 + 8.0*x - 94.0*x2 - 60.0*x3 + 280.0*x4 + 256.0*x5 + 11.0*x6) + z*z*(4.0 - 60.0*x - 84.0*x2 + 190.0*x3 + 150.0*x4 - 258.0*x5 - 185.0*x6)
		denominador = power(2.0*x2 + z*z, 4)*power(1.0 - x, 3)
		return -8.0*x2*numerador/(3.0*denominador)
	end

	function c5(x)
		z = zd(x)
		w = zs(x)
		x2 = x*x
		x3 = x*x2
		x4 = x2*x2
		x5 = x2*x3
		x6 = x3*x3
		x7 = x3*x4
		x8 = x4*x4
		x9 = x4*x5
		x10 = x5*x5
		numerador = 96.0 - 1248.0*x - 4960.0*x2 - 256.0*x3 + 14976.0*x4 + 7232.0*x5 - 24376.0*x6 - 25632.0*x7 - 724.0*x8 + 5348.0*x9 + 384.0*x10 - z*(72.0 - 96.0*x - 1152.0*x2 - 932.0*x3 + 3012.0*x4 + 4920.0*x5 + 604.0*x6 - 2286.0*x7 - 708.0*x8 + 211.0*x9) + z*z*(12.0 - 24.0*x - 110.0*x2 + 150.0*x3 + 522.0*x4 - 32.0*x5 - 774.0*x6 - 462.0*x7 - 11.0*x8)
		denominador = w*w*w*w*power(2.0*x2 + z*z, 4)*power(1.0 - x, 3)
		return -32.0*x4*numerador/denominador
	end

	function c6(x)
		z = zd(x)
		w = zs(x)
		x2 = x*x
		x3 = x*x2
		x4 = x2*x2
		x5 = x2*x3
		x6 = x3*x3
		x7 = x3*x4
		x8 = x4*x4
		x9 = x4*x5
		numerador = 72.0 - 96.0*x - 1152.0*x2 - 932.0*x3 + 3012.0*x4 + 4920.0*x5 + 604.0*x6 - 2286.0*x7 - 708.0*x8 + 211.0*x9 + z*(12.0 - 24.0*x - 110.0*x2 + 150.0*x3 + 522.0*x4 - 32.0*x5 - 774.0*x6 - 462.0*x7 - 11.0*x8)
		denominador = w*w*w*power(2.0*x2 + z*z, 4)*power(1.0 - x, 3)
		return 32.0*sqrt(3)*x4*numerador/denominador
	end

	function c7(x)
		z = zd(x)
		w = zs(x)
		x2 = x*x
		x3 = x*x2
		x4 = x2*x2
		x5 = x2*x3
		x6 = x3*x3
		numerador = 4.0*x*(34.0 + 18.0*x - 195.0*x2 - 332.0*x3 - 144.0*x4 + 69.0*x5 + 64.0*x6) + 2.0*z*(2.0 + 6.0*x + 30.0*x2 + 176.0*x3 + 174.0*x4 - 66.0*x5 - 79.0*x6) + 3.0*z*z*(4.0 - 14.0*x - 32.0*x2 + 32.0*x3 + 68.0*x4 + 23.0*x5)
		denominador =power(2.0*x2 + z*z, 3)*power(1.0 - x, 4)
		return 2.0*x2*numerador/(3.0*denominador)
	end

	function c8(x)
		z = zd(x)
		w = zs(x)
		x2 = x*x
		x3 = x*x2
		x4 = x2*x2
		x5 = x2*x3
		x6 = x3*x3
		x7 = x3*x4
		x8 = x4*x4
		numerador = 12.0 + 48.0*x + 280.0*x2 + 1224.0*x3 + 1776.0*x4 - 80.0*x5 - 1794.0*x6 - 816.0*x7 + 79.0*x8 + z*(36.0 - 88.0*x - 420.0*x2 + 72.0*x3 + 1172.0*x4 + 897.0*x5 - 63.0*x6 - 148.0*x7) - z*z*(68.0 + 48.0*x - 432.0*x2 - 760.0*x3 - 192.0*x4 + 342.0*x5 + 197.0*x6)
		denominador = w*w*power(2.0*x2 + z*z, 3)*power(1.0 - x, 4)
		return -8.0*x3*numerador/(3.0*denominador)
	end

	function c9(x)
		z = zd(x)
		w = zs(x)
		x2 = x*x
		x3 = x*x2
		x4 = x2*x2
		x5 = x2*x3
		x6 = x3*x3
		x7 = x3*x4
		x8 = x4*x4
		numerador = (4.0*x2 + z*z)*(36.0 - 88.0*x - 420.0*x2 + 72.0*x3 + 1172.0*x4 + 897.0*x5 - 63.0*x6 - 148.0*x7) - z*(12.0 + 48.0*x + 416.0*x2 + 1320.0*x3 + 912.0*x4 - 1600.0*x5 - 2178.0*x6 - 132.0*x7 + 473.0*x8)
		denominador = w*w*w*power(2.0*x2 + z*z, 3)*power(1.0 - x, 4)
		return 8.0*x3*numerador/(sqrt(3.0)*denominador)
	end

	function H3(x, r)
		return c1(x)*exp(A(x)*(r - 3.0)) + c2(x)*exp(B(x)*(r - 3.0))*cos(C(x)*(r - 3.0)) + c3(x)*exp(B(x)*(r - 3.0))*sin(C(x)*(r - 3.0)) + c4(x)*(r - 3.0)*exp(A(x)*(r - 3.0)) + c5(x)*(r - 3.0)*exp(B(x)*(r - 3.0))*cos(C(x)*(r - 3.0)) + c6(x)*(r - 3.0)*exp(B(x)*(r - 3.0))*sin(C(x)*(r - 3.0)) + c7(x)*(r - 3.0)*(r - 3.0)*exp(A(x)*(r - 3.0)) + c8(x)*(r - 3.0)*(r - 3.0)*exp(B(x)*(r - 3.0))*cos(C(x)*(r - 3.0)) + c9(x)*(r - 3.0)*(r - 3.0)*exp(B(x)*(r - 3.0))*sin(C(x)*(r - 3.0))
	end

	function g3(x, r)
		g = 0.0
		if r < 3.0
			return 0.0
		else
			return H3(x, r)/r
		end
	end

	return g1(x, r) + g2(x, r) + g3(x, r)
end
