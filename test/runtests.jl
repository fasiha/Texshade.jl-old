using Texshade
reload("Texshade")
using Base.Test

mkf(n) = ((0 : (n - 1)) - (n >> 1)) / n
alpha = 0.5

x = randn(17,17)
XF = fftshift(fft(x))

y = real(ifft(ifftshift(map(t->(t[1]*t[1] + t[2]*t[2]) .^ (alpha / 2), Base.product(mkf(size(XF,1)), mkf(size(XF,2)))) .* XF)))

@test all(x->isapprox(x,0), Texshade.core2(x, alpha) - y)
@test all(x->isapprox(x,0), Texshade.core(x, alpha) - y)

x2 = zeros(18,18)
x2[1:17,1:17] = x
y2 = real(ifft(ifftshift(map(t->(t[1]*t[1] + t[2]*t[2]) .^ (alpha / 2), Base.product(mkf(18), mkf(18))) .* fftshift(fft(x2)))))

@test all(x->isapprox(x,0), Texshade.texshade(x; alpha=alpha) - y2[1:17, 1:17])

@test all(x->isapprox(x,0), Texshade.core2(x2, alpha) - y2)

@test all(x->isapprox(x,0; atol=2 * eps()), Texshade.core3(x2, alpha) - y2)

x = randn(3000,3000)

using BenchmarkTools
@benchmark Texshade.core(x, alpha)
@benchmark Texshade.core2(x, alpha)
@benchmark Texshade.core3(x, alpha)

@time Texshade.core(x, alpha);
@time Texshade.core2(x, alpha);
@time Texshade.core3(x, alpha);
@time Texshade.innercore3!(x, alpha);
