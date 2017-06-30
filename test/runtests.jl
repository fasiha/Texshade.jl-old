using Texshade
using Base.Test

x = randn(17,17)
XF = fftshift(fft(x))

mkf(n) = ((0 : (n - 1)) - (n >> 1)) / n

alpha = 0.5

y = real(ifft(ifftshift(map(t->(t[1]*t[1] + t[2]*t[2]) .^ (alpha / 2), Base.product(mkf(size(XF,1)), mkf(size(XF,2)))) .* XF)))

x2 = zeros(18,18)
x2[1:17,1:17] = x

y2 = real(ifft(ifftshift(map(t->(t[1]*t[1] + t[2]*t[2]) .^ (alpha / 2), Base.product(mkf(32), mkf(32))) .* fftshift(fft(x2)))))

@test 1 == 2


texshade(x)-y2[1:17,1:17]