using Texshade
# write your own tests here

x = randn(17,17)
XF = fftshift(fft(x))

mkf(n) = ((0 : (n - 1)) - (n >> 1)) / n

alpha = 0.5

y = real(ifft(ifftshift(map(t->(t[1]*t[1] + t[2]*t[2]) .^ (alpha / 2), Base.product(mkf(size(XF,1)), mkf(size(XF,2)))) .* XF)))

x2 = zeros(18,18)
x2[1:17,1:17] = x

y2 = real(ifft(ifftshift(map(t->(t[1]*t[1] + t[2]*t[2]) .^ (alpha / 2), Base.product(mkf(18), mkf(18))) .* fftshift(fft(x2)))))

import PyPlot
PyPlot.imshow(y2[1:17,1:17] - y);PyPlot.colorbar()

PyPlot.imshow(y2[1:17,1:17] - texshade(x));PyPlot.colorbar()


texshade(x)

dir = "/Users/ahmed.fasih/Downloads/srtm1arcsecond/aster-v2/030407688661871/unzipped/vrt/"
dir = "/Users/ahmed.fasih/Downloads/srtm1arcsecond/unzipped/vrt/"
aster = read(dir * "output.envi", Int16, (3601, 7201))'

const plt = PyPlot

plt.imshow(aster); plt.colorbar()

t = texshade(aster)

plt.imshow(t); plt.colorbar()
typeof(t)


write(dir * "tex.envi", t')