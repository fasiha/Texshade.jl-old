module Texshade

export texshade;

mkf(n::Integer) = ((0 : (n - 1)) - (n >> 1)) / n
mkf2(n::Integer) = ifftshift(((0 : (n - 1)) - (n >> 1)) / n)

function texshade(x; alpha::Real=0.5)
  n = map(x -> nextprod([2,3,5,7], x), size(x))
  if n == size(x)
    return core(x, alpha)
  end

  x2 = zeros(n)
  x2[1:size(x,1), 1:size(x,2)] = x
  y2 = core(x2, alpha)
  y2[1:size(x,1), 1:size(x,2)]
end


function core(x, alpha)
  mapped = map(t -> (t[1]*t[1] + t[2]*t[2]) .^ (alpha * 0.5), Base.product(mkf(size(x, 1)), mkf(size(x, 2))))
  real(ifft(ifftshift(mapped .* fftshift(fft(x)))))
end

function core2(x, alpha)
  xf = fft(x)
  for (i, (f1, f2)) in Base.enumerate(Base.product(mkf2(size(x, 1)), mkf2(size(x, 2))))
    xf[i] *= (f1 * f1 + f2 * f2) .^ (alpha * 0.5)
  end
  real(ifft(xf))
end

function innercore3!(x, alpha)
  assert(all(x -> x % 2 == 0, size(x)))
  (n1, n2) = size(x)
  d1 = 1.0 / n1
  alpha = 0.5 * alpha

  for i2 = 1 : n2
    local f2sq = if (i2 <= (n2 ÷ 2))
      (i2 - 1) / n2 else
      (n2 - i2 + 1) / -n2 end # superfluous minus sign since it's going to be ^2
    f2sq *= f2sq

    local f1 = 0.0
    for i1 in 1 : (n1 ÷ 2)
      x[i1, i2] *= (f1 * f1 + f2sq) .^ alpha
      f1 += d1
    end

    f1 = -0.5
    for i1 in (n1 ÷ 2 + 1) : n1
      x[i1, i2] *= (f1 * f1 + f2sq) .^ alpha
      f1 += d1
    end
  end
end

function core3(x, alpha)
  xf = fft(x)
  innercore3!(xf, alpha)
  ifft!(xf)
  real(xf)
end

end # module
