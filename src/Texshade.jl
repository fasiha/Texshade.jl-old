module Texshade

export texshade;

mkf(n::Integer) = ((0 : (n - 1)) - (n >> 1)) / n

function texshade(x; alpha::Real=0.5)
  n = map(x -> nextprod([2,3,5,7], x), size(x))
  if n == size(x)
    return core(x; alpha)
  end

  x2 = zeros(n)
  x2[1:size(x,1), 1:size(x,2)] = x
  y2 = core(x2; alpha=alpha)
  return y2[1:size(x,1), 1:size(x,2)]
end


function core(x; alpha::Real=0.5)
  mapped = map(t -> (t[1]*t[1] + t[2]*t[2]) .^ (alpha * 0.5), Base.product(mkf(size(x, 1)), mkf(size(x, 2))))
  return real(ifft(ifftshift(mapped .* fftshift(fft(x)))))
end

end # module
