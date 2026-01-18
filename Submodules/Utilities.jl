module Utilities

#Define the Heaviside function.
function Heaviside(x)
  Dt = typeof(x)
  if x<0
    return zero(Dt)
  elseif x==0
    return convert(Dt, 1//2)
  else
    return one(Dt)
  end
end

end
