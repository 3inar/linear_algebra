# find the best polynomial in P_5(R) to approximate sin(x) in C_R[-\pi, \pi]
# with the inner product <f, g> = \int_{-\pi}^\pi f(x)g(x) dx
#
# We know that this will be the orthogonal projection on P_5(R), and we know
# that if we have an orthonormal basis e_0, ..., e_5 this projection will be
# simply P_U v = <v, e_0>e_0 + ... + <v, e_5>e_5

using SymPy
@vars x


# defines an inner product
function inner(p, q)
  integrate(p*q, (x, -PI, PI))
end

function norm(p)
  sqrt(inner(p, p))
end

basis = [1, x, x^2, x^3, x^4, x^5]
orthonormal = deepcopy(basis)

# does Gram-Schmidt
for i = 1:6
  orthonormal[i] = basis[i]
  if i > 0
    orth_p = 0;
    for j = 1:(i-1)
      orth_p += inner(basis[i], orthonormal[j])*orthonormal[j]
    end
    orthonormal[i] -= orth_p
  end
  orthonormal[i] = orthonormal[i]/norm(orthonormal[i])
end

v = sin(x)


# project v on P_5(R0
u = Sym(0)
for i=1:6
  u += inner(v, orthonormal[i])*orthonormal[i]
end

# get a numerical value instead of a tricy expression invarious powers of pi
u.evalf(5)
 
# Returns:
#            5            3            
# 0.0056441⋅x  - 0.15528⋅x  + 0.98789⋅x
