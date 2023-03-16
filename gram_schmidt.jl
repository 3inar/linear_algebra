# Find an orthonormal basis for the polynomials with degree < 3

using SymPy

@vars x

# defines an inner product
function inner(p, q)
  integrate(p*q, (x, 0, 1))
end

function norm(p)
  sqrt(inner(p, p))
end

basis = [1, x, x^2]

# Performs the Gram-Schmidt procedure on this basis
e1 = basis[1]

e2 = basis[2] - inner(basis[2], e1)*e1
e2 = e2/norm(e2)

e3 = basis[3] - inner(basis[3], e2)*e2 - inner(basis[3], e1)*e1
e3 = e3/norm(e3)

orthonormal_basis = [e1, e2, e3]

println("Orthonormal basis for P_2(R):")
println(orthonormal_basis[1])
println(orthonormal_basis[2])
println(orthonormal_basis[3])


# Riesz representation Theorem: given a linear functional phi there is a unique
# vector q such that phi(p) = <p, q> for all p, <> being the inner product.
# We can find q easily using an orthonormal basis

function phi(p)
  integrate(p*cos(pi*x), (x, 0, 1))
end

q = phi(e1)*e1 + phi(e2)*e2 + phi(e3)*e3

