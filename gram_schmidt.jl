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
