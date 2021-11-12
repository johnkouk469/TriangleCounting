cd(@__DIR__)
using Pkg
Pkg.activate(".")

function triangle_counting(A)
	A .* ( A * A )
end

using MatrixMarket

A = MatrixMarket.mmread("s12.mtx")

# compute triangles
C = triangle_counting( A )
e = ones( size(A,1) )
c = C * e / 2

print(c)