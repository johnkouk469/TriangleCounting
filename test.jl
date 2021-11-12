cd(@__DIR__)
using Pkg
Pkg.activate(".")

using HypertextLiteral
using GraphMakie
using CairoMakie
using LinearAlgebra
using SparseArrays
using Graphs
using SGtSNEpi
using MAT
using Colors
using BenchmarkTools
using FLANN
export FLANN

function plot_triangles( A, C, c;
	title = "#triangles",
	node_color = colorant"#0002",
	edge_color = colorant"#0002" )
	
	g = Graph( A )
	node_size = c
	node_size = node_size .- minimum(node_size)
	node_size = node_size ./ maximum(node_size)
	node_size = node_size .* 20 .+ 8 
	edge_width = [C[e.src, e.dst] + 1 for e in edges(g)]
	edge_width = log2.( edge_width .+ 1 )
	edge_width = edge_width .- minimum(edge_width)
	edge_width = edge_width ./ maximum(edge_width)
	edge_width = edge_width .* 4 .+ 0.1 
		
	
	f = Figure(); 
	ax = Axis( f[1,1]; 
		aspect = DataAspect(),
		title ); 
	hidedecorations!( ax )
	hidespines!( ax )
	
	Y = sgtsnepi( adjacency_matrix(g) )
	layout = _ -> Point.( zip( Y[:,1], Y[:,2] ) )
	
	graphplot!( ax, g; layout, node_size, edge_width, 
		node_color, edge_color )

end

function triangle_counting(A)
	A .* ( A * A )
end

using MatrixMarket

A = MatrixMarket.mmread("s12.mtx")

# compute triangles
C = triangle_counting( A )
e = ones( size(A,1) )
c = C * e / 2

# plot result
plot_triangles( A, C, c; 
    title = "#triangles incident with each edge/node on Karate club",
    edge_color = colorant"#0008", node_color = colorant"#0008" )