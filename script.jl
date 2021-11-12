function triangle_counting(A)
	A .* ( A * A )
end

 
A = adjacency_matrix( smallgraph( :karate ) )  # load the karate club
	# compute triangles
C = triangle_counting( A )
e = ones( size(A,1) )
c = C * e / 2

# plot result
plot_triangles( A, C, c; 
title = "#triangles incident with each edge/node on Karate club",
edge_color = colorant"#0008", node_color = colorant"#0008" )

# example of a graph from the web
A = tamudata("Arenas", "email")

# compute triangles
C = triangle_counting( A )
e = ones( size(A,1) )
c = C * e / 2

plot_triangles( A, C, c; 
    title = "#triangles incident with each edge/node on a real network from the web" )
	
