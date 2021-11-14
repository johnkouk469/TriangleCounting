function sequential_masked_triangle_counting(A)
    rows = rowvals(A)
    vals = nonzeros(A)
    m, n = size(A)
    C = zeros(n)
    for j = 1:n
        for k in nzrange(A, j)
            i = rows[k]
            v = vals[k]
            # perform operations with A[i,j] = v
            r = rows[nzrange(A, i)]
            c = rows[nzrange(A, j)]
            for l = 1:size(r,1)
                for w = 1:size(c,1)                    
                    if r[l] == c[w]
                        C[j] = C[j]+1
                        break
                    elseif r[l] < c[w]
                        break
                    end
                end
            end
        end
    end
    return C/2
end

using MatrixMarket

A = MatrixMarket.mmread("com-Youtube.mtx")
# print(sequential_masked_triangle_counting(A))
# rows = rowvals(A)
# vals = nonzeros(A)
# m, n = size(A)
# print(rows)
# for i = 1:n
#     println(nzrange(A, i))
# end
# for i =1:n
#     println((rows[nzrange(A, i)]))
# end

using BenchmarkTools
@benchmark sequential_masked_triangle_counting(A)