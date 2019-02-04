using LinearAlgebra
using Random
using BenchmarkTools
using Profile
using ASTInterpreter2

function gauss_elim(A::Matrix,b::Vector)
    # Create augmented matrix
    # make sure its a float otherwise we get float errors
    n = size(b)[1]
    A_Aug = float([A b])
    # For every column except the last one because pivots
    for k = 1:n-1
        # Find the largest pivot in this column
        location = findmax(broadcast(abs,A_Aug[k:end,k]))[2] + k - 1
        maxVal = A_Aug[location,k]
        println("Row $k maxval is $maxVal at $location")
        # Check to see if current pivot is the max,
        # if it isn't - swap
        if A_Aug[k,k] != maxVal
            for j in k:n+1
                temp = A_Aug[location,j]
                A_Aug[location,j] = A_Aug[k,j]
                A_Aug[k,j] = temp
            end
        end
        # Now perform gaussian elimination
        # For every row under the pivot
        for i = k+1:n
            # Normalize to the pivot and subtract from pivot row
            if A_Aug[i,k] != 0 # Ensures that we need to perform elimination
                scale = A_Aug[k,k] / A_Aug[i,k]
                for j in k:n+1 # Start at this column to the end including augmentation
                    println("Modifying [$i,$j]")
                    A_Aug[i,j] = A_Aug[i,j] * scale
                    A_Aug[i,j] = A_Aug[i,j] - A_Aug[k,j]
                    display(A_Aug)
                    if abs(A_Aug[i,j]) < eps()
                        A_Aug[i,j] = 0 # For stability against epsilon
                    end
                end
            end
        end
    end
    # Now back substitute to get to x
    x = zeros(Float64,n) # Create zeros for solution vector
    x[n] = A_Aug[end,end] / A_Aug[n,n] # Set the starting x
    # For every row, working backwards starting with the second from the bottom
    for i in n-1:-1:1
        x[i] = (A_Aug[i,end] - sum( [x[j] * A_Aug[i,j] for j in i+1:n] )) / A_Aug[i,i]
    end
    # Return solution vector x and U
    return x, A_Aug[:,1:end-1]
end

A = [2 6 10;1 3 3;3 14 28]
b = [0;2;-8]
x,U = gauss_elim(A,b)

function my_lu(A::Matrix)
    # Create copy
    # make sure its a float otherwise we get float errors
    n = size(A)[1]
    A_Copy = float(A)
    L = float(one(A))
    P = one(A)
    # For every column except the last one because pivots
    for k = 1:n-1
        # Find the largest pivot in this column
        location = findmax(broadcast(abs,A_Copy[k:end,k]))[2] + k - 1
        maxVal = A_Copy[location,k]
        # Check to see if current pivot is the max,
        # if it isn't - swap
        if A_Copy[k,k] != maxVal
            for j in k:size(A_Copy)[2]
                temp = A_Copy[location,j]
                A_Copy[location,j] = A_Copy[k,j]
                A_Copy[k,j] = temp
            end
            # And swap in P
            P[location,:] , P[k,:] = P[k,:] , P[location,:]
            # And swap columns in L
            L[:,location] , L[:,k] = L[:,k] , L[:,location]
        end

        # Now perform gaussian elimination
        # For every row under the pivot
        for i = k+1:n
            # Normalize to the pivot and subtract from pivot row
            if A_Copy[i,k] != 0 # Ensures that we need to perform elimination
                scale =  A_Copy[i,k] / A_Copy[k,k]
                thisL = float(one(A))
                thisL[i,k] = scale
                L = L * thisL
                for j in k:size(A_Copy)[2]
                    A_Copy[i,j] = A_Copy[i,j] - scale * A_Copy[k,j]
                    if abs(A_Copy[i,j]) < eps()
                        A_Copy[i,j] = 0 # For stability against epsilon
                    end
                end
            end
        end
    end
    # Return solution vector L and U
    return P*L,A_Copy,P
end

L,U,P = my_lu(A)

A = [1 2 1;3 8 1;0 4 1]

L,U,P = my_lu(A)

L,U = lu(A)
