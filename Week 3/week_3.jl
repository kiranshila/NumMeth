function gauss_elim(A::Matrix,b::Vector)
    # Create augmented matrix
    A_Aug = float([A b])
    # For every column except the last one because we augmented
    for k = 1:size(A_Aug)[2]-1
        # Find the largest pivot in this column
        maxVal, location = findmax(A_Aug[k:end,k])
        # Check to see if current pivot is the max,
        # if it isn't - swap
        if A_Aug[k,k] != maxVal
            tempRow = A_Aug[location,:] # Grab the row we want to swap
            A_Aug[location,:] = A_Aug[k,:] # Swap with the pivot row
            A_Aug[k,:] = tempRow # Place the swap row into the pivot row
        end
        # Now perform gaussian elimination
        # For every row under the pivot
        for i = k+1:size(A_Aug)[1]
            # Normalize to the pivot and subtract from pivot row
            A_Aug[i,:] = (A_Aug[k,k] / A_Aug[i,k]) .* A_Aug[i,:]
            A_Aug[i,:] -= A_Aug[k,:]
        end
    end
    # Return solution vector x and U
    return A_Aug[:,end], A_Aug[:,1:end-1]
end
