"""
    extract_data(filename,cols = [col1 col2 col3 ... coln],comment = "#")

Extracts the vectors column 1,column 2, .. and column n  from the file given by
filename. Ignores the lines prefixed with the comment.
cols and comment are keyword argumens with default values.
Returns Nx2 where N is the length of each column vector in the file
"""
function extract_data(filename::String;cols::Array{Int64,2} = Array{Int64,2}(undef,0,0),comment::String = "")
    # Open the file and parse it into lines
    file = open(filename)
    lines = readlines(file)

    # Create empty n length vector to fill with data horizontally
    # n is the numbers of columns requested
    returnArray = nothing

    # Iterate through all the lines
    for line in lines
        # Split on comment characters to only grab the part we want
        # The first chunk has our data
        # Only do this if we specified a comment character
        if comment != ""
            lineWithoutComments = split(line,comment)[1]
            # If this chunk is empty, the entire line was a comment
            if lineWithoutComments == ""
                continue
            end
            # Parse to Float64s from number strings seperated by spaces,
            # stripping trailing and leading whitespace
            # If this is the first one
            if returnArray == nothing
                returnArray = [parse(Float64,n) for n in split(strip(lineWithoutComments)," ")]'
            else
                returnArray = [returnArray;[parse(Float64,n) for n in split(strip(lineWithoutComments)," ")]']
            end
        else # For non-commented data input
            # If this is the first one
            if returnArray == nothing
                returnArray = [parse(Float64,n) for n in split(strip(line)," ")]'
            else
                returnArray = [returnArray;[parse(Float64,n) for n in split(strip(line)," ")]']
            end
        end
    end

    # If we are extracting certain columns, not the whole thing
    if length(cols) != 0
        outArray = nothing
        for col in cols
            if outArray == nothing
                outArray = returnArray[:,col] # Set the first column to the first requested column
            else
                outArray = hcat(outArray,returnArray[:,col]) # Append additional columns
            end
        end
        return outArray
    else # Else return the whole thing
        return returnArray
    end
end

"""
    bin_data(data,nbins)

Breaks up data into nbins equal bins and counts how many data points fall into
each bin for use in a histogram. bin_data returns bin midpoints and the counts
per bin.

nbins is optional argument, defaults to 10 bins
"""
function bin_data(data::Vector,nbins::Int64=10)
    # Determine range of data
    data_min = minimum(data)
    data_max = maximum(data)

    # These bins will be lower-inclusive, as in the data will be binned to
    # lower bin first
    range = data_max-data_min
    bin_size = range/nbins

    # Initialize bin counts to zero
    counts = zeros(Int64,nbins)

    for datum in data
        # Zero the data
        thisDatum = datum-data_min
        for i = 1:nbins
            if bin_size*i >= thisDatum
                counts[i] += 1
                break
            end
        end
    end

    # Calculate midpoints
    midpoints = []
    for i = 1:nbins
        append!(midpoints,data_min + (bin_size*i)/2)
    end

    return midpoints,counts
end

using Printf
function print_hist(bins,counts,ch::String="#",ch_N=-1)
    # If no ch_N is given, calculate optimum
    if ch_N == -1
        # I will use the smallest non-zero count to show maximum precision
        for count in sort(counts)
            if count != 0
                ch_N = count
                break
            end
        end
    end
    # ASCII Print
    println("Number of each $ch represents a count of $ch_N")
    println("\tCount\t|\tMidpoint\t|")
    println("-----------------------------")
    for (i,bin) in enumerate(bins)
        print("\t$(counts[i])\t\t|\t")
        @printf "%8.3f\t|" bin

        # Print number of characters/ch_N
        for i = 1:floor(Int64,counts[i]/ch_N)
            print(ch)
        end
        println("")
    end
    println("")
end

#### Script to run everything

# cd("/Users/kiranshila/Desktop/USF/USF Spring 2019/Numerical Methods/")
# This script should be run in the folder with the data, otherwise cd to the
# directory with the data in the REPL.
data = extract_data("data_w_comments.txt",cols = [1 3 5],comment = "#")
(midpoints,counts) = bin_data(data[:,1],25)
print_hist(midpoints,counts,"*",4)
