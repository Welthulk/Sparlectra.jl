# Author: Udo Schmitz (https://github.com/Welthulk)
# Date: 08.05.2023
# include-file nbi.jl

function getNBI(nodeNumberVec, branchTupleVec)
  num_nodes = length(nodeNumberVec)
  num_branches = length(branchTupleVec)

  #Initialize the NBI matrix with zeros
  NBI_matrix = zeros(Int, num_nodes, num_branches)

  #Fill the NBI matrix based on the branchTupleVec values
  for (j, branch) in enumerate(branchTupleVec)
    from_bus, to_bus = branch

    #Find the index positions of the buses involved
    from_bus_index = findfirst(x -> x == from_bus, nodeNumberVec)
    to_bus_index = findfirst(x -> x == to_bus, nodeNumberVec)

    #Set the appropriate entries in the NBI matrix
    NBI_matrix[from_bus_index, j] = -1
    NBI_matrix[to_bus_index, j] = 1
  end

  return NBI_matrix
end


# RCM-algorithm: Reverse Cuthill-McKeen algorithm (https://en.wikipedia.org/wiki/Cuthill%E2%80%93McKee_algorithm)
function mdoRCM(n, branchTupleVec)::Vector{Int}
  
  adjacency = zeros(Int, n, n)

  for (i, branch) in enumerate(branchTupleVec)
    from_bus, to_bus = branch
    adjacency[from_bus, to_bus] = 1
    adjacency[to_bus, from_bus] = 1
    #adjacency[to_bus, to_bus] = 1
    #adjacency[from_bus, from_bus] = 1
  end

  degrees = sum(adjacency, dims = 1)[:]

  # Initialize the elimination list
  eliminationList = collect(1:n)

  # Start the reverse Cuthill-McKee ordering
  currentLabel = 0
  order = zeros(Int, n)

  while length(eliminationList) > 0
    # Find the node with the minimum degree in the remaining graph
    minDegreeNode = argmin(degrees[eliminationList])
    currentLabel += 1
    currentNode = eliminationList[minDegreeNode]

    # Assign the label to the current node
    order[currentLabel] = currentNode

    # Update the degrees of the neighbors
    neighbors = findall(adjacency[currentNode, :] .== 1)
    for neighbor in neighbors
      degrees[neighbor] -= 1
    end

    # Remove the current node from the elimination list
    deleteat!(eliminationList, minDegreeNode)
  end

  return order
end
