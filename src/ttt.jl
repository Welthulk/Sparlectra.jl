# Author: Udo Schmitz (https://github.com/Welthulk)
# Date: 05.03.2024
# experimental arae

# outdated functions:
function setParallelBranches!(branches::Vector{Branch})
  branchTupleSet = Set{Tuple}()
  branchDict = Dict{Tuple{Integer,Integer},Vector{Branch}}()

  for b in branches
    tupple = (b.fromBus, b.toBus)

    if tupple in branchTupleSet
      existing_branches = branchDict[tupple]
      push!(existing_branches, b)
    else
      branchDict[tupple] = [b]
      push!(branchTupleSet, tupple)
    end
  end

  for (k, b_vec) in branchDict
    if length(b_vec) > 1
      sum_b_pu = 0.0
      sum_g_pu = 0.0
      sum_z = 0.0

      for b in b_vec
        b.isParallel = true
        b.skipYBus = true
        if b.status == 1          
          sum_b_pu += b.b_pu
          sum_g_pu += b.g_pu
          sum_z += (b.r_pu - b.x_pu * im) / (b.r_pu^2 + b.x_pu^2)
        end
      end

      z_total = 1.0 / sum_z
      r_pu = real(z_total)
      x_pu = imag(z_total)

      last_b = b_vec[end]
      for (i,b) in enumerate(b_vec)
        if b.status == 1
          adjP = AdjElecParams(r_pu = r_pu, x_pu = x_pu, b_pu = sum_b_pu, g_pu = sum_g_pu)
          setAdjElecParam!(adjP, b)
          @debug "branch (2): $(b)"
          if i==1
            @debug "last parallel element (3): $(b)"
            b.skipYBus = false
          end
        end
      end
    end
  end
end
