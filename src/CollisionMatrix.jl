using Combinatorics
using DataStructures
using Random

"""
    S_matrix(ts, window)

Returns S matrix containing all segments of length 'window' in 'ts'.
Done with a sliding window.
"""
function S_matrix(ts, window)
    S = Array{typeof(ts[1]),2}(undef, length(ts) - window, window)
    for index in 1:length(ts) - window
        S[index,:] = ts[index:index+window-1]
    end
    return S
end

"""
Returns all possible masks of length d among the possible w positions.
positive masks: the indexes represent the columns that are KEPT.
w : motif length
t : length of masked motifs
"""
function get_masks(w, t)
    if t >= w
        error("mask length must be smaller than w")
    end
    return collect(combinations(collect(1:w), t))
end

"""
Returns a masked version of the S matrix.
"""
function apply_mask(S, mask)
    masked_S = Array{typeof(S[1, 1]),2}(undef, size(S)[1], length(mask))
    for (index, value) in enumerate(mask)
        masked_S[:, index] = S[:, value]
    end
    return masked_S
end

"""
    update_c_matrix!(collision_matrix, masked_S)

Updates the collision matrix by adding 1 in the row where repetitions of a motif are found.
The row are the index of the repetitions, the column represent the first found motif of this type.
"""
function update_c_matrix!(collision_matrix, masked_S)
    masked_S_list = [masked_S[index,:] for index in 1:size(masked_S)[1]]
    count = counter(masked_S_list)
    for motif in keys(count)
        if count[motif] > 1
            positions = findall(x -> x == motif, masked_S_list)
            for p in positions
                collision_matrix[p, positions[1]] += 1
            end
        end
    end
end



"""
    collision_matrix(ts, w, d)

Constructs and returns the collision matrix of time-series 'ts'.
inputs (Int):
ts : input time-series
w : motif size (window size)
d : # of allowed errors between motifs
t : length of projections after applying mask. Defaults to w - d.
returns (Int):
C : collision matrix
"""
function collision_matrix(ts, w, d, t = w - d)
    C = zeros(length(ts) - w, length(ts) - w)
    S = S_matrix(ts, w)
    masked_S = Array{typeof(S[1, 1]),2}(undef, size(S)[1], t)
    for mask in get_masks(w, t)
        masked_S = apply_mask(S, mask)
        update_c_matrix!(C, masked_S)
    end
    return C
end

test = ["a","b","c","d","a","b","c","a","b","c","a","b","c","a","b","c"]
c = collision_matrix(test, 4, 2)
