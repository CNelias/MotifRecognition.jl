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
"""
function get_masks(w, d)
    if d >= w
        error("d must be smaller than w")
    end
    return collect(combinations(collect(1:w), w - d))
end

"""
Returns a masked version of the S matrix.
"""
function apply_mask(mask, S)
    masked_S = Array{typeof(S[1, 1]),2}(undef, size(S)[1], length(mask))
    for (index, value) in enumerate(mask)
        masked_S[:, index] = S[:, value]
    end
    return masked_S
end

test = ["a","b","c","d","a","b","c","a","b","c","a","b","c","a","b","c"]
s = S_matrix(test, 4)
