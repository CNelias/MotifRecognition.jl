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


test = ["a","b","c","d","a","b","c","a","b","c","a","b","c","a","b","c"]
s = S_matrix(test, 4)
