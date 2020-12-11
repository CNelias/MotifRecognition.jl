using Combinatorics
using Distributions
using Random


"""
Returns dictionnary containing the probability of each symbol found in ts.
"""
function symbol_probabilities(ts)
    proba = values(counter(ts))
    proba = proba./sum(proba)
end

"""
    match_probability(l, d, p, iter)

Returns the probability of two random generated words of length 'w' from alphabet having probabilities 'p'
to match up to 'd' error. 't' is the size of the masked words.
Example:
Alphabet = {a,b,c} with probabilities [0.25, 0.25. 0.5], l = 8 and d = 3.
We sample from this distribution 'iter' words of length l = 8 and look how many of them match up to d = 3 errors.
"""
function match_probability(w, d, p, t)
    p_symbol_match = sum(p.*p)
    p_symbol_mismatch = 1 - p_symbol_match
    match_proba = 0
    for i in 1:d
        match_proba += binomial(w, i)*p_symbol_match^(w-i)*p_symbol_mismatch^i*(1-i/w)^t
    end
    return match_proba
end


"""
    significance_threshold(l, p, w, d, t)

Computes the expectation value of collision matrix entries for random words.
This is a conservative estimate that assumes all symbols are equiprobable
inputs (Int):
l : total length (size(S)[1])
p : probabilities of different symbols in alphabet.
w : motif size (window size)
d : # of allowed errors between motifs
t : length of projections after applying mask. Defaults to w - d.
returns (Int):
E : expecation value of collision entires.
"""
function significance_threshold(l, p, w, d, t = w - d)
    k = binomial(l, 2)
    average_match = k*match_probability(w, d, p, t)
    return floor(average_match)
end
