using Test
using DelimitedFiles
using MotifRecognition

cd(@__DIR__)
data = readdlm("confirmation")
notes = readdlm(path)
pitch = mod.(notes, 12)
intervals = pitch[2:end] .- pitch[1:end-1]
m = detect_motifs(intervals, 7, 1; iters = 700, tolerance = 0.7)
@test m[1].shape == [-1.0, -2.0, 10.0, -10.0, 2.0, 3.0, 5.0]
