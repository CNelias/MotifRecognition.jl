using Test
using DelimitedFiles
using MotifRecognition

cd(@__DIR__)
#input containing a single timeseries
data = readdlm("confirmation")
pitch = mod.(data, 12)
intervals = pitch[2:end] .- pitch[1:end-1]
m = detect_motifs(intervals, 7, 1; iters = 700, tolerance = 0.7)
consensus_shape = m[1].shape
similar_motifs = find_motifs(intervals, consensus_shape, 1)
@test m[1].shape == [-1.0, -2.0, 10.0, -10.0, 2.0, 3.0, 5.0]
@test similar_motifs.instances[1] ==  [-1.0, -2.0, 10.0, -10.0, 2.0, 3.0, 5.0]

#input including several timeseries
data = [intervals[1:540], intervals[540:1080], intervals[1080:1620]]
m = detect_motifs(intervals, 7, 1; iters = 700, tolerance = 0.7)

