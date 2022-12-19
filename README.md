# Motif recognition

 | **Appveyor** |
|:---------------------:|
|[![Build status](https://ci.appveyor.com/api/projects/status/lafds8w3gv04t39v?svg=true)](https://ci.appveyor.com/project/johncwok/motifrecognition-jl)

This module proposes a detection algorithm based on JEREMY BUHLER and MARTIN TOMPA's paper "*Finding Motifs Using Random Projections*" which is well adapted for categorical time-series.  It is now part of the package [CategoricalTimeSeries.jl](https://github.com/johncwok/CategoricalTimeSeries.jl)<br/>
The algorithm although very precise is not exact. Therefore, when you are done detecting potential motifs with the `detect_motifs` function, you can refine your results with `find_motifs` for an exact search.
<br/>

For a quick overview, have a look at the **examples** provided **at the end of this readme**.

The main functions return instances of a class called **pattern**:
- - -
**pattern — Class**
- - -
A class storing useful information about found motifs in a time-series. An array of `pattern` instances is returned when the searching algorithm is done running.
>**Attributes**:

- **shape** ([Array{Any,1}](https://docs.julialang.org/en/v1/base/arrays/)): Array containing the shape (or contour) of the first found repetition of the motif.
- **repetitions** ([Array{Array{Any,1},1}](https://docs.julialang.org/en/v1/base/arrays/)): all the different shapes from the motif's repetitions, they can vary a bit from one to the next.
- **positions** ([Array{Int,1}](https://docs.julialang.org/en/v1/base/arrays/)): the positions at which the different repetitions of the motif were found.

## Main functions

**detect_motifs**

```Julia
detect_motifs(ts, w, d, t = w - d; iters = 1000, tolerance = 0.95)
```
Detects all motifs of length 'w' occuring more often than chance, being identical to each other up to 'd' differences inside of imput time-series 'ts'. 'ts' can be a single vector containing one timeseries, a vector of vectors or a matrix containing several timeseries. 
Returns an array of `pattern`, inside of which the patterns are classified by how frequently they are observed. The first elements is therefore the most frequently observed motif, and so on.
> **Parameters**:

>>* **ts** ([Array{Any,1}, Array{Array{Any,1}, 1}, Array{Any, 2}](https://docs.julialang.org/en/v1/base/arrays/)): input time-series in which motifs are searched for.
>>* **w** ([Int](https://docs.julialang.org/en/v1/manual/integers-and-floating-point-numbers/)): length of motifs to look for.
>>* **d** ([Int](https://docs.julialang.org/en/v1/manual/integers-and-floating-point-numbers/)): allowed errors (differences) between motifs repetitions.
>>* **t = w - d** ([Int](https://docs.julialang.org/en/v1/manual/integers-and-floating-point-numbers/)): size of the masks to use for random projection in the detection (defaults to w - d).
>>* **iters = 1000** ([Int](https://docs.julialang.org/en/v1/manual/integers-and-floating-point-numbers/)): the numbers of iterations for the random projection process (defaults to 1000)
>>* **tolerance = 0.95** ([Float](https://docs.julialang.org/en/v1/manual/integers-and-floating-point-numbers/)): threshold of motif identification. If set to 1, only matrix entries that are strictly superior to the (probabilistic) threshold are taken into account. Defaults to 0.7, meaning that matrix entries need to be bigger than 0.7*threshold.

> **Returns** :
>>* **motifs** : list of `pattern` instances sorted by frequency of occurence. motifs[1] is therefore the most frequent motif, motifs[2] the second most observed and so on.


**find_motifs**

```Julia
find_motifs(ts, shape, d)
```
Given a motif of shape 'shape' (array{any,1}), looks for all the repetitions of it which differ only up to 'd' differences inside of the input time-series 'ts'.
Input:

>**Parameters**:
>>* **ts** ([Array{Any,1}](https://docs.julialang.org/en/v1/base/arrays/)) : time-series in which to look for motifs
>>* **shape** ([Array{Any,1}](https://docs.julialang.org/en/v1/base/arrays/)): shape (aray{any,1}) of the motif to look for.
>>* **d** ([Int](https://docs.julialang.org/en/v1/manual/integers-and-floating-point-numbers/)): allowed errors (differences) between motifs

>* **Returns** :

>>* **motif** : an instance of `pattern` containing the found repetition of the input 'shape'.

## Plotting
To help visualize results, two simple plotting functions are provided.

**plot_motif**

```Julia
plot_motif(m::pattern)
```
Plots all repetitions of an input `pattern` instance on top of each other to see how similar they are to each other.
> **Parameters**:

>>* **m** : Instance of the `pattern` class


**plot_motif**

```Julia
plot_motif(m::pattern, ts)
```
Plots all repetitions of an input `pattern` instance on top of the input time-series 'ts' to better visualize their repartition in time.
> **Parameters**:

>>* **m** : Instance of the `pattern` class
>>* **ts** ([Array{Any,1}](https://docs.julialang.org/en/v1/base/arrays/)): Input time-series


## Example
From Michael Brecker's improvisation over the piece ["confirmation"](https://github.com/johncwok/MotifRecognition.jl/tree/main/test), we extract a time-series of pitch intervals (difference from one note to the next).
A spectral envelope analysis reveals a peak at period 6~7, so we look for motifs of length 7 and allow for 1 error between them.
After detection, we visualize the most frequent motif:
```julia
using MotifRecognition
using DelimitedFiles
using Plots

data = readdlm(joinpath(pkgdir(MotifRecognition), "test", "confirmation"))
pitch = mod.(data, 12) #Removing octave position: not needed
intervals = pitch[2:end] .- pitch[1:end-1] #getting interval time-series.
m = detect_motifs(intervals, 7, 1; iters = 700, tolerance = 0.7)
plot(m[1]) #plotting most frequent motif
```

<img src=https://user-images.githubusercontent.com/34754896/104308882-9c2c9e80-54d1-11eb-8882-cc31b7b2af8b.PNG width = "500">

We notice that the motif `[-1, -2, 10, -10, 2, 3, 5]` seems to be the underlying (consensus) shape. In musical notation, this motif would look like this (written in C major):

<img src=https://user-images.githubusercontent.com/34754896/104315350-1ca3cd00-54db-11eb-864d-3a1da9d5efeb.PNG width = "500">

We do an exact search with 1 error allowed to check if our previous detection missed any repetitions, and plot the found motif on top of each other:
```julia
consensus_shape = [-1, -2, 10, -10, 2, 3, 5]
motif = find_motifs(intervals, consensus_shape, 1)
plot(motif)
```
<img src=https://user-images.githubusercontent.com/34754896/104308882-9c2c9e80-54d1-11eb-8882-cc31b7b2af8b.PNG width = "500">

Here, we obtain the same plot as before but this is not necessarily always the case. Knowing the consensus motif usually allows to find its repetitions more efficiently.

Now, we visualize the repetitions of the motif in the time-series:
```julia
plot(motif, intervals)
```
<img src=https://user-images.githubusercontent.com/34754896/104313663-a1411c00-54d8-11eb-9854-70bd5ed9ba2f.PNG width = "600">


### Citing
If you used this module in a scientific publication, please consider citing the package it came from:
```Bibtex
@article{nelias2021categoricaltimeseries,
  title={CategoricalTimeSeries. jl: A toolbox for categorical time-series analysis},
  author={Nelias, Corentin},
  journal={Journal of Open Source Software},
  volume={6},
  number={67},
  pages={3733},
  year={2021}
}
```
