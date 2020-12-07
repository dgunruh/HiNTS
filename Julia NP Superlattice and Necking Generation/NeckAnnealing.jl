include("NPSuperlattice.jl")
using Revise

using .NPSuperlattice
using CSV
using KernelDensity
using Profile
using Distributions
using Random

global referencefreqarray = zeros(Int64, 7)

#Get the user defined input parameters
function getparameters()
    println()
    parameters = Dict{String,Any}()
    print("Number of NPs in the x direction: ")
    parameters["x_num"] = parse(Int64, readline())
    print("Number of NPs in the y direction: ")
    parameters["y_num"] = parse(Int64, readline())
    print("Number of NPs in the z direction: ")
    parameters["z_num"] = parse(Int64, readline())
    print("Expected average NP diameter (nm): ")
    parameters["np_diameter"] = parse(Float64, readline())
    print("NP-NP center-to-center spacing (nm): ")
    parameters["lattice_constant"] = parse(Float64, readline())
    print("NP diameter size disorder (nm): ")
    parameters["np_diameter_std"] = parse(Float64, readline())
    print("NP location disorder (nm): ")
    parameters["np_location_std"] = parse(Float64, readline())
    print("Number of samples: ")
    parameters["nsamples"] = parse(Int64, readline())

    parameters["n_nps"] =
        parameters["x_num"] * parameters["y_num"] * parameters["z_num"]

    return parameters
end

function updatereferencearray!(
    reference::Array{Int64,1},
    number_nanoparticles::Int64,
)

    #First create frequency histogram
    #In Xiaolei's case:
    #Grain 1: Gaussian with mean 3.85 and std 1.39
    #Grain 2: Gaussian with mean 4.1 and std 1.24
    #Grain 3: Gaussian with mean 3.55 and std 1.32
    print("Mean of the NP neck degree distribution: ")
    neck_degree_mean = parse(Float64, readline())
    print("Standard deviation of the NP neck degree distribution: ")
    neck_degree_std = parse(Float64, readline())

    gaus = Normal(neck_degree_mean, neck_degree_std)

    #We need the gaussian to be scaled so that the sum of all columns adds up to
    #the number of NPs
    neck_degrees = [0.0:6.0...]
    unscaled_frequencies = pdf.(gaus, neck_degrees)
    scale_factor = number_nanoparticles / sum(unscaled_frequencies)
    reference[:] = round.(Int64, scale_factor * unscaled_frequencies)
end

function createwidthregression()::Array{Float64}
    #Now create the linear regression which relates neck width
    #to NP-NP separation
    #Xiaolei's numbers: rms = 0.55, slope = -.17, intercept = 6.61
    #Here we will assume the same rms and slope, but will assume that the intercept
    #is the NP-NP spacing, equal to diameter + 2*(average ligand length)
    neckLinearRegression = [-0.17, 6.7, 0.55]
    return neckLinearRegression
end

"""
Create a distribution of neck widths which will be used to assign
neck widths to NP-NP pairs later.
This distribution can be drawn from data, but here we will assume
a Gaussian distribution of neck widths. Xiaolei found no correlation
between neck width and NP diameter, so random assignment of neck widths
is fine for our purposes.

Inputs:
-----------
a:  Float64
    The lattice constant of the NP-NP array.
"""
function createneckwidthdistribution(a::Float64)
    #Here we will assume that the neck width distribution is gaussian
    #Xiaolei observed two types: One heavily skewed to the right, and
    #one that was fat and centered
    print("Mean of the neck width distribution: ")
    μ = parse(Float64, readline())
    print("Standard deviation of the NP neck width distribution: ")
    σ = parse(Float64, readline())
    print("Number of samples to draw from: ")
    n = parse(Int8, readline())

    #Truncate the gaussian at the lattice spacing, and an assumed min size
    gaus = truncated(Normal(μ, σ), 1.0, a)
    neck_widths = rand(gaus, n)

    return neck_widths
end

#given a cumulative distribution, return a neck width probabilistically
function returnneckwidth(data_grid, cum_data)
    r = rand()
    for (index, value) in enumerate(cum_data)
        if value >= r
            return data_grid[index]
        end
    end
end

#calculate the monte-carlo energy of a frequency array compared to the reference.
#If it is the reference, the energy is zero
function energy(
    target_distribution::Array{Float64},
    current_distribution::Array{Float64},
)
    difference = target_distribution - current_distribution
    difference_squared = difference .^ 2
    return sqrt(sum(difference_squared))
end

#normalize an array. Generic function
norm(a::Array{Float64}) = a / sum(a)
norm(a::Array{Int64}) = a / sum(a)

#create an artifical starting distribution of necks
#the nanoparticles are distributed in a cubical lattice with dimensions specified by the user
#approximately 50% of possible necks are filled
function createdistribution(parameters)
    npx = parameters["x_num"]
    npz = parameters["z_num"]
    npy = parameters["y_num"]
    nptotal = parameters["n_nps"]

    yneighbor = npx * npz
    xneighbor = npz
    nbor = 1

    adjacency_array = zeros(Int64, nptotal, nptotal)
    edgearray = []
    for i = 1:nptotal
        rightz = false
        topx = false
        topy = false

        #go through all not possible connections
        #right z edge
        if i % xneighbor == 0
            rightz = true
        end
        #top x edge
        if i % (yneighbor) >= (yneighbor - xneighbor + 1) ||
           i % (yneighbor) == 0
            topx = true
        end
        #top y edge
        if i >= (nptotal - yneighbor + 1)
            topy = true
        end

        #add all connections
        #First add connection to the right-z neighbor. Include periodic bounds
        nbor = i + 1
        if rightz
            nbor = i + 1 - xneighbor
        end
        fillnecktables(i, nbor, adjacency_array, edgearray)

        #add connection to the above-x neighbor. Include periodic bounds
        nbor = i + xneighbor
        if topx
            nbor = i + xneighbor - yneighbor
        end
        fillnecktables(i, nbor, adjacency_array, edgearray)

        #add connection to the above-y neighbor. No periodic bounds
        nbor = i + yneighbor
        if !topy
            #nbor = i - (npy - 1) * yneighbor
            fillnecktables(i, nbor, adjacency_array, edgearray)
        end
    end

    necknumbers = sum(adjacency_array, dims = 2)
    thisfrequencylist = zeros(7)

    for k in necknumbers
        necknum = trunc(Int, k)
        thisfrequencylist[necknum+1] += 1
    end

    return thisfrequencylist, edgearray, adjacency_array
end

function fillnecktables(
    i::Int64,
    nbor::Int64,
    adjacency_array::Array{Int64,2},
    edgearray::Array{Any,1},
)
    add = rand()
    if add <= 0.5
        adjacency_array[i, nbor] = 1
        adjacency_array[nbor, i] = 1
        push!(edgearray, [(i, nbor), 1])
    else
        push!(edgearray, [(i, nbor), 0])
    end
end

#running the actual simulated annealing steps
function montecarlosteps(
    reffreqarray,
    thisfreqarray,
    edgearray,
    adjacencyarray,
    minbeta = 1000.0,
    maxbeta = 10000.0,
    betastep = 100.0,
    mcsteps = 1000,
)
    totaledges = length(edgearray)

    #create normalized reference array and current array
    normed_reffreqarray = norm(reffreqarray)
    normed_freqarray = norm(thisfreqarray)
    totalfreq = sum(thisfreqarray)

    #set current energy
    current_energy = energy(normed_reffreqarray, normed_freqarray)

    #set containers for lowest energy state
    #lowest_energy = current_energy
    #lowest_energy_adjacencyarray = similar(adjacencyarray)
    #lowest_energy_edgearray = similar(edgearray)
    proposed_freqarray = Vector{Float64}(undef, 7)
    degreetable = sum(adjacencyarray, dims = 1)
    change = (-1, 1)
    numbersteps =
        trunc(Int, (maxbeta - minbeta + betastep) * mcsteps / betastep)
    edgerandoms = rand(1:totaledges, numbersteps)
    montecarlorandoms = rand(Float64, numbersteps)
    riterator = 1

    for beta = minbeta:betastep:maxbeta
        for i = 1:mcsteps

            #select edge
            #edgenum = rand(1:totaledges)
            edgenum = edgerandoms[riterator]

            edgestate = edgearray[edgenum][2]
            node1 = edgearray[edgenum][1][1]
            node2 = edgearray[edgenum][1][2]
            node1degree = degreetable[node1]
            node2degree = degreetable[node2]


            #             proposed_fnd = node1degree - change[edgestate + 1]
            #             proposed_snd = node2degree - change[edgestate + 1]
            #             proposed_edgestate = 1 - edgestate

            if edgestate == 0
                proposed_fnd = node1degree + 1
                proposed_snd = node2degree + 1
                proposed_edgestate = 1
            else
                proposed_fnd = node1degree - 1
                proposed_snd = node2degree - 1
                proposed_edgestate = 0
            end

            #look at change in frequency array
            @inbounds thisfreqarray[node1degree+1] -= 1
            @inbounds thisfreqarray[node2degree+1] -= 1
            @inbounds thisfreqarray[proposed_fnd+1] += 1
            @inbounds thisfreqarray[proposed_snd+1] += 1
            normed_freqarray = norm(thisfreqarray)

            #calculate resulting change in energy
            new_energy = energy(normed_reffreqarray, normed_freqarray)
            delta_energy = new_energy - current_energy

            #calculate probability of flip
            probabilityflip = exp(-beta * delta_energy)

            #make flip if possible
            r = montecarlorandoms[riterator]
            if r <= probabilityflip
                #                 adjacencyarray[edgearray[edgenum][1][1],edgearray[edgenum][1][2]] = proposed_edgestate
                #                 adjacencyarray[edgearray[edgenum][1][2],edgearray[edgenum][1][1]] = proposed_edgestate
                #                 degreetable[edgearray[edgenum][1][1]] += change[proposed_edgestate + 1]
                #                 degreetable[edgearray[edgenum][1][2]] += change[proposed_edgestate + 1]
                adjacencyarray[node1, node2] = proposed_edgestate
                adjacencyarray[node2, node1] = proposed_edgestate
                degreetable[node1] += change[proposed_edgestate+1]
                degreetable[node2] += change[proposed_edgestate+1]
                edgearray[edgenum][2] = proposed_edgestate
                current_energy = new_energy

                #                 if current_energy < lowest_energy
                #                    lowest_energy = current_energy
                #                    lowest_energy_adjacencyarray .= adjacencyarray
                #                    lowest_energy_edgearray .= edgearray
                #                 end
            else
                @inbounds thisfreqarray[node1degree+1] += 1
                @inbounds thisfreqarray[node2degree+1] += 1
                @inbounds thisfreqarray[proposed_fnd+1] -= 1
                @inbounds thisfreqarray[proposed_snd+1] -= 1
            end
            riterator += 1
        end
    end

    #return lowest_energy_edgearray, lowest_energy_adjacencyarray
    return thisfreqarray, edgearray, adjacencyarray
end

#run the simulated annealing steps to create a distribution of necks
#use that to create a frequency array
#then assign neck widths to each neck
function simulatedannealing(
    referencefreqarray,
    thisfreqarray,
    edgearray,
    adjacencyarray,
    neckwidtharray,
)

    finalfreqarray, lowest_energy_edgearray, lowest_energy_adjacencyarray =
        montecarlosteps(
            referencefreqarray,
            thisfreqarray,
            edgearray,
            adjacencyarray,
        )

    k = kde(neckwidtharray)
    data_grid = k.x
    kernel_data = k.density
    sum_kernel_data = sum(kernel_data)
    cum_kernel_data = cumsum(kernel_data / sum_kernel_data)

    #assign neck widths to each edge
    for (n, i) in enumerate(lowest_energy_edgearray)
        pair = i[1]
        neckexists = (i[2] == 1)

        if neckexists
            neckwidth = returnneckwidth(data_grid, cum_kernel_data)
            lowest_energy_edgearray[n][2] = neckwidth
        end
    end

    return finalfreqarray, lowest_energy_edgearray, lowest_energy_adjacencyarray
end

function createneckingstructures(parameters::Dict{String,Any}, folder::String)
    #pairneckdict = readinneckfile();
    #npdict, necknumberdict = readinnpfile();
    #fillreferencearray(false, pairneckdict, necknumberdict)
    #parameters = getparameters()
    #referencefreqarray = returnreferencearray(parameters["n_nps"])
    updatereferencearray!(referencefreqarray, parameters["n_nps"])
    neck_linear_regression = createwidthregression()
    neckwidtharray = createneckwidthdistribution(parameters["lattice_constant"])

    println("Reference neck degree array:")
    println(referencefreqarray)

    stored_edgearrays = Array{Any,1}[]

    t0 = time_ns()
    for i = 0:(parameters["nsamples"]-1)
        thisfreqarray, edgearray, adjacencyarray =
            createdistribution(parameters)
        finalfreqarray, lowest_energy_edgearray, lowest_energy_adjacencyarray =
            simulatedannealing(
                referencefreqarray,
                thisfreqarray,
                edgearray,
                adjacencyarray,
                neckwidtharray,
            )
        println(finalfreqarray)

        filepath = "$(folder)neckSample$i.inp"
        open(filepath, "w") do io
            #create(io)
            write(io, "id1 id2 type necking_width\n")

            for j = 1:length(lowest_energy_edgearray)
                id1 = lowest_energy_edgearray[j][1][1] - 1
                id2 = lowest_energy_edgearray[j][1][2] - 1
                neck_width = lowest_energy_edgearray[j][2]

                if neck_width > 0
                    write(io, "$id1 $id2 1 $neck_width\n")
                end
            end
        end

        push!(stored_edgearrays, lowest_energy_edgearray)
    end
    tf = time_ns()
    println((tf - t0) / 1.0e9)

    return stored_edgearrays
end

function createnpstructures(
    parameters::Dict{String,Any},
    neckinginfo::Array{Array{Any,1},1},
    folder::String,
)
    #pairneckdict = readinneckfile();
    #npdict, necknumberdict = readinnpfile();
    #fillreferencearray(false, pairneckdict, necknumberdict)
    #parameters = getparameters()
    samples = Array{Any,1}[]
    neck_linear_regression = createwidthregression()
    for i in neckinginfo
        sample = createtricliniclattice(
            parameters,
            i,
            neck_linear_regression,
        )
        np_array = sample.np_array
        flipped_np_array = sample.flipped_np_array
        cellsize_array = sample.cellsize_array

        push!(samples, [np_array, cellsize_array])
        push!(samples, [flipped_np_array, cellsize_array])
    end

    t0 = time_ns()
    for i = 1:(parameters["nsamples"]*2)
        filepath = "$(folder)npSample$(i-1).inp"
        open(filepath, "w") do io
            write(io, "BANDS\n")
            write(io, "8 8\n")
            write(io, "NANOPARTICLES\n")
            write(io, "cell(1), $(samples[i][2][1])\n")#X  center-center distance leftmost-rightmost + desired diameter
            write(io, "cell(2), $(samples[i][2][2])\n")#Y  Defines height (layers)
            write(io, "cell(3), $(samples[i][2][3])\n")#Z  Transport direction
            write(io, "id type layer x y z diameter\n")

            for j = 1:length(samples[i][1])
                np_type = 1
                np_layer = 1
                write(
                    io,
                    "$(j-1) $np_type $np_layer $(samples[i][1][j][1]) $(samples[i][1][j][2]) $(samples[i][1][j][3]) $(samples[i][1][j][4])\n",
                )
            end
        end
    end
    tf = time_ns()
    println((tf - t0) / 1.0e9)
end

function runsimulation()
    parameters = getparameters()
    x_num = parameters["x_num"]
    y_num = parameters["y_num"]
    z_num = parameters["z_num"]
    neck_extension = "neckDistributions/$(x_num)x$(y_num)x$(z_num)/"
    np_extension = "npSolids/$(x_num)x$(y_num)x$(z_num)/"
    reference_folder = "c:/Users/Davis/Research/Superlattice Project/"
    currentpath = pwd()
    cd(reference_folder)
    mkpath(neck_extension)
    stored_neckstructures = createneckingstructures(parameters, neck_extension)
    mkpath(np_extension)
    createnpstructures(parameters, stored_neckstructures, np_extension)
    cd(currentpath)
end

runsimulation()
