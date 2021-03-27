include("NPSuperlattice.jl")
using Revise

using .NPSuperlattice
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
    parameters["np_diameter_std"] = parse(Float64, readline()) #Currently a percentage
    print("NP location disorder (nm): ")
    parameters["np_location_std"] = parse(Float64, readline())
    print("Number of samples: ")
    parameters["nsamples"] = parse(Int64, readline())

    parameters["n_nps"] =
        parameters["x_num"] * parameters["y_num"] * parameters["z_num"]

    return parameters
end

#create an empty starting distribution of necks
function createdistributions(parameters)
    npx = parameters["x_num"]
    npz = parameters["z_num"]
    npy = parameters["y_num"]
    nptotal = parameters["n_nps"]

    yneighbor = npx * npz
    xneighbor = npz
    nbor = 1

    adjacency_array = zeros(Int64, nptotal, nptotal)
    stored_edgearrays = Array{Any,1}[]
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

    for i = 0:(parameters["nsamples"]-1)
        push!(stored_edgearrays, edgearray)
    end

    return stored_edgearrays
end

function fillnecktables(
    i::Int64,
    nbor::Int64,
    adjacency_array::Array{Int64,2},
    edgearray::Array{Any,1},
)
    push!(edgearray, [(i, nbor), 0])
end

function createnpstructures(
    parameters::Dict{String,Any},
    folder::String,
)
    #pairneckdict = readinneckfile();
    #npdict, necknumberdict = readinnpfile();
    #fillreferencearray(false, pairneckdict, necknumberdict)
    #parameters = getparameters()
    samples = Array{Any,1}[]
    neckinginfo = createdistributions(parameters)
    neck_linear_regression = [0.0, 0.0,0.0]
    for i in neckinginfo
        sample = NPSuperlattice.createtricliniclattice(
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
    for d in 1:25
        disorder = 0.01*d
        parameters["np_diameter_std"] = disorder
        np_extension = "450_3.5nm_$(disorder)D/"
        reference_folder = "c:/Users/dunru/GitHub/HiNTS/data/"
        currentpath = pwd()
        cd(reference_folder)
        mkpath(np_extension)
        createnpstructures(parameters, np_extension)
        cd(currentpath)
    end
end

runsimulation()
