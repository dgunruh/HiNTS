module NPSuperlattice

using Revise
using Distributions
using Random
using LinearAlgebra

export superlattice, createtricliniclattice

mutable struct superlattice
    structure::Tuple{Vararg{Int64,3}}
    angles::Tuple{Vararg{Float64,3}}
    basisvectors::Tuple{Vararg{Array{Float64},3}}
    neckmap::Array{Array{Any,1},1}
    neckcenter_regression::Array{Float64}
    np_array::Array{Any,1}
    flipped_np_array::Array{Any,1}
    cellsize_array::Array{Float64}
end

#normalize an array. Generic function
normalize(a::Array{Float64}) = a / sum(a)
normalize(a::Array{Int64}) = a / sum(a)



function createtricliniclattice(
    parameters::Dict{String,Any},
    neckmap::Array{Any,1},
    neckcenter_regression::Array{Float64},
)::superlattice
    x_num = parameters["x_num"]
    y_num = parameters["y_num"]
    z_num = parameters["z_num"]
    nptotal = parameters["n_nps"]
    diameter_μ = parameters["np_diameter"]
    diameter_σ = parameters["np_diameter_std"]
    ligandlength = (parameters["lattice_constant"] - diameter_μ) / 2

    #The basis lengths will be ordered as x, y, z
    start = [
        parameters["lattice_constant"] / 2,
        parameters["lattice_constant"] / 2,
        parameters["lattice_constant"] / 2,
    ]
    lengths = 2 * start
    angles, vectors = setlatticeinformation(lengths)
    cellx = x_num * vectors[2][1]
    celly = y_num * vectors[3][2]
    cellz = z_num * vectors[1][3]
    cellarray = [cellx, celly, cellz]

    #To prevent nanoparticles from overlapping vertically once shifted by necks,
    #cap the initial radius at half of the lattice spacing in the y-direction
    diameter_max = vectors[3][2]
    #diameter_min = diameter_μ - 2 * ligandlength
    diameter_min = 3*diameter_μ/4
    diameter_gaussian =
        truncated(Normal(diameter_μ, diameter_σ), diameter_min, diameter_max)

    sl = superlattice(
        (x_num, y_num, z_num),
        angles,
        vectors,
        neckmap,
        neckcenter_regression,
        [],
        [],
        cellarray,
    )

    np_id = 1
    np_diameter = diameter_μ

    base_x = start[1]
    base_z = start[3]
    base_z_y = start[3]
    np_loc = start

    for y = 1:y_num
        for x = 1:x_num
            for z = 1:z_num
                np_diameter = rand(diameter_gaussian)
                appendNP!(np_loc, np_diameter, np_id, sl)
                np_id += 1

                np_loc[1] += vectors[1][1]
                np_loc[2] += vectors[1][2]
                np_loc[3] += vectors[1][3]

            end

            base_z += vectors[2][3]

            np_loc[3] = base_z
            np_loc[1] += vectors[2][1]
            np_loc[2] += vectors[2][2]
        end

        base_x += vectors[3][1]
        base_z_y += vectors[3][3]
        base_z = base_z_y

        np_loc[1] = base_x
        np_loc[3] = base_z

        np_loc[2] += vectors[3][2]
    end

    #Now iterate through each array in a random walk and allow NP upper bounds
    #to extend in size
    prime = 37
    offset = 0
    for index in 1:nptotal
        i = (index * prime + offset)%nptotal + 1
        xlower = i - z_num
        if (i-1)%(x_num*z_num) < z_num
            xlower = i + (x_num - 1) * z_num
        end
        xupper = i + z_num
        if (i-1)%(x_num*z_num) >= (x_num - 1) * z_num
            xupper = i - (x_num - 1) * z_num
        end
        ylower = i - x_num*z_num
        yupper = i + x_num*z_num
        zlower = i - 1
        if i%z_num == 1
            zlower = i + z_num - 1
        end
        zupper = i + 1
        if i %z_num == 0
            zupper = i - (z_num - 1)
        end
        currentdiam = sl.np_array[i][4]

        lowxdiam = sl.np_array[xlower][4]
        upxdiam = sl.np_array[xupper][4]
        if !(ylower < 1)
            lowydiam = sl.np_array[ylower][4]
        else
            lowydiam = 0.0
        end
        if !(yupper > nptotal)
            upydiam = sl.np_array[yupper][4]
        else
            upydiam = 0.0
        end
        lowzdiam = sl.np_array[zlower][4]
        upzdiam = sl.np_array[zupper][4]

        maxneighdiam = max(lowxdiam, upxdiam, lowydiam, upydiam, lowzdiam, upzdiam)
        alloweddiam = (vectors[3][2] - maxneighdiam/2.0) * 2
        diameter_gaussian =
            truncated(Normal(diameter_μ, diameter_σ), diameter_min, alloweddiam)
        np_diameter = rand(diameter_gaussian)
        sl.np_array[i][4] = np_diameter
        sl.flipped_np_array[i][4] = np_diameter
    end

    return sl
end

function setlatticeinformation(basislengths::Array{Float64})
    #Triclinic angles. You can only have periodic necking in a cubic lattice,
    #or a triclinic lattice with a carefully curated number of NPs
    α = π * 103 / 180.0
    β = π * 95 / 180.0
    γ = π * 96 / 180.0

    #α = π * 90 / 180.0
    #β = π * 90 / 180.0
    #γ = π * 90 / 180.0

    a₁ = [0, 0, basislengths[3]]
    a₂ = [basislengths[1] * sin(γ), 0, basislengths[1] * cos(γ)]

    x₃ = basislengths[2] * (cos(β) - cos(α) * cos(γ)) / sin(γ)
    z₃ = basislengths[2] * cos(α)
    y₃ = sqrt(basislengths[2]^2 - x₃^2 - z₃^2)
    a₃ = [x₃, y₃, z₃]

    return (α, β, γ), (a₁, a₂, a₃)
end

function movedistance(
    neckwidth::Float64,
    sidelength::Float64,
    ncr::Array{Float64})::Float64
    return (sidelength - rand(truncated(
            Normal(neckwidth * ncr[1] + ncr[2], ncr[3]),
            sidelength * 3 / 4,
            sidelength,
        ),))/2
end

function movedistance(
    neckwidth::Int64,
    sidelength::Float64,
    ncr::Array{Float64})::Float64
    return 0.0
end

function appendNP!(
    location::Array{Float64},
    diameter::Float64,
    np_id::Int64,
    sl::superlattice,
)

    pairs = [i for i in sl.neckmap if np_id in i[1]]
    jitter = [0.0, 0.0, 0.0]
    (x_num, y_num, z_num) = sl.structure
    x_sidelength = norm(sl.basisvectors[1])
    y_sidelength = norm(sl.basisvectors[2])
    z_sidelength = norm(sl.basisvectors[3])
    # sidelength = x_sidelength
    # slope = sl.neckcenter_regression[1]
    # intercept = sl.neckcenter_regression[2]
    # rms = sl.neckcenter_regression[3]
    difference = 0.0

    for i in pairs
        if np_id == i[1][1]
            other_id = i[1][2]
        else
            other_id = i[1][1]
        end
        #draw predicted center to center distance from normal distribution centered on line
        #predicted_center_center = np.random.normal(neckWidth*slope + intercept, rMS)

        #z neighbors. Don't forget periodic boundary conditions!
        # if other_id == np_id - 1 ||
        #    other_id == np_id + 1 ||
        #    other_id == np_id + (z_num - 1) ||
        #    other_id == np_id - (z_num - 1)
        #     sidelength = z_sidelength
        # end
        #
        # #x neighbors
        # if other_id == np_id + z_num ||
        #    other_id == np_id - z_num ||
        #    other_id == np_id + z_num * (x_num - 1) ||
        #    other_id == np_id - z_num * (x_num - 1)
        #     sidelength = x_sidelength
        # end
        #
        # #y neighbors
        # if other_id == np_id + z_num * x_num ||
        #    other_id == np_id - z_num * x_num ||
        #    other_id == np_id + z_num * x_num * (y_num - 1) ||
        #    other_id == np_id - z_num * x_num * (y_num - 1)
        #     sidelength = y_sidelength
        # end
        # if (other_id == np_id - 1 && !(np_id%z_num == 1)) || (other_id == np_id + z_num && (np_id%z_num == 1))
        #     sidelength = z_sidelength
        # end
        #
        # #z-right neighbor
        # if (other_id == np_id + 1 && !(np_id%z_num == 0)) || (other_id == np_id - (z_num - 1) && (np_id%z_num == 0))
        #     sidelength = z_sidelength
        # end
        #
        # #x-bottom neighbor
        # if (other_id == np_id - z_num && !(np_id%(z_num*x_num + 1) <= z_num)) || (other_id == np_id + z_num*(x_num - 1) && (np_id%(z_num*x_num + 1) <= z_num))
        #     sidelength = x_sidelength
        # end
        #
        # #x-top neighbor
        # if (other_id == np_id + z_num && !(np_id%(z_num*x_num + 1) > (z_num*(x_num - 1)))) || (other_id == np_id + z_num && !(np_id%(z_num*x_num) > (z_num*(x_num - 1))))
        #     sidelength = x_sidelength
        # end
        #
        # #y-bottom neighbor
        # if other_id == np_id - z_num*z_num
        #     sidelength = y_sidelength
        # end
        #
        # #y-top neighbor
        # if other_id == np_id + x_num*z_num
        #     sidelength = y_sidelength
        # end

        #draw predicted center to center distance from normal distribution centered on line
        # if i[2] == 0
        #     predicted_center_center = sidelength
        # else
        #     predicted_center_center = rand(truncated(
        #         Normal(i[2] * slope + intercept, rms),
        #         sidelength * 3 / 4,
        #         sidelength,
        #     ),)
        # end

        #difference = (sidelength - predicted_center_center) / 2

        #z-left neighbor
        #if other_id == np_id - 1 && !(np_id%z_num == 0)
        if other_id == np_id - 1 || other_id == np_id + z_num - 1
            difference = movedistance(i[2], z_sidelength, sl.neckcenter_regression)
            jitter[3] -= difference
        end

        #z-right neighbor
        #if other_id == np_id + 1 && !(np_id%z_num == (z_num - 1))
        if other_id == np_id + 1 || other_id == np_id - (z_num - 1)
            difference = movedistance(i[2], z_sidelength, sl.neckcenter_regression)
            jitter[3] += difference
        end

        #x-bottom neighbor
        #if other_id == np_id - self.zNum && !(np_id%(self.zNum*self.xNum) < self.zNum)
        if other_id == np_id - z_num || other_id == np_id + z_num * (x_num - 1)
            difference = movedistance(i[2], x_sidelength, sl.neckcenter_regression)
            jitter[1] -= difference * sl.basisvectors[2][1] / x_sidelength
            jitter[3] -= difference * sl.basisvectors[2][3] / x_sidelength
        end

        #x-top neighbor
        #if other_id == (np_id + self.zNum) && !(np_id%(self.zNum*self.xNum) > (self.zNum*self.xNum - self.zNum))
        if other_id == np_id + z_num || other_id == np_id - z_num * (x_num - 1)
            difference = movedistance(i[2], x_sidelength, sl.neckcenter_regression)
            jitter[1] += difference * sl.basisvectors[2][1] / x_sidelength
            jitter[3] += difference * sl.basisvectors[2][3] / x_sidelength
        end

        #y-bottom neighbor
        if other_id == np_id - x_num * z_num ||
           other_id == np_id + x_num * z_num * (y_num - 1)
            difference = movedistance(i[2], y_sidelength, sl.neckcenter_regression)
            jitter .-= difference * normalize(sl.basisvectors[3])
        end

        #y-top neighbor
        if other_id == np_id + x_num * z_num ||
           other_id == np_id - x_num * z_num * (y_num - 1)
            difference = movedistance(i[2], y_sidelength, sl.neckcenter_regression)
            jitter .+= difference * normalize(sl.basisvectors[3])
        end
    end

    #append nanoparticle
    flipped_information = push!(location + jitter, diameter)
    flipped_information[3] += (sl.cellsize_array[3] - 2 * location[3] - 2*jitter[3])
    push!(sl.np_array, push!(location + jitter, diameter))
    push!(sl.flipped_np_array, flipped_information)
end
end
