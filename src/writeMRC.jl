
writeMRC(file::MRCImg, filename::AbstractString) = writeMRC(data(file), properties(file), filename)

function writeMRC(map::AbstractArray{T, 2}, head::Dict{Symbol, Any}, filename::AbstractString) where {T<:Real}
    map = permutedims(map, [2, 1])
    head = convertdictionary(head)
	f = writeMRCheader(map, head, filename)
    write(f, map)
    close(f)
    nothing
end

function writeMRC(map::AbstractArray{T, 3}, head::Dict{Symbol, Any}, filename::AbstractString) where {T<:Real}
    map = permutedims(map, [2, 1, 3])
    head = convertdictionary(head)
	f = writeMRCheader(map, head, filename)
    write(f, map)
    close(f)
    nothing
end

#=
    Convert the properties of the file from type Dict{Symbol, Any} to Dict{String, Any}.
=#
function convertdictionary(head::Dict{Symbol, Any})
    newHead = Dict{String, Any}()
    for (k, v) in head
        k = String(k)
        newHead[k] = v
    end
    return newHead
end

#=
    Determine the MODE as specified by the type.
=#
function writemode(meta::Dict{String, Any}, T::Type)
    # 0, 8-bit signed integer, -128 to 127.
    if T == Int8
        meta["MODE"] = 0
    # 1, 16-bit signed integer.
    elseif T == Int16
        meta["MODE"] = 1
    # 2, 32-bit signed real.
    elseif T == Float32
        meta["MODE"] = 2
    # 3, transform: complex 16-bit integer.
    elseif T == Complex{Int16}
        meta["MODE"] = 3
    # 4, transform: complex 32-bit reals.
    elseif T == ComplexF32
        meta["MODE"] = 4
    # 6, 16-bit unsigned integer.
    elseif T == UInt16
        meta["MODE"] = 6
    # 7, 32-bit signed integer.
    elseif T == Int32
        meta["MODE"] = 7
	else
		error("This data type is not supported!")
    end
    return meta
end

#=
    If we want to write a plane array to MRC file.
=#
function writeMRC(map::AbstractArray{T}, filename::AbstractString, flag::AbstractString) where {T<:Real}
    meta = MRCmeta()
    meta["NX"] = size(map, 1)
    meta["NY"] = size(map, 2)
    meta["NZ"] = size(map, 3)

    meta = writemode(meta, T)

    if flag == "image"        # single image
        meta["ISPG"] = 0
        meta["MZ"]   = 1
    elseif flag == "stack"    # image stack
        meta["ISPG"] = 0
        meta["MZ"]   = size(map, 3)
    elseif flag == "volume"    # single volume
        meta["ISPG"] = 1
        meta["MZ"]   = size(map, 3)
    else
        error("$flag is not supported for now!")
    end
    writeMRC(map, meta, filename)
end

#####################################################################################
#   Add padding to an Integer-array                                                 #
#####################################################################################
#   An array of type `IntX` may lose elements when reinterpreted to an array of     #
#   type `IntY` in cases cases where:                                               #
#       Y > X   AND     (Y modulus X) != 0                                          #
#                                                                                   #
#   e.g. Reinterpreting `Array{Int8,1}(collect(1:9)))` as Array{Int32} will result  #
#       in the loss of value `9`                                                    #
#                                                                                   #
#   This function resolves the issue by padding the array end of the array with 0's #
#   such that:                                                                      #
#       Y modulus X = 0                                                             #
#                                                                                   #
#  NOTE                                                                             #
#       This method is only configured for use with `Int32` return values.          #
#####################################################################################
function buffer(arr::Array{})
    len = length(arr)

    if isa(arr, Array{Int8,1})
        return reinterpret(Int32, vcat(arr[1:len], zeros(Int8, (4 - mod(len,4)))))
    elseif isa(arr, Array{UInt8,1})
        return reinterpret(Int32, vcat(arr[1:len], zeros(UInt8, (4 - mod(len,4)))))
    elseif isa(arr, Array{Int16,1})
        return reinterpret(Int32, vcat(arr[1:len], zeros(Int16, (2 - mod(len,2)))))
    elseif isa(arr, Array{UInt16,1})
        return reinterpret(Int32, vcat(arr[1:len], zeros(UInt16, (2 - mod(len,2)))))
    elseif isa(arr, Array{Int32,1})
        return arr
    elseif isa(arr, Array{UInt32,1})
        return reinterpret(Int32, arr)
    else
        error("Cannot buffer array of type \'", typeof(arr), "\' for conversion to Int32.")
    end
end

function stats(map::AbstractArray{T}) where {T<:Real}
    μ   = mean(map)
    σ   = std(map)
    Max = maximum(map)
    Min = minimum(map)
    return μ, σ, Max, Min
end

#=
    Record the layout of the MRC image.
=#
function writelayout(hdr::Vector{Int32}, meta::Dict{String, Any})
    hdr[1]  = meta["NX"]        # Number of image columns
    hdr[2]  = meta["NY"]        # Number of image rows
    hdr[3]  = meta["NZ"]        # Number of image sections
    hdr[4]  = meta["MODE"]      # Pixel-encoding signifier
    hdr[5]  = meta["NXSTART"]   # The number of the first column in the map (Default: 0)
    hdr[6]  = meta["NYSTART"]   # The number of the first row in the map
    hdr[7]  = meta["NZSTART"]   # The number of the first section in the map
    hdr[8]  = meta["MX"]        # The number of intervals along the X-axis
    hdr[9]  = meta["MY"]        # The number of intervals along the Y-axis
    hdr[10] = meta["MZ"]        # The number of intervals along the Z-axis
    return hdr
end

#=
    Record the dimensions, angles, axes, and density values of the MRC image.
=#
function writeaxes(hdr::Vector{Int32}, meta::Dict{String, Any}, Min::Float32, Max::Float32, μ::Float32)
    # Cell dimensions in ångströms
    hdr[11] = reinterpret(Int32, meta["CELL_A"])
    hdr[12] = reinterpret(Int32, meta["CELL_B"])
    hdr[13] = reinterpret(Int32, meta["CELL_C"])
    # Cell angles in degrees
    hdr[14] = reinterpret(Int32, meta["CELL_ALPHA"])
    hdr[15] = reinterpret(Int32, meta["CELL_BETA"])
    hdr[16] = reinterpret(Int32, meta["CELL_GAMMA"])
    # Axis alignment
    hdr[17] = reinterpret(Int32, meta["MAPC"])
    hdr[18] = reinterpret(Int32, meta["MAPR"])
    hdr[19] = reinterpret(Int32, meta["MAPS"])
    # Density values
    hdr[20:22] = reinterpret(Int32, Vector{Float32}([Min, Max, μ]))
    return hdr
end

#=
    Record the space, size of the extended header in bytes, type of the extended header, 
    and the version of the MRC format.
=#
function writespace(hdr::Vector{Int32}, meta::Dict{String, Any})
    hdr[23] = (meta["ISPG"] !== nothing) ? meta["ISPG"] : 0         # Space group (default value is 0)
    hdr[24] = meta["NSYMBT"]                                        # Length of the extended header (in bytes)
    hdr[25:26] = meta["EXTRA_SPACE1"]                               # Extra Space (misc data)
    hdr[27] = reinterpret(Int32, Vector{UInt8}(meta["EXTTYPE"]))[1] # Convert the 4-character extended data identifier (EXTTYPE) from String to Int32
    hdr[28] = meta["NVERSION"]                                      # The version of the MRC format used for data encoding
    hdr[29:49] = meta["EXTRA_SPACE2"]                               # Extra Space (misc data)
    return hdr
end

#=
    Record the phase origin (pixels) or origin of subvolume (A).
=#
function writeorigin(hdr::Vector{Int32}, meta::Dict{String, Any})
    hdr[50] = reinterpret(Int32, meta["ORIGIN_X"])
    hdr[51] = reinterpret(Int32, meta["ORIGIN_Y"])
    hdr[52] = reinterpret(Int32, meta["ORIGIN_Z"])
    return hdr
end

#= 
    Record the endianness of the machine. This helps to ensure that the image data 
    is read correctly on other machines.
=#
function writeendian(hdr::Vector{Int32})
    machineLE = (ENDIAN_BOM == 0x04030201) # True for Little-Endian machine
    if machineLE
        hdr[53] = reinterpret(Int32, b"MAP ")[1]
        hdr[54] = reinterpret(Int32, UInt8[68, 65, 0, 0])[1]
    else
        hdr[53] = reinterpret(Int32, b" PAM")[1]
        hdr[54] = reinterpret(Int32, UInt8[0, 0, 65, 68])[1]
    end
    return hdr
end

#=
    Record the label to the output header.
=#           
function writelabel(hdr::Vector{Int32}, meta::Dict{String, Any})
    # Labels are stored as an array of 10 80-character strings.
    tmp = zeros(UInt8, 10*80)
    for i in 1:length(meta["LABEL"])
        tmp[i] = UInt8(meta["LABEL"][i])
    end
    hdr[57:end] = reinterpret(Int32, tmp)
    return hdr
end

function dimensionerror(map::AbstractArray{T}, meta::Dict{String, Any}) where {T<:Real}
    if ndims(map) == 2
        size(map) == (meta["NX"], meta["NY"]) || error("Map dimensions do not match that of the input MRC file. This is not supported at the time.")
    else
        size(map) == (meta["NX"], meta["NY"], meta["NZ"]) || error("Map dimensions do not match that of the input MRC file. This is not supported at the time.")
    end
end

function writeMRCheader(map::AbstractArray{T}, meta::Dict{String, Any}, filename::AbstractString) where {T<:Real}
    # Use an array of Int32 types for the output header.
    # This is simply to allocate 32-bits for each metadata value,
    # actual data components will be of varying types.
    hdr = fill(Int32(0),256)
    
    dimensionerror(map, meta)

    # Get statistics
    μ, σ, Max, Min = stats(map)
    meta["NX"] = size(map, 1)
    meta["NY"] = size(map, 2)
    meta["NZ"] = size(map, 3)

    meta = writemode(meta, T)
    
    # Write header
    hdr = writelayout(hdr, meta)
    hdr = writeaxes(hdr, meta, Min, Max, μ)
    hdr = writespace(hdr, meta)
    hdr = writeorigin(hdr, meta)
    hdr = writeendian(hdr)
    hdr[55] = reinterpret(Int32, Float32(σ))    # RMS: rms deviation of map from mean density.
    hdr[56] = meta["NLABL"]                     # NLABL: The number of labels being used
    hdr = writelabel(hdr, meta)
    
    # If applicable, write extended header
    if meta["NSYMBT"] > 0
        hdr = vcat(hdr[1:end], buffer(meta["EXTHEAD"]))
    end

    handle=open(filename, "w")
    write(handle, hdr)
    return handle
end
