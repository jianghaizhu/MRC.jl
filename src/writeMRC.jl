
writeMRC(file::MRCImg, filename::AbstractString) = writeMRC(data(file), properties(file), filename)

function writeMRC(map::AbstractArray{T, 2}, head::Dict{String, Any}, filename::AbstractString) where {T<:Real}
	map = permutedims(map, [2, 1])
	f = writeMRCheader(map, head, filename)
    write(f, map)
    close(f)
    nothing
end

function writeMRC(map::AbstractArray{T, 3}, head::Dict{String, Any}, filename::AbstractString) where {T<:Real}
	map = permutedims(map, [2, 1, 3])
	f = writeMRCheader(map, head, filename)
    write(f, map)
    close(f)
    nothing
end

# if we want to write a plane array to mrc file
function writeMRC(map::AbstractArray{T}, filename::AbstractString, flag::AbstractString) where {T<:Real}
    meta = MRCmeta()
    meta["NX"] = size(map, 1)
    meta["NY"] = size(map, 2)
    meta["NZ"] = size(map, 3)
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
    elseif isa(arr, Array{Unt8,1})
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

function writeMRCheader(map::AbstractArray{T}, meta::Dict{String, Any}, filename::AbstractString) where {T<:Real}
    machineLE = (ENDIAN_BOM == 0x04030201) # True for Little-Endian machine
    #= Use an array of Int32 types for the output header.
        This is simply to allocate 32-bits for each metadata value,
        actual data components will be of varying types.                 =#
    hdr = fill(Int32(0),256)

    size(map) == (meta["NX"], meta["NY"], meta["NZ"]) || error("Map dimensions do not match that of the input MRC file. This is not supported at the time.")

    # Get statistics
    μ   = mean(map)
    σ   = std(map)
    Max = maximum(map)
    Min = minimum(map)
    meta["NX"] = size(map, 1)
    meta["NY"] = size(map, 2)
    meta["NZ"] = size(map, 3)

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
    end

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
    #= Cell dimensions in ångströms: =#
    hdr[11] = reinterpret(Int32, meta["CELL_A"])
    hdr[12] = reinterpret(Int32, meta["CELL_B"])
    hdr[13] = reinterpret(Int32, meta["CELL_C"])
    #=  Cell angles in degrees: =#
    hdr[14] = reinterpret(Int32, meta["CELL_ALPHA"])
    hdr[15] = reinterpret(Int32, meta["CELL_BETA"])
    hdr[16] = reinterpret(Int32, meta["CELL_GAMMA"])
    #= Axis Alignment: =#
    hdr[17] = reinterpret(Int32, meta["MAPC"])
    hdr[18] = reinterpret(Int32, meta["MAPR"])
    hdr[19] = reinterpret(Int32, meta["MAPS"])
    #= Density Values: =#
    hdr[20:22] = reinterpret(Int32, Vector{Float32}([Min, Max, μ]))
    # Space group (default value is 0)
    hdr[23] = (meta["ISPG"] != nothing) ? meta["ISPG"] : 0
    # Length of the extended header (in bytes):
    hdr[24] = meta["NSYMBT"]
    # Extra Space (misc data)
    hdr[25:26] = meta["EXTRA_SPACE1"]
    # Convert the 4-character extended data identifier (EXTTYPE)
    #   from String to Int32
    hdr[27] = reinterpret(Int32, Vector{UInt8}(meta["EXTTYPE"]))[1]
    # The version of the MRC format used for data encoding
    hdr[28] = meta["NVERSION"]
    # Extra Space (misc data)
    hdr[29:49] = meta["EXTRA_SPACE2"]
    hdr[50] = reinterpret(Int32, meta["ORIGIN_X"])
    hdr[51] = reinterpret(Int32, meta["ORIGIN_Y"])
    hdr[52] = reinterpret(Int32, meta["ORIGIN_Z"])
    # Record the endianness of the machine
    #   This helps to ensure that the image data is read correctly on other machines
    if machineLE
        hdr[53] = reinterpret(Int32, b"MAP ")[1]
        hdr[54] = reinterpret(Int32, UInt8[68, 65, 0, 0])[1]
    else
        hdr[53] = reinterpret(Int32, b" PAM")[1]
        hdr[54] = reinterpret(Int32, UInt8[0, 0, 65, 68])[1]
    end
    hdr[55] = reinterpret(Int32, Float32(σ))

    # NLABL:    The number of labels being used
    # LABEL:    Label values
    #           Labels are stored as an array of 10 80-character strings.

    hdr[56] = meta["NLABL"]

    # Add the label to the output header
    tmp = zeros(UInt8, 10*80)
    for i in 1:length(meta["LABEL"])
      tmp[i] = UInt8(meta["LABEL"][i])
    end
    hdr[57:end] = reinterpret(Int32, tmp)

    # Extended Header:
    if meta["NSYMBT"] > 0
        hdr = vcat(hdr[1:end], buffer(meta["EXTHEAD"]))
    end

    handle=open(filename, "w")
    write(handle, hdr)
    return handle
end
