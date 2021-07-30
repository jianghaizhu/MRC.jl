function readMRC(filename::AbstractString, startSlice::Int=1, numSlices::Int=0)
    occursin(r".+\.mrcs?$", filename) || error("The input file does not have an expected MRC extension.")

    file = open(filename)

    # MRC header
    head = MRCmeta()
    head["Filename"] = filename
    head = readMRCheader(file, head)
    head = readMRCextendedheader(file, head)

    # MRC image data
    img = readMRCimage(file, head, startSlice, numSlices)

    close(file)

    return MRCImg(img, head)
end

#=
    Create an array of type arrayType and length arrayLength.
=#
createarray(arrayType::Type, arrayLength) = return Array{arrayType}(undef, arrayLength)

#=
    Read the first 10 entries of the MRC header, which are integers.
=#
function readlayout(file::IOStream, head::Dict{String, Any})
    #  1 (1-4 bytes):   NX (number of columns, fastest changing in map)
    #  2 (4-8):         NY (number of rows)
    #  3 (9-12):        NZ (number of sections, slowest changing in map)
    #  4 (13-16):       MODE (0, 8-bit signed integer, -128 to 127; 1, 16-bit signed integer; 2, 32-bit signed real;
    #                        3, transform: complex 16-bit integer; 4, transform: complex 32-bit reals;
    #                        6, 16-bit unsigned integer; 7, 32-bit signed integer; 16, RGB data in 3 unsigned bytes)
    #  5 (17-20):       NXSTART (number of first column in map, Default = 0)
    #  6 (21-24):       NYSTART (number of first row in map)
    #  7 (24-28):       NZSTART (number of first section in map)
    #  8 (29-32):       MX (number of intervals along X)
    #  9 (33-36):       MY (number of intervals along Y)
    # 10 (37-40):       MZ (number of intervals along Z; In EM, where there is no unit cell. MZ represents the number
    #                        of sections in a single volume. For a volume stack, NZ/MZ will the number of volumes in the stack)
    layout = createarray(Int32, 10)
    read!(file, layout)

    head["NX"]       = layout[1]
    head["NY"]       = layout[2]
    head["NZ"]       = layout[3]
    head["MODE"]     = layout[4]
    head["NXSTART"]  = layout[5]
    head["NYSTART"]  = layout[6]
    head["NZSTART"]  = layout[7]
    head["MX"]       = layout[8]
    head["MY"]       = layout[9]
    head["MZ"]       = layout[10]

    return head
end

#=
    Read the next 12 entries of the MRC header, which are floats.
=#
function readaxes(file::IOStream, head::Dict{String, Any})
    # 11,12,13 (41-52): the first three are the cell dimensions in angstroms: CELLA (A B C)
    # 14,15,16 (53-64): the next three are cell angles in degrees: CELLB (α β γ)
    # 17 (65-68):       MAPC (axis corresponding to columns: 1,2,3 for X,Y,Z)
    # 18 (69-72):       MAPR (axis corresponding to rows: 1,2,3 for X,Y,Z)
    # 19 (73-76):       MAPS (axis corresponding to sections: 1,2,3 for X,Y,Z)
    # 20 (77-80):       DMIN (minimum density value)
    # 21 (81-84):       DMAX (maximum density value)
    # 22 (85-88):       DMEAN (mean density value)
    axes = createarray(Float32, 12)
    read!(file, axes)

    head["CELL_A"]     = axes[1]
    head["CELL_B"]     = axes[2]
    head["CELL_C"]     = axes[3]
    head["CELL_ALPHA"] = axes[4]
    head["CELL_BETA"]  = axes[5]
    head["CELL_GAMMA"] = axes[6]
    head["MAPC"]       = axes[7]
    head["MAPR"]       = axes[8]
    head["MAPS"]       = axes[9]
    head["DMIN"]       = axes[10]
    head["DMAX"]       = axes[11]
    head["DMEAN"]      = axes[12]

    return head
end

#=
    Read the next 27 entries of the MRC header, which are integers.
=#
function readspace(file::IOStream, head::Dict{String, Any})
    # 23 (89-92):       ISPG, space group number. Spacegroup 0 implies a 2D image or image stack.
    #                   For crystallography, ISPG represents the actual spacegroup. For single volumes
    #                   from EM/ET, the spacegroup should be 1. For volume stacks, we adopt the convention
    #                   that ISPG is the spacegroup number + 400, which in EM/ET will typically be 401.
    # 24 (93-96):       NSYMBT, NSYMBT specifies the size of the extended header in bytes, whether it
    #                   contains symmetry records (as in the original format definition) or any other
    #                   kind of additional metadata.
    # 25-49 (97-196):   extra space used for anything - 0 by default.
    # 27 (105-108):     EXTTYPE, (A 4-byte string); code for the type of extended header.
    # 28 (109-112):     NVERSION, version of the MRC format. The version of the MRC format that the file
    #                   adheres to, specified as: Year * 10 + version within the year (base 0). At the time
    #                   of this writing, the current format would have the value: 20140.
    space = createarray(Int32, 4)
    read!(file, space)

    head["ISPG"]         = space[1]
    head["NSYMBT"]       = space[2]
    head["EXTRA_SPACE1"] = space[3:4]
    
    # Read in the C-style chars as UInt8 and convert them to an String
    head["EXTTYPE"]      = String(read(file, 4))
    
    space = createarray(Int32, 22)
    read!(file, space)

    head["NVERSION"]     = space[1]
    head["EXTRA_SPACE2"] = space[2:22]

    return head
end

#=
    Read the next 3 entries of the MRC header, which are floats.
=#
function readorigin(file::IOStream, head::Dict{String, Any})
    # 50-52 (197-208):  ORIGIN (X Y Z), phase origin (pixels) or origin of subvolume (A).
    #                   For transforms (Mode 3 or 4):
    #                       ORIGIN is the phase origin of the transformed image in pixels, e.g. as used in
    #                       helical processing of the MRC package. For a transform of a padded image, this
    #                       value corresponds to the pixel position in the padded image of the center of
    #                       the unpadded image.
    #
    #                   For other modes:
    #                       ORIGIN specifies the real space location of a subvolume taken from a larger
    #                       volume. In the (2-dimensional) example shown above, the header of the map
    #                       containing the subvolume (red rectangle) would contain ORIGIN = 100, 120 to
    #                       specify its position with respect to the original volume (assuming the original
    #                       volume has its own ORIGIN set to 0, 0).
    origin = createarray(Float32, 3)
    read!(file, origin)

    head["ORIGIN_X"] = origin[1]
    head["ORIGIN_Y"] = origin[2]
    head["ORIGIN_Z"] = origin[3]

    return head
end

#=
    Read main header of MRC file.
=#
function readMRCheader(file::IOStream, head::Dict{String, Any})
    # The main header of MRC file is 1024 bytes in length.

    head = readlayout(file, head)   # Read the first 10 entries, which are integers.
    head = readaxes(file, head)     # Read the next 12 entries, which are floats.
    head = readspace(file, head)    # Read the next 27 entries, which are integers.
    head = readorigin(file, head)   # Read the next 3 entries, which are floats.

    # Read entries 53 and 54, which are text.
    head["MAP"]      = String(read(file, 4))    # 53 (209-212): MAP, Character string ‘MAP ’ to identify file type.
    head["MACHST"]   = String(read(file, 4))    # 54 (213-216): MACHST, Machine stamp.

    # Read entry 55, which is a float.
    head["RMS"]   = Float32(read(file, 4)[1])   # 55 (217-220): RMS, rms deviation of map from mean density.

    # Read entry 56, which is an integer.
    head["NLABL"] = Int32(read(file, 4)[1])     # 56 (221-224): NLABL, number of labels being used.

    # Read entries 57 to 256, which are text.
    head["LABEL"] = String(read(file, 10*80))   # 57-256 (225-1024): LABEL(20,10), 10 80-character text labels.

    head["DIM"]   = DIM(head)

    return head
end

#=
    Read extended header of MRC file.
=#
function readMRCextendedheader(file::IOStream, head::Dict{String, Any})
    # Length of extended header is specified by the `NSYMBT` in terms of number of bytes.
    if head["NSYMBT"] > 0
        head["EXTHEAD"] = read(file, head["NSYMBT"])
    end
    return head
end

#=
    Determine pixbytes as specified by `head["MODE"]`.
=#
function modepixbytes(head::Dict{String, Any})
    # 0, 8-bit signed integer, -128 to 127.
    if head["MODE"] == 0
        pixbytes = 1
    # 1, 16-bit signed integer.
    elseif head["MODE"] == 1
        pixbytes = 2
    # 2, 32-bit signed real.
    elseif head["MODE"] == 2
        pixbytes = 4
    # 3, transform: complex 16-bit integer.
    elseif head["MODE"] == 3
        pixbytes = 4
    # 4, transform: complex 32-bit reals.
    elseif head["MODE"] == 4
        pixbytes = 8
    # 6, 16-bit unsigned integer.
    elseif head["MODE"] == 6
        pixbytes = 2
    # 7, 32-bit signed integer.
    elseif head["MODE"] == 7
        pixbytes = 4
    # 16, RGB data in 3 unsigned bytes
    elseif head["MODE"] == 16
        error("readMRC: Sorry, RGB image processing is not currently supported.")
    else
        error("readMRC: Unknown Data Mode: $(head["MODE"]).")
    end
    return pixbytes
end

#=
    If a user wants to read an MRC image starting at a certain slice (startSlice),
    then the data before the StartSlice needs to be skipped. The size of each
    slice is NX * NY * pixbytes.
    If a user only wants to read a certain number of slices, then the number of
    Slices (numSlices) needs to be defined.
    If the number of slices a user wants to read is larger than the existing slices,
    then read the existing slices.
=#
function slices(file::IOStream, head::Dict{String, Any}, startSlice::Int, numSlices::Int)
    pixbytes = modepixbytes(head)
    zSlices = head["NZ"]
    skipBytes = 0

    if startSlice > 1
        skipBytes = (startSlice - 1) * head["NX"] * head["NY"] * pixbytes
        zSlices = head["NZ"] - (startSlice - 1)
        skip(file, skipBytes)
    end

    if numSlices > 0
        # Choose the smaller number between the existing slices and what a user wants.
        zSlices = min((head["NZ"] - (startSlice - 1)), numSlices)
    end

    zSlices = Int(zSlices)
    skipBytes = Int(skipBytes)

    return zSlices, skipBytes
end

#=
    Determine the length of the array that holds the image data.
=#
function lengtharray(file::IOStream, head::Dict{String, Any}, skipBytes::Int)
    lengthArray = Int((filesize(file) - head["NSYMBT"] - skipBytes - 1024)/4)
    return lengthArray
end

#=
    Create an array to hold the MRC image data as specified by `head["MODE"]`.
=#
function imagearray(head::Dict{String, Any}, lengthArray::Int)
    if head["MODE"] == 0
        img = createarray(UInt8, lengthArray)
    elseif head["MODE"] == 1
        img = createarray(Int16, lengthArray)
    elseif head["MODE"] == 2
        img = createarray(Float32, lengthArray)
    elseif head["MODE"] == 3
        img = createarray(Int16, lengthArray)
    elseif head["MODE"] == 4
        img = createarray(Pair{Float32}, lengthArray)
    elseif head["MODE"] == 6
        img = createarray(UInt16, lengthArray)
    end
    return img
end

#=
    Read image data in the format specified by `head["MODE"]`.
=#
function readMRCimage(file::IOStream, head::Dict{String, Any}, startSlice::Int, numSlices::Int)    
    zSlices, skipBytes = slices(file, head, startSlice, numSlices)
    lengthArray = lengtharray(file, head, skipBytes)

    # Error occurs if the number of bytes skipped is greater than the number of bytes for the image data.
    skipBytes ≤ lengthArray || error("readMRC: Sorry, cannot start at slice $startSlice.")
    # Error occurs if the dimensions do not equal the number of bytes for the image data.
    lengthArray == head["NX"] * head["NY"] * zSlices || error("readMRC: Sorry, cannot read $numSlices slices.")
    
    img = imagearray(head, lengthArray)
    read!(file, img)

    if head["NZ"] == 1
        # MRC is row major, permute dims 1 and 2 to make it column major
        img = reshape(img, head["NY"], head["NX"])
        img = permutedims(img, [2, 1])
    else
        # MRC is row major, permute dims 1 and 2 to make it column major
        img = reshape(img, head["NY"], head["NX"], zSlices)
        img = permutedims(img, [2, 1, 3])
    end
    return img
end