
function readMRC(filename::AbstractString, startSlice::Int=1, numSlices::Int=0)
    occursin(r".+\.mrcs?$", filename) || error("The input file does not have an expected MRC extension.")

    isStack = occursin(r".+\.mrcs$", filename)

    file = open(filename)
    head = MRCmeta()

    head["Filename"] = filename

    #The main header of MRC file is 1024 bytes in length.

    # Read the first 10 entries, which are integers.
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
    layout  = read(file, Int32, 10)

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

    # Read the next 12 entries, which are floats.
    # 11,12,13 (41-52): the first three are the cell dimensions in angstroms: CELLA (A B C)
    # 14,15,16 (53-64): the next three are cell angles in degrees: CELLB (α β γ)
    # 17 (65-68):       MAPC (axis corresponding to columns: 1,2,3 for X,Y,Z)
    # 18 (69-72):       MAPR (axis corresponding to rows: 1,2,3 for X,Y,Z)
    # 19 (73-76):       MAPS (axis corresponding to sections: 1,2,3 for X,Y,Z)
    # 20 (77-80):       DMIN (minimum density value)
    # 21 (81-84):       DMAX (maximum density value)
    # 22 (85-88):       DMEAN (mean density value)
    axes = read(file, Float32, 12)

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

    # Read the next 27 entires, which are integers.
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
    space = read(file, Int32, 4)
    head["ISPG"]         = space[1]
    head["NSYMBT"]       = space[2]
    head["EXTRA_SPACE1"] = space[3:4]
    head["EXTTYPE"]      = String(read(file, UInt8, 4)) # Read in the C-style chars as UInt8 and convert them to an String
    space = read(file, Int32, 22)
    head["NVERSION"]     = space[1]
    head["EXTRA_SPACE2"] = space[2:22]

    # Read the next 3 entires, which are floats.
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
    origin = read(file, Float32, 3)
    head["ORIGIN_X"] = origin[1]
    head["ORIGIN_Y"] = origin[2]
    head["ORIGIN_Z"] = origin[3]

    # Read entries 53 and 54, which are text.
    # 53 (209-212):   MAP, Character string ‘MAP ’ to identify file type.
    # 54 (213-216):   MACHST, Machine stamp.
    head["MAP"]    = String(read(file, UInt8, 4))
    head["MACHST"] = String(read(file, UInt8, 4))

    # Read entry 55, which is a float.
    # 55 (217-220):   RMS, rms deviation of map from mean density.
    head["RMS"]   = read(file, Float32, 1)[1]

    # Read entry 56, which is an integer.
    # 56 (221-224):   NLABL, number of labels being used.
    head["NLABL"] = read(file, Int32, 1)[1]

    # Read entries 57 to 256, which are text.
    # 57-256 (225-1024):    LABEL(20,10), 10 80-character text labels.
    head["LABEL"] = String(read(file, UInt8, 10*80))

    head["DIM"]   = DIM(head)

    # Read any extended header data.
    # Length of extended header is specified by the `NSYMBT` in terms of number of bytes.
    if head["NSYMBT"] > 0
        head["EXTHEAD"] = read(file, UInt8, head["NSYMBT"])
    end

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
        error("readMRC: Unknown Data Mode: ", head["MODE"])
    end

    SkipBytes = 0
    #=
        If a user wants to read an MRC image starting at a certain slice (startSlice),
        then the data before the StartSlice needs to be skipped. The size of each
        slice is NX * NY * pixbytes.
        If a user only wants to read a certain number of slices, then the number of
        Slices (numSlices) needs to be defined.
        If the number of slices a user wants to read is larger than the existing slices,
        then read the existing slices.
    =#
    if startSlice > 1
        skipBytes = (startSlice - 1) * head["NX"] * head["NY"] * pixbytes
        skip(file, skipBytes)
    end

    if numSlices == 0
        zSlices = head["NZ"]
    else
        #choose the smaller number between the existing slices and what a user wants.
        zSlices = min((head["NZ"] - (startSlice - 1)), numSlices)
    end

    # ndata = Int64(head["NX"] * head["NY"] * head["NZ"])

    #=
        Read image data in the format specified by `head["MODE"]`
    =#
    if head["NZ"] == 1
        if head["MODE"] == 0
            img = read(file, UInt8, head["NX"], head["NY"])
        elseif head["MODE"] == 1
            img = read(file, Int16, head["NX"], head["NY"])
        elseif head["MODE"] == 2
            img = read(file, Float32, head["NX"], head["NY"])
        elseif head["MODE"] == 3
            img = read(file, Int16, head["NX"], head["NY"])
        elseif head["MODE"] == 4
            img = read(file, Pair{Float32}, head["NX"], head["NY"])
        elseif head["MODE"] == 6
            img = read(file, UInt16, head["NX"], head["NY"])
        end
        # MRC is row major, permute dims 1 and 2 to make it column major
        img = permutedims(img, [2, 1])
    else
        if head["MODE"] == 0
            img = read(file, UInt8, head["NX"], head["NY"], zSlices)
        elseif head["MODE"] == 1
            img = read(file, Int16, head["NX"], head["NY"], zSlices)
        elseif head["MODE"] == 2
            img = read(file, Float32, head["NX"], head["NY"], zSlices)
        elseif head["MODE"] == 3
            img = read(file, Int16, head["NX"], head["NY"], zSlices)
        elseif head["MODE"] == 4
            img = read(file, Pair{Float32}, head["NX"], head["NY"], zSlices)
        elseif head["MODE"] == 6
            img = read(file, UInt16, head["NX"], head["NY"], zSlices)
        end
        # MRC is row major, permute dims 1 and 2 to make it column major
        img = permutedims(img, [2, 1, 3])
    end

    close(file)

#    if isStack && (head["DIM"] != 1)
#        println("WARNING: Metadata indicates that this is NOT an image stack.")
#        println("However, the file is using the extension '.mrcs'. Please be sure to remedy the situation.")
#    elseif !isStack && (head["DIM"] == 1)
#        println("WARNING: Metadata indicates that this is an image stack.")
#        print("However, the file is NOT using the extension '.mrcs'.")
#        println("Please be sure to remedy the situation.")
#    end
    return MRCImg(img, head)
end
