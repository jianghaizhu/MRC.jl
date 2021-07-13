
# Returns a template for MRC metadata
function MRCmeta()
    meta = Dict{String, Any}(
        "Filename" => "",

        "NX" => Int32(1),   # The number of columns
        "NY" => Int32(1),   # The number of rows
        "NZ" => Int32(1),   # The number of sections

        #= Mode
            0 == 8-bit  signed integer        (Int8)
            1 == 16-bit signed integer        (Int16)
            2 == 32-bit signed real           (Float32)
            3 == 16-bit unsigned integer      (UFloat16)
            4 == 32-bit unsigned reals        (UFloat32)
            6 == 16-bit unsigned integer      (UInt16)
            7 == 32-bit signed integer        (Int32)
           16 == RGB data in 3 unsigned bytes (RGB)         =#
        "MODE" => Int32(0),

        "NXSTART" => Int32(0),  # The number of the first column in the map (Default: 0)
        "NYSTART" => Int32(0),  # The number of the first row in the map
        "NZSTART" => Int32(0),  # The number of the first section in the map
        "MX" => Int32(1),       # The number of intervals along the X-axis
        "MY" => Int32(1),       # The number of intervals along the Y-axis
        "MZ" => Int32(1),       # The number of intervals along the Z-axis

        #=  Cell dimensions in ångströms (A,B,C):   =#
        "CELL_A" => Float32(0),
        "CELL_B" => Float32(0),
        "CELL_C" => Float32(0),

        #=  Cell angles in degrees: =#
        "CELL_ALPHA" => Float32(90),
        "CELL_BETA"  => Float32(90),
        "CELL_GAMMA" => Float32(90),

        #= Axis Alignment: =#
        "MAPC" => Float32(1),   # Axis corresponding to colmns:  (1, 2, 3) for: (X, Y. Z)
        "MAPR" => Float32(2),   # Axis corresponding to row:     (1, 2, 3) for: (X, Y. Z)
        "MAPS" => Float32(3),   # Axis corresponding to section: (1, 2, 3) for: (X, Y. Z)

        #= Density Values: =#
        "DMIN"  => Float32(0),   # Minimum density value
        "DMAX"  => Float32(0),   # Maximum density value
        "DMEAN" => Float32(0),   # Mean density value

        #= Space group number
            0:      implies a 2D image/ image stack
            1:      implies a single 3D volume
            401+:   implies a 3D volume stack      =#
        "ISPG" => Int32(0),

        # Length of the extended header  whether it contains symmetry records or additional metadata
        "NSYMBT"  => Int32(0),

        # Extra space used for anything - 0 by default
        "EXTRA_SPACE1" => zeros(Int32, 2),

        # Identifies the extended metadata type
        "EXTTYPE" => String("    "),

        #= Version of the MRC format.
            The version of the MRC format that the file adheres to
            Specified as:
                Year * 10 + version within the year (base 0).       =#
        "NVERSION" => Int32(0),

        # Extra space used for anything - 0 by default
        "EXTRA_SPACE2" => zeros(Int32, 21),

        #= ORIGIN (X Y Z)
            For transforms (Mode 3 or 4):
                ORIGIN is the phase origin of the transformed image in pixels, e.g. as used in
                helical processing of the MRC package. For a transform of a padded image, this
                value corresponds to the pixel position in the padded image of the center of
                the unpadded image.

            For other modes:
                ORIGIN specifies the real space location of a subvolume taken from a larger
                volume. In the (2-dimensional) example shown above, the header of the map
                containing the subvolume (red rectangle) would contain ORIGIN = 100, 120 to
                specify its position with respect to the original volume (assuming the original
                volume has its own ORIGIN set to 0, 0).                                             =#
        "ORIGIN_X" => Float32(0),
        "ORIGIN_Y" => Float32(0),
        "ORIGIN_Z" => Float32(0),

        #= Character string 'MAP '
            Used to determine the endianness of the the image data. =#
        "MAP"    => "    ",
        "MACHST" => "    ",    # Machine stamp

        "RMS"    => Float32(0),    # RMS deviation of map from mean density

        "NLABL"  => Int32(0),    # The number of labels being used
        "LABEL"  => repeat(" ", 10*80),  # 10 80-character text labels

        #= Dimensional description of the MRC file.

            0       Indicates that the MRC variable is a single 2D image
            1       Indicates that the MRC variable is a stack of 2D images
            2       Indicates that the MRC variable is a single 3D volume
            3       Indicates that the MRC variable is a stack of 3D volumes
        =#
        "DIM" => Int8(0),

        "EXTHEAD" => Array{UInt8}([])    # extended header
    )
    machineLE = (ENDIAN_BOM == 0x04030201) # True for Little-Endian machine
    if machineLE
        meta["MAP"]    = reinterpret(Int32, b"MAP ")[1]
        meta["MACHST"] = reinterpret(Int32, UInt8[68, 65, 0, 0])[1]
    else
        meta["MAP"]    = reinterpret(Int32, b" PAM")[1]
        meta["MACHST"] = reinterpret(Int32, UInt8[0, 0, 65, 68])[1]
    end
    return meta
end

#=
    Description:
        Return the MRC dimensional description.

    Returns:
        0       Meaning that the MRC variable is a single 2D image
        1       Meaning that the MRC variable is a stack of 2D images
        2       Meaning that the MRC variable is a single 3D volume
        3       Meaning that the MRC variable is a stack of 3D volumes
=#
function DIM(Data::Dict)
    if Data["ISPG"] == 0
        if Data["MZ"] == 1
            return 0
        elseif  Data["MZ"] > 1
            return 1
        else
            println("ISPG = ", Data["ISPG"])
            println("MZ = ", Data["MZ"])
            println("Unknown dimensional specification.")
        end
    elseif Data["ISPG"] == 1
        return 2
    elseif Data["ISPG"] > 400
        return 3
    else
        println("ISPG = ", Data["ISPG"])
        println("MZ = ", Data["MZ"])
        println("Unknown dimensional specification.")
    end
end

#=
    Description:
        determine the image format origin based on EXTTYPE

    EXTTYPE:
        CCP4:       CCP4
        MRCO:       Original MRCO
        SERI:       SerialEM
        AGAR:       Agard
        FEI1:       FEI software, e.g. EPU and Xplore3D, Amira, Avizo
=#
function formatorigin(Data::Dict)
    if Data["EXTTYPE"] == "CCP4"
        return "CCP4 Format"
    elseif Data["EXTTYPE"] == "MRCO"
        return "Original MRCO Format"
    elseif Data["EXTTYPE"] == "SERI"
        return "SerialEM Format"
    elseif Data["EXTTYPE"] == "AGAR"
        return "Agard Format"
    elseif Data["EXTTYPE"] == "FEI1"
        return "FEI Format"
    else
        println("Unknown Format Origin!")
    end
end
