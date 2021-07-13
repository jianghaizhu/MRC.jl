# This code reads in an MRC file
module MRC
using ImageMetadata
const MRCImg = ImageMeta
export MRCImg
export readMRC, writeMRC
export MRCmeta, DIM, formatorigin

include("mrc_h.jl")
include("readMRC.jl")
include("writeMRC.jl")
end
