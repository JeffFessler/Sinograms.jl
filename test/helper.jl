# helper.jl
# utilities used in multiple test files

macro NOTinferred(ex) # flag where @inferred fails
    :($(esc(ex)))
end
