module Testing

export src_dir
export test_dir

src_dir = @__DIR__ # current directory
test_dir = joinpath(src_dir, "..", "test") 


end
