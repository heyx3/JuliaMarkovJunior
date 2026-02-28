@markovjunior 'b' begin
    # Pick a source pixel.
    @rewrite 1 b => R

    # Do a backtracking random walk to carve out maze paths.
    @rewrite begin
		PRIORITIZE(earliest)
        Rbb => GGR
        GGR => Rww
    end

	@rewrite R => w
end