@markovjunior 'b' begin
    # Pick a source pixel.
    @rewrite 1 b => R

    # Do a backtracking random walk to carve out maze paths.
    @sequence repeat begin
        @rewrite Rbb => GGR
        @rewrite GGR => Rww
    end

	@rewrite R => w
end