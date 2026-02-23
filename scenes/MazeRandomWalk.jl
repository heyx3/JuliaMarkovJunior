@markovjunior 'b' begin
    # Pick a source pixel.
    @rewrite (area/500) b => R

    # Do a backtracking random walk to carve out maze paths.
    #TODO: Nest in a repeating sequence
    @rewrite Rbb => GGR
    @rewrite GGR => Rww
    @rewrite R => w
end