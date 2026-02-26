@markovjunior 'b' begin
    # Pick a source cell.
    @rewrite 1 b => w

	# Carve out walls until none are left.
    # Mark every other cell beige instead of white,
    #   so that tunnels can't be carved right next to each other.
	@rewrite bbw => wEw

    # Clean up the beige.
	@rewrite E => w
end