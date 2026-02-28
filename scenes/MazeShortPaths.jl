@markovjunior 'b' begin
    # Pick one source cell.
    @rewrite 1 b => w

	# Carve out walls from the source until none are left.
    # Mark alternating spaces Green instead of white,
    #   so that tunnels can't be carved right next to each other.
	@rewrite bbw => wGw

    # Clean up the Green.
	@rewrite G => w
end