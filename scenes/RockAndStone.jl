@markovjunior 'I' begin
	# Place cave seeds.
	@draw_box 'R' uv(min=0.1, size=0)
	@draw_box 'R' uv(min=0.5, size=0)
	@draw_box 'R' uv(min=0.9, size=0)

	# Grow the cave seeds and leave some veins behind.
    @rewrite (area*6) begin
        R[bI] => bR
		IRb => G__
    end

	# Clean up the veins.
	@rewrite begin
		PRIORITIZE(earliest)
		R => b
		GG => SS
		G => b
    end

	# Turn the veins into real minerals.
	@sequence repeat begin
		# Mark a vein as either Gold or Nitra.
		@rewrite  1   S => {YR}
		# Flesh out that vein.
		@rewrite begin
			YS => YY
			RS => RR  *3
		end
	end
end