@markovjunior 'I' begin
	# Place a seed.
	@fill 'b' uv(center=0.5, size=0)

	# Randomly carve outward, and
	#    randomly carve outward with partial (Slate) rock.
    @rewrite area*5 begin
        bI => Sb
		SI => bS
		Sb => bb
		Sb => SS /3
	end

	# Add a dripping grey shadow under the Slate rock.
	@rewrite begin
		Sb => gb %1.0 \[+y]
		# Drip, and optionally forbid further drip using a Green pixel.
		gbbb => ggg{Gb} %0.4 \[+y]
	end
	@rewrite G => b

	# Finalize colors.
	@rewrite 0 g => S %0.7

	# Add walls around the edges.
	@rewrite [bg]I => _S
end