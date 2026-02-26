@markovjunior begin
    @rewrite area/3 b=>w # Turn one third of the pixels to white
    @rewrite ww => gg \[ x ] # Any horizontal pairs of white turn to grey
    @rewrite gg => II \[ -y ] # Any vertical pairs of grey turn to Indigo
end