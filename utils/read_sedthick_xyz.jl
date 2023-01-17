function read_sedthick_xyz()
    f = open( "laske.sedmap.xyz")
    lines = readlines(f)
    field = fill(0.,nx,ny)
    for line in lines
        words = split(line)
        longt = parse(Float64,words[1]); lat = parse(Float64,words[2])
        if abs( longt - floor( longt ) - 0.5 ) < 0.01 && abs( lat - floor( lat ) - 0.5 ) < 0.01
            ix = Int(longt + 0.5 ) + 180; iy = 90 + Int(lat + 0.5)
            if iy == 0 
                error("[longt,lat][ix,iy]",longt," ",lat," ",ix," ",iy)
            end
            if words[3] != "NaN"
                field[ix,iy] = parse(Float64,words[3])
            end
        end
    end
    return field
end
