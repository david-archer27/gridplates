function read_globsed_xyz()
    f = open("../GlobSed-v3.xyz")
    field = fill(0.,nx,ny)
    lines = readlines(f)
    for iline in 1:length(lines)
        line = lines[iline]
        words = split(line)
        longt = parse(Float64,words[1])
        lat = parse(Float64,words[2])
        #println(longt, " ", lat)
        if longt == floor(longt) + 0.5 && lat == floor(lat) + 0.5
            ix = Int(longt - 0.5 + 181)
            iy = Int(lat - 0.5 + 91)
            thickness = parse(Float64,words[3])
            if thickness != thickness
                thickness = 0.
            end
            field[ix,iy] = thickness
            #println(longt, " ", lat," ",ix," ",iy," ",thickness)
        end
    end
    return field
end