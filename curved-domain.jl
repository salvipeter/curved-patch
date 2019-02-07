module CurvedDomain

import Graphics
using Gtk

curves = []
density = 0.1
nrlines = 10
mesh = nothing
bbox = nothing
sd = nothing
scale = 1.0
scaling = nothing

affine_combine(p, x, q) = p * (1 - x) + q * x

function generate_curves()
    global curves = []
    vert = mesh[1]
    tri = mesh[2]
    for param in 1:length(sd)
        !sd[param] && continue
        function slice(i, j)
            x, y = vert[i][param], vert[j][param]
            if y > x
                x, y = y, x
                i, j = j, i
            end
            q1 = div(x, density)
            q2 = div(y, density)
            q1 - q2 == 1 && q1 <= nrlines && affine_combine(vert[j][end-1:end],
                                                            (density * q1 - y) / (x - y),
                                                            vert[i][end-1:end])
        end
        foreach(tri) do t
            ab = slice(t[1], t[2])
            ac = slice(t[1], t[3])
            bc = slice(t[2], t[3])
            lst = filter(x -> x != false, [ab, ac, bc])
            if length(lst) == 2
                push!(curves, (lst[1], lst[2]))
            end
        end
    end
end

function draw_polygon(ctx, poly, closep = false)
    if isempty(poly)
        return
    end
    Graphics.new_path(ctx)
    Graphics.move_to(ctx, poly[1][1], poly[1][2])
    for p in poly[2:end]
        Graphics.line_to(ctx, p[1], p[2])
    end
    if closep && length(poly) > 2
        Graphics.close_path(ctx)
        Graphics.fill(ctx)
    else
        Graphics.stroke(ctx)
    end
end

function draw_mesh(ctx)
    v = mesh[1]
    for tri in mesh[2]
        draw_polygon(ctx, map(p -> scaling(p[end-1:end]), [v[tri[1]], v[tri[2]], v[tri[3]]]), true)
    end
end

function draw_segments(ctx, segments)
    for s in segments
        Graphics.new_path(ctx)
        Graphics.move_to(ctx, s[1][1], s[1][2])
        Graphics.line_to(ctx, s[2][1], s[2][2])
        Graphics.stroke(ctx)
    end
end

@guarded function draw_callback(canvas)
    ctx = Graphics.getgc(canvas)
    width = Graphics.width(canvas)
    height = Graphics.height(canvas)
    global scale = min(width, height)

    # White background
    Graphics.rectangle(ctx, 0, 0, width, height)
    Graphics.set_source_rgb(ctx, 1, 1, 1)
    Graphics.fill(ctx)

    # Input polygon
    Graphics.set_source_rgb(ctx, 0.9, 0.9, 0.9)
    draw_mesh(ctx)

    # Generated curves
    Graphics.set_source_rgb(ctx, 0, 0, 0)
    Graphics.set_line_width(ctx, 1.0)
    draw_segments(ctx, map(s -> map(scaling, s), curves))
end

function readOBJ(filename)
    vertices = []
    triangles = []
    open(filename) do f
        for line in eachline(f)
            parts = split(line)
            if parts[1] == "v"
                p = map(s -> parse(Float64, s), parts[2:end])
                push!(vertices, p)
            elseif parts[1] == "f"
                tri = map(s -> parse(Int, s), parts[2:end])
                push!(triangles, tri)
            end
        end
    end
    (vertices, triangles)
end

function run(filename)
    global mesh = readOBJ(filename)

    bboxu = extrema(map(p -> p[end-1], mesh[1]))
    bboxv = extrema(map(p -> p[end], mesh[1]))
    bblength = max(bboxu[2] - bboxu[1], bboxv[2] - bboxv[1])
    start = [bboxu[1], bboxv[1]]
    global scaling = p -> (p - start) * scale / bblength

    win = GtkWindow("Parameterization Test")
    vbox = GtkBox(:v)

    canvas = GtkCanvas(600, 600)
    draw(draw_callback, canvas)

    push!(win, vbox)
    push!(vbox, canvas)
    hbox = GtkBox(:h)
    push!(vbox, hbox)

    n = Int(length(mesh[1][1]) / 2 - 1)
    global sd = [false for _ in 1:2n]
    for i in 1:n
        cb = GtkCheckButton("s$i")
        signal_connect(cb, :toggled) do cb
            global sd[2i-1] = get_gtk_property(cb, :active, Bool)
            generate_curves()
            draw(canvas)
        end
        push!(hbox, cb)
        cb = GtkCheckButton("d$i")
        signal_connect(cb, :toggled) do cb
            global sd[2i] = get_gtk_property(cb, :active, Bool)
            generate_curves()
            draw(canvas)
        end
        push!(hbox, cb)
    end

    hbox = GtkBox(:h)
    push!(vbox, hbox)

    push!(hbox, GtkLabel("Density: "))
    sb = GtkSpinButton(0.01, 0.2, 0.05)
    set_gtk_property!(sb, :value, density)
    signal_connect(sb, :value_changed) do sb
        global density = get_gtk_property(sb, :value, Float64)
        generate_curves()
        draw(canvas)
    end
    push!(hbox, sb)

    push!(hbox, GtkLabel("# of lines: "))
    sb = GtkSpinButton(1, 100, 1)
    set_gtk_property!(sb, :value, nrlines)
    signal_connect(sb, :value_changed) do sb
        global nrlines = get_gtk_property(sb, :value, Int)
        generate_curves()
        draw(canvas)
    end
    push!(hbox, sb)

    showall(win)
    generate_curves()
    draw(canvas)
end

if length(ARGS) == 1
    run(ARGS[1])
end

end # module
