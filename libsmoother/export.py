from .quarry import Quarry

try:
    import drawSvg

    HAS_DRAW_SVG = True
except:
    HAS_DRAW_SVG = False


def __write_header(out_file):
    out_file.write("##libSmoother Version: " + Quarry.get_libSmoother_version() + "\n")
    out_file.write("##libSps Version: " + Quarry.get_libSps_version() + "\n")


def __conf_session(session):
    ret = Quarry()
    ret.copy_from(session)

    ret.set_value(["settings", "interface", "add_draw_area", "val"], 0)

    return ret


def export_tsv(session, print_callback=lambda s: None):
    session = __conf_session(session)
    with open(
        session.get_value(["settings", "export", "prefix"]) + ".heatmap.tsv", "w"
    ) as out_file:
        __write_header(out_file)
        out_file.write("#chr_x\tstart_x\tend_x\tchr_y\tstart_y\tend_y\tscore\n")
        for tup in session.get_heatmap_export(print_callback):
            out_file.write("\t".join([str(x) for x in tup]) + "\n")

    if session.get_value(["settings", "interface", "show_hide", "raw"]):
        for x_axis, suff in [(True, "x"), (False, "y")]:
            with open(
                session.get_value(["settings", "export", "prefix"])
                + ".track."
                + suff
                + ".tsv",
                "w",
            ) as out_file:
                __write_header(out_file)
                out_file.write("#chr_x\tstart_x\tend_x")
                for track in session.get_track_export_names(x_axis, print_callback):
                    out_file.write("\t" + track)
                out_file.write("\n")

                for tup in session.get_track_export(x_axis, print_callback):
                    out_file.write(
                        "\t".join(
                            [
                                "\t".join([str(y) for y in x])
                                if isinstance(x, list)
                                else str(x)
                                for x in tup
                            ]
                        )
                        + "\n"
                    )


def __conv_coords(
    p, is_bottom, idx, cds, w_cds, o_cds, w_plane, o_plane, none_for_oob=False
):
    x = w_plane * (p - o_cds) / w_cds + o_plane
    if none_for_oob and (x < o_plane or x >= o_plane + w_plane):
        return None
    return max(min(x, o_plane + w_plane), o_plane)


def __cat_coords(cp, is_bottom, idx, cds, w_plane, o_plane, categories):
    c, p = cp
    return max(
        min(
            w_plane
            * (
                categories[::-1].index(c)
                + p
                + (-1 if is_bottom else 1) * cds["size"][idx] / 2
                + 0.5
            )
            / len(categories)
            + o_plane,
            o_plane + w_plane,
        ),
        o_plane,
    )


def __draw_tick_lines(
    d,
    tick_lines,
    transform,
    start,
    end,
    axis,
    stroke="black",
    stroke_width=2,
    conv_coords=__conv_coords,
    opacity=1,
    labels=None,
):
    for idx, line in enumerate(tick_lines):
        p = conv_coords(line, None, None, None, *transform, none_for_oob=True)
        if not p is None:
            if axis:
                xf = p
                xt = p
                xl = p
                yf = start if labels is None else end - 8
                yt = end
                yl = end - 15
            else:
                xf = start if labels is None else end - 8
                xt = end
                yf = p
                yt = p
                yl = p
                xl = end - 15
            if stroke_width > 0:
                d.append(
                    drawSvg.Line(
                        xf,
                        yf,
                        xt,
                        yt,
                        stroke=stroke,
                        stroke_width=stroke_width,
                        stroke_opacity=opacity,
                    )
                )
            if (
                not labels is None
                and p > transform[3]
                and p < transform[3] + transform[2]
            ):
                label = labels[idx]
                d.append(
                    drawSvg.Text(
                        label,
                        14,
                        xl,
                        yl,
                        font_family="Consolas, sans-serif",
                        transform="rotate(-90," + str(xl) + "," + str(-yl) + ")"
                        if axis
                        else "rotate(0)",
                        text_anchor="end",
                        dominant_baseline="middle",
                    )
                )


def __to_readable_pos(x, dividend, genome_end, contig_starts):
    x = int(x * dividend)
    if x >= genome_end * dividend:
        x -= contig_starts[-1] * dividend
    else:
        for start, end in zip(contig_starts, contig_starts[1:] + [genome_end]):
            if x >= start * dividend and x < end * dividend:
                x -= start * dividend
                break

    if x == 0:
        label = "0 bp"
    elif x % 1000000 == 0:
        label = "{:,}".format(x // 1000000) + " mbp"
    elif x % 1000 == 0:
        label = "{:,}".format(x // 1000) + " kbp"
    else:
        label = "{:,}".format(x) + " bp"
    return label


def __adaptive_ticker(
    d,
    transform,
    start,
    end,
    axis,
    base=10,
    mantissas=[1, 2, 5],
    desired_num_ticks=6,
    num_minor_ticks=5,
    stroke="black",
    conv_coords=__conv_coords,
    opacity=1,
    label_major=False,
    to_readable_pos=str,
):
    w_cds, o_cds, _, _ = transform
    tick_size = 0
    for n, m in [(n, m) for n in range(100) for m in mantissas]:
        tick_size = base**n * m
        if w_cds / tick_size <= desired_num_ticks:
            break

    def draw():
        fst_tick = o_cds - o_cds % tick_size
        if fst_tick < o_cds:
            fst_tick += tick_size
        ticks = []
        labels = []
        while fst_tick < w_cds + o_cds:
            if label_major:
                ticks.append(fst_tick)
                labels.append(to_readable_pos(fst_tick))
            else:
                ticks.append(fst_tick)
            fst_tick += tick_size
        __draw_tick_lines(
            d,
            ticks,
            transform,
            start,
            end,
            axis,
            stroke=stroke,
            stroke_width=1,
            conv_coords=conv_coords,
            opacity=opacity,
            labels=labels if label_major else None,
        )

    draw()
    if num_minor_ticks > 0:
        tick_size /= num_minor_ticks
        start = end - 4
        label_major = False
        draw()


def __draw_rectangles(
    d,
    cds,
    x_transform,
    y_transform,
    left="l",
    right="r",
    bottom="b",
    top="t",
    color="c",
    conv_coords_x=__conv_coords,
    conv_coords_y=__conv_coords,
):
    for idx, (l, r, b, t, c) in enumerate(
        zip(cds[left], cds[right], cds[bottom], cds[top], cds[color])
    ):
        l, r, b, t = [
            conv_coords(p, is_bottom, idx, cds, *transform)
            for p, is_bottom, transform, conv_coords in [
                (l, True, x_transform, conv_coords_x),
                (r, False, x_transform, conv_coords_x),
                (b, True, y_transform, conv_coords_y),
                (t, False, y_transform, conv_coords_y),
            ]
        ]
        d.append(drawSvg.Rectangle(l, b, r - l, t - b, fill=c))


def __draw_lines(
    d,
    cds,
    x_transform,
    y_transform,
    x="x",
    y="y",
    color="c",
    conv_coords_x=__conv_coords,
    conv_coords_y=__conv_coords,
):
    for (
        xs,
        ys,
        c,
    ) in zip(cds[x], cds[y], cds[color]):
        xs2, ysy = [
            [conv_coords(p, is_bottom, None, None, *transform) for p in ps]
            for ps, is_bottom, transform, conv_coords in [
                (xs, True, x_transform, conv_coords_x),
                (ys, True, y_transform, conv_coords_y),
            ]
        ]
        for xf, yf, xt, yt in zip(xs2[:-1], ysy[:-1], xs2[1:], ysy[1:]):
            if not float("NaN") in [xf, yf, xt, yt]:
                d.append(drawSvg.Line(xf, yf, xt, yt, stroke=c, stroke_width=2))


def __get_transform(session, w_plane, x_plane, h_plane, y_plane):
    x_transform = [
        session.get_value(["area", "x_end"]) - session.get_value(["area", "x_start"]),
        session.get_value(["area", "x_start"]),
        w_plane,
        x_plane,
    ]
    y_transform = [
        session.get_value(["area", "y_end"]) - session.get_value(["area", "y_start"]),
        session.get_value(["area", "y_start"]),
        h_plane,
        y_plane,
    ]

    return x_transform, y_transform


def __draw_heatmap(session, d, sizes, print_callback=lambda s: None):
    if sizes["show_heat"]:
        offset = 0
        if sizes["show_anno"]:
            offset += sizes["annotation"] + sizes["margin"]
        if sizes["show_secondary"]:
            offset += sizes["secondary"] + sizes["margin"]
        if sizes["show_coords"]:
            offset += sizes["coords"]
        if sizes["show_contigs"]:
            offset += sizes["contigs"]
        if sizes["show_axis"] and (sizes["show_secondary"] or sizes["show_anno"]):
            offset += sizes["axis"]

        d.append(
            drawSvg.Rectangle(
                offset,
                offset,
                sizes["heatmap"],
                sizes["heatmap"],
                fill=session.get_background_color(print_callback),
            )
        )

        x_transform, y_transform = __get_transform(
            session, sizes["heatmap"], offset, sizes["heatmap"], offset
        )

        __draw_rectangles(
            d,
            session.get_heatmap(print_callback),
            x_transform,
            y_transform,
            left="screen_left",
            right="screen_right",
            bottom="screen_bottom",
            top="screen_top",
            color="color",
        )

        if sizes["show_grid_lines"]:
            __adaptive_ticker(
                d,
                x_transform,
                offset,
                offset + sizes["heatmap"],
                True,
                num_minor_ticks=0,
                stroke="lightgrey",
                opacity=0.3,
            )
            __adaptive_ticker(
                d,
                y_transform,
                offset,
                offset + sizes["heatmap"],
                False,
                num_minor_ticks=0,
                stroke="lightgrey",
                opacity=0.3,
            )
        if sizes["show_contig_borders"]:
            __draw_tick_lines(
                d,
                session.get_tick_list(True, print_callback),
                x_transform,
                offset,
                offset + sizes["heatmap"],
                True,
                "lightgrey",
                opacity=0.5,
            )
            __draw_tick_lines(
                d,
                session.get_tick_list(False, print_callback),
                y_transform,
                offset,
                offset + sizes["heatmap"],
                False,
                "lightgrey",
                opacity=0.5,
            )
        if sizes["show_ident_line"]:
            w, x, _, _ = x_transform
            h, y, _, _ = y_transform
            bx, by, tx, ty = [
                __conv_coords(v, None, None, None, *transform)
                for v, transform in [
                    (max(x, y), x_transform),
                    (max(x, y), y_transform),
                    (min(x + w, y + h), x_transform),
                    (min(x + w, y + h), y_transform),
                ]
            ]
            print(bx, by, tx, ty)
            d.append(
                drawSvg.Line(
                    bx,
                    by,
                    tx,
                    ty,
                    stroke="lightgrey",
                    stroke_width=2,
                    stroke_opacity=0.5,
                )
            )


def __draw_annotation(session, d, sizes, print_callback=lambda s: None):
    if sizes["show_anno"]:
        offset = 0
        if sizes["show_coords"]:
            offset += sizes["coords"]
        if sizes["show_contigs"]:
            offset += sizes["contigs"]

        offset_x = offset + sizes["annotation"] + sizes["margin"]
        if sizes["show_secondary"]:
            offset_x += sizes["secondary"] + sizes["margin"]
        if sizes["show_axis"]:
            offset_x += sizes["axis"]

        x_transform, y_transform = __get_transform(
            session, sizes["heatmap"], offset_x, sizes["heatmap"], offset_x
        )
        active_anno_x = [
            sizes["annotation"],
            offset,
            session.get_displayed_annos(True, print_callback),
        ]
        active_anno_y = [
            sizes["annotation"],
            offset,
            session.get_displayed_annos(False, print_callback),
        ]

        __draw_rectangles(
            d,
            session.get_annotation(True, print_callback),
            x_transform,
            active_anno_x,
            left="screen_start",
            right="screen_end",
            bottom="anno_name",
            top="anno_name",
            color="color",
            conv_coords_y=__cat_coords,
        )
        __draw_rectangles(
            d,
            session.get_annotation(False, print_callback),
            active_anno_y,
            y_transform,
            bottom="screen_start",
            top="screen_end",
            left="anno_name",
            right="anno_name",
            color="color",
            conv_coords_x=__cat_coords,
        )

        if sizes["show_grid_lines"]:
            __adaptive_ticker(
                d,
                x_transform,
                offset,
                offset + sizes["annotation"],
                True,
                num_minor_ticks=0,
                stroke="lightgrey",
                opacity=0.3,
            )
            __adaptive_ticker(
                d,
                y_transform,
                offset,
                offset + sizes["annotation"],
                False,
                num_minor_ticks=0,
                stroke="lightgrey",
                opacity=0.3,
            )
        if sizes["show_contig_borders"]:
            __draw_tick_lines(
                d,
                session.get_tick_list(True, print_callback),
                x_transform,
                offset,
                offset + sizes["annotation"],
                True,
                "lightgrey",
                opacity=0.5,
            )
            __draw_tick_lines(
                d,
                session.get_tick_list(False, print_callback),
                y_transform,
                offset,
                offset + sizes["annotation"],
                False,
                "lightgrey",
                opacity=0.5,
            )
        if sizes["show_axis"]:
            __draw_tick_lines(
                d,
                [x + 0.5 for x in range(len(active_anno_x[2]))],
                [len(active_anno_x[2]), 0, sizes["annotation"], offset],
                offset_x - sizes["axis"],
                offset_x,
                False,
                labels=active_anno_x[2][::-1],
            )
            d.append(
                drawSvg.Line(
                    offset_x,
                    offset,
                    offset_x,
                    offset + sizes["annotation"],
                    stroke="black",
                    stroke_width=2,
                )
            )

            __draw_tick_lines(
                d,
                [x + 0.5 for x in range(len(active_anno_y[2]))],
                [len(active_anno_y[2]), 0, sizes["annotation"], offset],
                offset_x - sizes["axis"],
                offset_x,
                True,
                labels=active_anno_y[2][::-1],
            )
            d.append(
                drawSvg.Line(
                    offset,
                    offset_x,
                    offset + sizes["annotation"],
                    offset_x,
                    stroke="black",
                    stroke_width=2,
                )
            )


def __draw_secondary(session, d, sizes, print_callback=lambda s: None):
    if sizes["show_secondary"]:
        offset = 0
        if sizes["show_anno"]:
            offset += sizes["annotation"] + sizes["margin"]
        if sizes["show_coords"]:
            offset += sizes["coords"]
        if sizes["show_contigs"]:
            offset += sizes["contigs"]

        offset_x = offset + sizes["secondary"] + sizes["margin"]
        if sizes["show_axis"]:
            offset_x += sizes["axis"]

        x_transform, y_transform = __get_transform(
            session, sizes["heatmap"], offset_x, sizes["heatmap"], offset_x
        )
        min_x, max_x = session.get_min_max_tracks(True, print_callback)
        min_y, max_y = session.get_min_max_tracks(False, print_callback)
        active_anno_x = [max_x - min_x, min_x, sizes["secondary"], offset]
        active_anno_y = [max_y - min_y, min_y, sizes["secondary"], offset]

        __draw_lines(
            d,
            session.get_tracks(True, print_callback),
            x_transform,
            active_anno_x,
            x="screen_pos",
            y="values",
            color="colors",
        )
        __draw_lines(
            d,
            session.get_tracks(False, print_callback),
            active_anno_y,
            y_transform,
            y="screen_pos",
            x="values",
            color="colors",
        )

        if sizes["show_grid_lines"]:
            __adaptive_ticker(
                d,
                x_transform,
                offset,
                offset + sizes["secondary"],
                True,
                num_minor_ticks=0,
                stroke="lightgrey",
                opacity=0.3,
            )
            __adaptive_ticker(
                d,
                y_transform,
                offset,
                offset + sizes["secondary"],
                False,
                num_minor_ticks=0,
                stroke="lightgrey",
                opacity=0.3,
            )
        if sizes["show_contig_borders"]:
            __draw_tick_lines(
                d,
                session.get_tick_list(True, print_callback),
                x_transform,
                offset,
                offset + sizes["secondary"],
                True,
                "lightgrey",
                opacity=0.5,
            )
            __draw_tick_lines(
                d,
                session.get_tick_list(False, print_callback),
                y_transform,
                offset,
                offset + sizes["secondary"],
                False,
                "lightgrey",
                opacity=0.5,
            )
        if sizes["show_axis"]:
            __adaptive_ticker(
                d,
                active_anno_x,
                offset_x - sizes["axis"],
                offset_x,
                False,
                label_major=True,
            )
            d.append(
                drawSvg.Line(
                    offset_x,
                    offset,
                    offset_x,
                    offset + sizes["annotation"],
                    stroke="black",
                    stroke_width=2,
                )
            )

            __adaptive_ticker(
                d,
                active_anno_y,
                offset_x - sizes["axis"],
                offset_x,
                True,
                label_major=True,
            )
            d.append(
                drawSvg.Line(
                    offset,
                    offset_x,
                    offset + sizes["annotation"],
                    offset_x,
                    stroke="black",
                    stroke_width=2,
                )
            )


def __draw_coordinates(session, d, sizes, print_callback=lambda s: None):
    if sizes["show_coords"]:
        offset = 0
        if sizes["show_contigs"]:
            offset += sizes["contigs"]

        offset_heat = 0
        if sizes["show_anno"]:
            offset_heat += sizes["annotation"] + sizes["margin"]
        if sizes["show_secondary"]:
            offset_heat += sizes["secondary"] + sizes["margin"]
        if sizes["show_coords"]:
            offset_heat += sizes["coords"]
        if sizes["show_contigs"]:
            offset_heat += sizes["contigs"]
        if sizes["show_axis"] and (sizes["show_secondary"] or sizes["show_anno"]):
            offset_heat += sizes["axis"]

        x_transform, y_transform = __get_transform(
            session, sizes["heatmap"], offset_heat, sizes["heatmap"], offset_heat
        )
        contig_starts_x = session.get_tick_list(True, print_callback)
        contig_starts_y = session.get_tick_list(False, print_callback)

        def to_readable_pos_x(x):
            return __to_readable_pos(
                x,
                session.get_value(["dividend"]),
                contig_starts_x[-1],
                contig_starts_x[:-1],
            )

        def to_readable_pos_y(x):
            return __to_readable_pos(
                x,
                session.get_value(["dividend"]),
                contig_starts_y[-1],
                contig_starts_y[:-1],
            )

        __adaptive_ticker(
            d,
            x_transform,
            offset,
            offset + sizes["coords"],
            True,
            to_readable_pos=to_readable_pos_x,
            label_major=True,
        )
        __adaptive_ticker(
            d,
            y_transform,
            offset,
            offset + sizes["coords"],
            False,
            to_readable_pos=to_readable_pos_y,
            label_major=True,
        )
        d.append(
            drawSvg.Line(
                offset_heat,
                offset + sizes["coords"],
                offset_heat + sizes["heatmap"],
                offset + sizes["coords"],
                stroke="black",
                stroke_width=2,
            )
        )
        d.append(
            drawSvg.Line(
                offset + sizes["coords"],
                offset_heat,
                offset + sizes["coords"],
                offset_heat + sizes["heatmap"],
                stroke="black",
                stroke_width=2,
            )
        )


def minmax(v, mi, ma):
    return min(max(v, mi), ma)


def __draw_contigs(session, d, sizes, print_callback=lambda s: None):
    if sizes["show_contigs"]:
        offset = 0

        offset_heat = 0
        if sizes["show_anno"]:
            offset_heat += sizes["annotation"] + sizes["margin"]
        if sizes["show_secondary"]:
            offset_heat += sizes["secondary"] + sizes["margin"]
        if sizes["show_coords"]:
            offset_heat += sizes["coords"]
        if sizes["show_contigs"]:
            offset_heat += sizes["contigs"]
        if sizes["show_axis"] and (sizes["show_secondary"] or sizes["show_anno"]):
            offset_heat += sizes["axis"]

        x_transform, y_transform = __get_transform(
            session, sizes["heatmap"], offset_heat, sizes["heatmap"], offset_heat
        )
        w, o, _, _ = x_transform
        contig_starts_x = session.get_tick_list(True, print_callback)
        contig_centers_x = [
            (minmax(e, o, w + o) + minmax(s, o, w + o)) / 2
            for s, e in zip(contig_starts_x[:-1], contig_starts_x[1:])
        ]
        contig_starts_y = session.get_tick_list(False, print_callback)
        w, o, _, _ = y_transform
        contig_centers_y = [
            (minmax(e, o, w + o) + minmax(s, o, w + o)) / 2
            for s, e in zip(contig_starts_y[:-1], contig_starts_y[1:])
        ]
        contig_names_x = session.get_contig_ticks(True, print_callback)["contig_names"]
        contig_names_y = session.get_contig_ticks(False, print_callback)["contig_names"]
        __draw_tick_lines(
            d,
            contig_centers_x,
            x_transform,
            offset,
            offset + sizes["contigs"],
            True,
            stroke_width=0,
            labels=contig_names_x,
        )
        __draw_tick_lines(
            d,
            contig_centers_y,
            y_transform,
            offset,
            offset + sizes["contigs"],
            False,
            stroke_width=0,
            labels=contig_names_y,
        )

        y_label, x_label = session.get_value(
            ["settings", "interface", "axis_lables"]
        ).split("_")

        d.append(
            drawSvg.Text(
                x_label,
                18,
                offset_heat + sizes["heatmap"] / 2,
                offset,
                font_family="Consolas, sans-serif",
                text_anchor="middle",
                dominant_baseline="bottom",
            )
        )
        d.append(
            drawSvg.Text(
                y_label,
                18,
                offset,
                offset_heat + sizes["heatmap"] / 2,
                font_family="Consolas, sans-serif",
                transform="rotate(-90,"
                + str(offset)
                + ","
                + str(-(offset_heat + sizes["heatmap"] / 2))
                + ")",
                text_anchor="middle",
                dominant_baseline="hanging",
            )
        )


def __get_sizes(session):
    return {
        "show_heat": True,
        "show_coords": session.get_value(
            ["settings", "interface", "show_hide", "coords"]
        ),
        "show_contigs": session.get_value(
            ["settings", "interface", "show_hide", "regs"]
        ),
        "show_ident_line": session.get_value(
            ["settings", "interface", "show_hide", "indent_line"]
        ),
        "show_axis": session.get_value(["settings", "interface", "show_hide", "axis"]),
        "coords": session.get_value(["settings", "export", "coords", "val"]),
        "contigs": session.get_value(["settings", "export", "contigs", "val"]),
        "axis": session.get_value(["settings", "export", "axis", "val"]),
        "show_contig_borders": session.get_value(
            ["settings", "interface", "show_hide", "contig_borders"]
        ),
        "show_grid_lines": session.get_value(
            ["settings", "interface", "show_hide", "grid_lines"]
        ),
        "heatmap": session.get_value(["settings", "export", "size", "val"]),
        "margin": session.get_value(["settings", "export", "margins", "val"]),
        "show_anno": session.get_value(
            ["settings", "interface", "show_hide", "annotation"]
        ),
        "annotation": session.get_value(["settings", "interface", "anno_size", "val"]),
        "show_secondary": session.get_value(
            ["settings", "interface", "show_hide", "raw"]
        ),
        "secondary": session.get_value(["settings", "interface", "raw_size", "val"]),
    }


def __make_drawing(session, sizes):
    size = 0

    if sizes["show_heat"]:
        size += sizes["heatmap"] + sizes["margin"]
    if sizes["show_anno"]:
        size += sizes["annotation"] + sizes["margin"]
    if sizes["show_secondary"]:
        size += sizes["secondary"] + sizes["margin"]
    if sizes["show_coords"]:
        size += sizes["coords"]
    if sizes["show_contigs"]:
        size += sizes["contigs"]
    if sizes["show_axis"] and (sizes["show_secondary"] or sizes["show_anno"]):
        size += sizes["axis"]

    size -= sizes["margin"]

    d = drawSvg.Drawing(size, size, displayInline=False)

    d.append(
        drawSvg.Rectangle(
            0,
            0,
            size,
            size,
            fill="#ffffff",
        )
    )

    return d


def __draw(session):
    session = __conf_session(session)
    sizes = __get_sizes(session)
    d = __make_drawing(session, sizes)
    __draw_heatmap(session, d, sizes)
    __draw_annotation(session, d, sizes)
    __draw_secondary(session, d, sizes)
    __draw_coordinates(session, d, sizes)
    __draw_contigs(session, d, sizes)
    return d


def assert_has_draw_svg():
    if not HAS_DRAW_SVG:
        raise RuntimeError(
            "could not import drawSvg, is it installed? (try pip install drawSvg==1.9.0)"
        )


def export_png(session):
    assert_has_draw_svg()
    __draw(session).savePng(
        session.get_value(["settings", "export", "prefix"]) + ".png"
    )


def export_svg(session):
    assert_has_draw_svg()
    __draw(session).saveSvg(
        session.get_value(["settings", "export", "prefix"]) + ".svg"
    )
