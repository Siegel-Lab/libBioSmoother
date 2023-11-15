from .quarry import Quarry
import os

try:
    import drawSvg

    HAS_DRAW_SVG = True
except:
    HAS_DRAW_SVG = False


def __get_header(session):
    return (
        "##libBioSmoother Version: "
        + Quarry.get_libBioSmoother_version()
        + "\n"
        + "##libSps Version: "
        + Quarry.get_libSps_version()
        + "\n"
        + "## "
        + session.get_readable_area()
        + "\n"
    )


def __conf_session(session):
    ret = Quarry()
    ret.copy_from(session)

    ret.set_value(["settings", "interface", "add_draw_area", "val"], 0)

    return ret


def __conv_coords(
    p, is_bottom, idx, cds, w_cds, o_cds, w_plane, o_plane, none_for_oob=False
):
    x = w_plane * (p - o_cds) / w_cds + o_plane
    if none_for_oob and (x < o_plane or x > o_plane + w_plane):
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
        label = "{:,}".format(x // 1000000) + " Mbp"
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
    stroke_width=2,
    stroke_linecap="round",
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
                d.append(
                    drawSvg.Line(
                        xf,
                        yf,
                        xt,
                        yt,
                        stroke=c,
                        stroke_width=stroke_width,
                        stroke_linecap=stroke_linecap,
                    )
                )


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


def __draw_region(session, d, sizes, print_callback=lambda s: None):
    if sizes["show_region"]:
        offset_y = 0
        if sizes["show_anno_y"]:
            offset_y += sizes["annotation"] + sizes["margin"]
        if sizes["show_secondary_y"]:
            offset_y += sizes["secondary"] + sizes["margin"]
        if sizes["show_coords_y"]:
            offset_y += sizes["coords"]
        if sizes["show_contigs_y"]:
            offset_y += sizes["contigs"]
        if sizes["show_axis"] and (sizes["show_secondary_x"] or sizes["show_anno_x"]):
            offset_y += sizes["axis"]
        if sizes["show_heat"]:
            offset_y += sizes["heatmap"]
            offset_y += sizes["margin"]

        d.append(
            drawSvg.Text(
                session.get_readable_area(),
                14,
                d.width / 2,
                offset_y,
                font_family="Consolas, sans-serif",
                text_anchor="middle",
                dominant_baseline="bottom",
            )
        )


def __draw_heatmap(session, d, sizes, print_callback=lambda s: None):
    if sizes["show_heat"]:
        offset_x = 0
        if sizes["show_anno_x"]:
            offset_x += sizes["annotation"] + sizes["margin"]
        if sizes["show_secondary_x"]:
            offset_x += sizes["secondary"] + sizes["margin"]
        if sizes["show_coords_x"]:
            offset_x += sizes["coords"]
        if sizes["show_contigs_x"]:
            offset_x += sizes["contigs"]
        if sizes["show_axis"] and (sizes["show_secondary_y"] or sizes["show_anno_y"]):
            offset_x += sizes["axis"]
        offset_y = 0
        if sizes["show_anno_y"]:
            offset_y += sizes["annotation"] + sizes["margin"]
        if sizes["show_secondary_y"]:
            offset_y += sizes["secondary"] + sizes["margin"]
        if sizes["show_coords_y"]:
            offset_y += sizes["coords"]
        if sizes["show_contigs_y"]:
            offset_y += sizes["contigs"]
        if sizes["show_axis"] and (sizes["show_secondary_x"] or sizes["show_anno_x"]):
            offset_y += sizes["axis"]

        d.append(
            drawSvg.Rectangle(
                offset_x,
                offset_y,
                sizes["heatmap"],
                sizes["heatmap"],
                fill=session.get_background_color(print_callback),
            )
        )

        x_transform, y_transform = __get_transform(
            session, sizes["heatmap"], offset_x, sizes["heatmap"], offset_y
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
                offset_y,
                offset_y + sizes["heatmap"],
                True,
                num_minor_ticks=0,
                stroke="lightgrey",
                opacity=0.3,
            )
            __adaptive_ticker(
                d,
                y_transform,
                offset_x,
                offset_x + sizes["heatmap"],
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
                offset_y,
                offset_y + sizes["heatmap"],
                True,
                "lightgrey",
                opacity=0.5,
            )
            __draw_tick_lines(
                d,
                session.get_tick_list(False, print_callback),
                y_transform,
                offset_x,
                offset_x + sizes["heatmap"],
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
    if sizes["show_anno_y"]:
        offset = 0
        if sizes["show_coords_y"]:
            offset += sizes["coords"]
        if sizes["show_contigs_y"]:
            offset += sizes["contigs"]

        offset_x = 0
        if sizes["show_coords_x"]:
            offset_x += sizes["coords"]
        if sizes["show_contigs_x"]:
            offset_x += sizes["contigs"]
        if sizes["show_anno_x"]:
            offset_x += sizes["annotation"] + sizes["margin"]
        if sizes["show_secondary_x"]:
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
        if sizes["show_axis"]:
            __draw_tick_lines(
                d,
                [x + 0.5 for x in range(len(active_anno_x[2]))],
                [len(active_anno_x[2]), 0, sizes["annotation"], offset],
                offset_x - sizes["axis"],
                offset_x - sizes["spacing"],
                False,
                labels=active_anno_x[2][::-1],
            )
            d.append(
                drawSvg.Line(
                    offset_x - sizes["spacing"],
                    offset,
                    offset_x - sizes["spacing"],
                    offset + sizes["annotation"],
                    stroke="black",
                    stroke_width=2,
                )
            )

    if sizes["show_anno_x"]:
        offset = 0
        if sizes["show_coords_x"]:
            offset += sizes["coords"]
        if sizes["show_contigs_x"]:
            offset += sizes["contigs"]

        offset_x = 0
        if sizes["show_coords_y"]:
            offset_x += sizes["coords"]
        if sizes["show_contigs_y"]:
            offset_x += sizes["contigs"]
        if sizes["show_anno_y"]:
            offset_x += sizes["annotation"] + sizes["margin"]
        if sizes["show_secondary_y"]:
            offset_x += sizes["secondary"] + sizes["margin"]
        if sizes["show_axis"]:
            offset_x += sizes["axis"]

        x_transform, y_transform = __get_transform(
            session, sizes["heatmap"], offset_x, sizes["heatmap"], offset_x
        )
        active_anno_y = [
            sizes["annotation"],
            offset,
            session.get_displayed_annos(False, print_callback),
        ]

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
                [x + 0.5 for x in range(len(active_anno_y[2]))],
                [len(active_anno_y[2]), 0, sizes["annotation"], offset],
                offset_x - sizes["axis"],
                offset_x - sizes["spacing"],
                True,
                labels=active_anno_y[2][::-1],
            )
            d.append(
                drawSvg.Line(
                    offset,
                    offset_x - sizes["spacing"],
                    offset + sizes["annotation"],
                    offset_x - sizes["spacing"],
                    stroke="black",
                    stroke_width=2,
                )
            )


def __draw_secondary(session, d, sizes, print_callback=lambda s: None):
    stroke_width = session.get_value(
        ["settings", "export", "secondary_stroke_width", "val"]
    )

    if sizes["show_secondary_y"]:
        offset = 0
        if sizes["show_coords_y"]:
            offset += sizes["coords"]
        if sizes["show_contigs_y"]:
            offset += sizes["contigs"]

        offset_x = 0
        if sizes["show_coords_x"]:
            offset_x += sizes["coords"]
        if sizes["show_contigs_x"]:
            offset_x += sizes["contigs"]
        if sizes["show_anno_x"]:
            offset_x += sizes["annotation"] + sizes["margin"]
        if sizes["show_secondary_x"]:
            offset_x += sizes["secondary"] + sizes["margin"]
        if sizes["show_axis"]:
            offset_x += sizes["axis"]

        if sizes["show_anno_y"]:
            offset += sizes["annotation"] + sizes["margin"]

        x_transform, y_transform = __get_transform(
            session, sizes["heatmap"], offset_x, sizes["heatmap"], offset_x
        )
        min_x, max_x = session.get_min_max_tracks(True, print_callback)
        if sizes["secondary_x_start"] != sizes["secondary_x_end"]:
            min_x = sizes["secondary_x_start"]
            max_x = sizes["secondary_x_end"]
        min_y, max_y = session.get_min_max_tracks(False, print_callback)
        if sizes["secondary_y_start"] != sizes["secondary_y_end"]:
            min_y = sizes["secondary_y_start"]
            max_y = sizes["secondary_y_end"]
        active_anno_x = [
            max_x - min_x,
            min_x,
            sizes["secondary"] - stroke_width,
            offset + stroke_width / 2,
        ]
        active_anno_y = [
            max_y - min_y,
            min_y,
            sizes["secondary"] - stroke_width,
            offset + stroke_width / 2,
        ]

        tracks = session.get_tracks(True, print_callback)
        __draw_lines(
            d,
            tracks,
            x_transform,
            active_anno_x,
            x="screen_pos",
            y="values",
            color="colors",
            stroke_width=stroke_width,
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
        if sizes["show_axis"]:
            __adaptive_ticker(
                d,
                active_anno_x,
                offset_x - sizes["axis"],
                offset_x - sizes["spacing"],
                False,
                label_major=True,
            )
            d.append(
                drawSvg.Line(
                    offset_x - sizes["spacing"],
                    offset,
                    offset_x - sizes["spacing"],
                    offset + sizes["annotation"],
                    stroke="black",
                    stroke_width=2,
                )
            )

    if sizes["show_secondary_x"]:
        offset = 0
        if sizes["show_coords_x"]:
            offset += sizes["coords"]
        if sizes["show_contigs_x"]:
            offset += sizes["contigs"]

        offset_x = 0
        if sizes["show_coords_y"]:
            offset_x += sizes["coords"]
        if sizes["show_contigs_y"]:
            offset_x += sizes["contigs"]
        if sizes["show_anno_y"]:
            offset_x += sizes["annotation"] + sizes["margin"]
        if sizes["show_secondary_y"]:
            offset_x += sizes["secondary"] + sizes["margin"]
        if sizes["show_axis"]:
            offset_x += sizes["axis"]

        if sizes["show_anno_x"]:
            offset += sizes["annotation"] + sizes["margin"]

        x_transform, y_transform = __get_transform(
            session, sizes["heatmap"], offset_x, sizes["heatmap"], offset_x
        )
        min_x, max_x = session.get_min_max_tracks(True, print_callback)
        if sizes["secondary_x_start"] != sizes["secondary_x_end"]:
            min_x = sizes["secondary_x_start"]
            max_x = sizes["secondary_x_end"]
        min_y, max_y = session.get_min_max_tracks(False, print_callback)
        if sizes["secondary_y_start"] != sizes["secondary_y_end"]:
            min_y = sizes["secondary_y_start"]
            max_y = sizes["secondary_y_end"]
        active_anno_x = [
            max_x - min_x,
            min_x,
            sizes["secondary"] - stroke_width,
            offset + stroke_width / 2,
        ]
        active_anno_y = [
            max_y - min_y,
            min_y,
            sizes["secondary"] - stroke_width,
            offset + stroke_width / 2,
        ]

        __draw_lines(
            d,
            session.get_tracks(False, print_callback),
            active_anno_y,
            y_transform,
            y="screen_pos",
            x="values",
            color="colors",
            stroke_width=stroke_width,
        )

        if sizes["show_grid_lines"]:
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
                active_anno_y,
                offset_x - sizes["axis"],
                offset_x - sizes["spacing"],
                True,
                label_major=True,
            )
            d.append(
                drawSvg.Line(
                    offset,
                    offset_x - sizes["spacing"],
                    offset + sizes["annotation"],
                    offset_x - sizes["spacing"],
                    stroke="black",
                    stroke_width=2,
                )
            )


def __draw_coordinates(session, d, sizes, print_callback=lambda s: None):
    offset_x = 0
    if sizes["show_contigs_x"]:
        offset_x += sizes["contigs"]
    offset_y = 0
    if sizes["show_contigs_y"]:
        offset_y += sizes["contigs"]

    offset_heat_x = 0
    offset_heat_y = 0
    if sizes["show_anno_x"]:
        offset_heat_x += sizes["annotation"] + sizes["margin"]
    if sizes["show_secondary_x"]:
        offset_heat_x += sizes["secondary"] + sizes["margin"]
    if sizes["show_coords_x"]:
        offset_heat_x += sizes["coords"]
    if sizes["show_contigs_x"]:
        offset_heat_x += sizes["contigs"]
    if sizes["show_axis"] and (sizes["show_secondary_y"] or sizes["show_anno_y"]):
        offset_heat_x += sizes["axis"]

    if sizes["show_anno_y"]:
        offset_heat_y += sizes["annotation"] + sizes["margin"]
    if sizes["show_secondary_y"]:
        offset_heat_y += sizes["secondary"] + sizes["margin"]
    if sizes["show_coords_y"]:
        offset_heat_y += sizes["coords"]
    if sizes["show_contigs_y"]:
        offset_heat_y += sizes["contigs"]
    if sizes["show_axis"] and (sizes["show_secondary_x"] or sizes["show_anno_x"]):
        offset_heat_y += sizes["axis"]

    x_transform, y_transform = __get_transform(
        session, sizes["heatmap"], offset_heat_x, sizes["heatmap"], offset_heat_y
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

    if sizes["show_coords_y"]:
        __adaptive_ticker(
            d,
            x_transform,
            offset_y,
            offset_y + sizes["coords"] - sizes["spacing"],
            True,
            to_readable_pos=to_readable_pos_x,
            label_major=True,
        )
        d.append(
            drawSvg.Line(
                offset_heat_x,
                offset_y + sizes["coords"] - sizes["spacing"],
                offset_heat_x + sizes["heatmap"],
                offset_y + sizes["coords"] - sizes["spacing"],
                stroke="black",
                stroke_width=2,
            )
        )
    if sizes["show_coords_x"]:
        __adaptive_ticker(
            d,
            y_transform,
            offset_x,
            offset_x + sizes["coords"] - sizes["spacing"],
            False,
            to_readable_pos=to_readable_pos_y,
            label_major=True,
        )
        d.append(
            drawSvg.Line(
                offset_x + sizes["coords"] - sizes["spacing"],
                offset_heat_y,
                offset_x + sizes["coords"] - sizes["spacing"],
                offset_heat_y + sizes["heatmap"],
                stroke="black",
                stroke_width=2,
            )
        )


def minmax(v, mi, ma):
    return min(max(v, mi), ma)


def __draw_contigs(session, d, sizes, print_callback=lambda s: None):
    offset = 0

    offset_heat_x = 0
    offset_heat_y = 0
    if sizes["show_anno_x"]:
        offset_heat_x += sizes["annotation"] + sizes["margin"]
    if sizes["show_secondary_x"]:
        offset_heat_x += sizes["secondary"] + sizes["margin"]
    if sizes["show_coords_x"]:
        offset_heat_x += sizes["coords"]
    if sizes["show_contigs_x"]:
        offset_heat_x += sizes["contigs"]
    if sizes["show_axis"] and (sizes["show_secondary_y"] or sizes["show_anno_y"]):
        offset_heat_x += sizes["axis"]

    if sizes["show_anno_y"]:
        offset_heat_y += sizes["annotation"] + sizes["margin"]
    if sizes["show_secondary_y"]:
        offset_heat_y += sizes["secondary"] + sizes["margin"]
    if sizes["show_coords_y"]:
        offset_heat_y += sizes["coords"]
    if sizes["show_contigs_y"]:
        offset_heat_y += sizes["contigs"]
    if sizes["show_axis"] and (sizes["show_secondary_x"] or sizes["show_anno_x"]):
        offset_heat_y += sizes["axis"]

    x_transform, y_transform = __get_transform(
        session, sizes["heatmap"], offset_heat_x, sizes["heatmap"], offset_heat_y
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

    y_label, x_label = session.get_value(
        ["settings", "interface", "axis_lables"]
    ).split("_")

    if sizes["show_contigs_y"]:
        __draw_tick_lines(
            d,
            contig_starts_x,
            x_transform,
            offset + sizes["contigs"] - 8,
            offset + sizes["contigs"] - sizes["spacing"],
            True,
            labels=None,
        )
        d.append(
            drawSvg.Line(
                offset_heat_x,
                offset + sizes["contigs"] - sizes["spacing"],
                offset_heat_x + sizes["heatmap"],
                offset + sizes["contigs"] - sizes["spacing"],
                stroke="black",
                stroke_width=2,
            )
        )

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

        d.append(
            drawSvg.Text(
                x_label,
                18,
                offset_heat_x + sizes["heatmap"] / 2,
                offset,
                font_family="Consolas, sans-serif",
                text_anchor="middle",
                dominant_baseline="bottom",
            )
        )

    if sizes["show_contigs_x"]:
        __draw_tick_lines(
            d,
            contig_starts_y,
            y_transform,
            offset + sizes["contigs"] - 8,
            offset + sizes["contigs"] - sizes["spacing"],
            False,
            labels=None,
        )
        d.append(
            drawSvg.Line(
                offset + sizes["contigs"] - sizes["spacing"],
                offset_heat_y,
                offset + sizes["contigs"] - sizes["spacing"],
                offset_heat_y + sizes["heatmap"],
                stroke="black",
                stroke_width=2,
            )
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

        d.append(
            drawSvg.Text(
                y_label,
                18,
                offset,
                offset_heat_y + sizes["heatmap"] / 2,
                font_family="Consolas, sans-serif",
                transform="rotate(-90,"
                + str(offset)
                + ","
                + str(-(offset_heat_y + sizes["heatmap"] / 2))
                + ")",
                text_anchor="middle",
                dominant_baseline="hanging",
            )
        )


def __get_sizes(session):
    return {
        "show_heat": session.get_value(
            ["settings", "interface", "show_hide", "heatmap"]
        ),
        "show_region": session.get_value(["settings", "export", "print_region"]),
        "show_coords_x": session.get_value(
            ["settings", "interface", "show_hide", "coords"]
        )
        and session.get_value(["settings", "interface", "show_hide", "heatmap"]),
        "show_coords_y": session.get_value(
            ["settings", "interface", "show_hide", "coords"]
        ),
        "show_contigs_x": session.get_value(
            ["settings", "interface", "show_hide", "regs"]
        )
        and session.get_value(["settings", "interface", "show_hide", "heatmap"]),
        "show_contigs_y": session.get_value(
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
        "show_anno_x": session.get_value(
            ["settings", "interface", "show_hide", "annotation"]
        )
        and session.get_value(["settings", "interface", "show_hide", "heatmap"])
        and len(session.get_annotation(False, lambda x: None)["anno_name"]) > 0,
        "show_anno_y": session.get_value(
            ["settings", "interface", "show_hide", "annotation"]
        )
        and len(session.get_annotation(True, lambda x: None)["anno_name"]) > 0,
        "annotation": session.get_value(["settings", "interface", "anno_size", "val"]),
        "show_secondary_x": session.get_value(
            ["settings", "interface", "show_hide", "raw"]
        )
        and session.get_value(["settings", "interface", "show_hide", "heatmap"])
        and len(session.get_tracks(False, lambda x: None)["values"]) > 0,
        "show_secondary_y": session.get_value(
            ["settings", "interface", "show_hide", "raw"]
        )
        and len(session.get_tracks(True, lambda x: None)["values"]) > 0,
        "secondary": session.get_value(["settings", "interface", "raw_size", "val"]),
        "spacing": session.get_value(["settings", "export", "spacing", "val"]),
        "white_background": session.get_value(
            ["settings", "export", "white_background"]
        ),
        "secondary_x_start": session.get_value(
            ["settings", "export", "secondary_x_range", "val_min"]
        ),
        "secondary_x_end": session.get_value(
            ["settings", "export", "secondary_x_range", "val_max"]
        ),
        "secondary_y_start": session.get_value(
            ["settings", "export", "secondary_y_range", "val_min"]
        ),
        "secondary_y_end": session.get_value(
            ["settings", "export", "secondary_y_range", "val_max"]
        ),
    }


def __make_drawing(session, sizes):
    size_x = 0
    size_y = 0

    size_x += sizes["heatmap"] + sizes["margin"]
    if sizes["show_anno_x"]:
        size_x += sizes["annotation"] + sizes["margin"]
    if sizes["show_secondary_x"]:
        size_x += sizes["secondary"] + sizes["margin"]
    if sizes["show_coords_x"]:
        size_x += sizes["coords"]
    if sizes["show_contigs_x"]:
        size_x += sizes["contigs"]
    if sizes["show_axis"] and (sizes["show_secondary_y"] or sizes["show_anno_y"]):
        size_x += sizes["axis"]

    size_x -= sizes["margin"]

    if sizes["show_heat"]:
        size_y += sizes["heatmap"] + sizes["margin"]
    if sizes["show_anno_y"]:
        size_y += sizes["annotation"] + sizes["margin"]
    if sizes["show_secondary_y"]:
        size_y += sizes["secondary"] + sizes["margin"]
    if sizes["show_coords_y"]:
        size_y += sizes["coords"]
    if sizes["show_contigs_y"]:
        size_y += sizes["contigs"]
    if sizes["show_axis"] and (sizes["show_secondary_x"] or sizes["show_anno_x"]):
        size_y += sizes["axis"]
    if sizes["show_region"]:
        size_y += 16 * 14 / 12 + sizes["margin"]

    size_y -= sizes["margin"]

    d = drawSvg.Drawing(size_x, size_y, displayInline=False)

    if sizes["white_background"]:
        d.append(
            drawSvg.Rectangle(
                0,
                0,
                size_x,
                size_y,
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
    __draw_region(session, d, sizes)
    return d


def assert_has_draw_svg():
    if not HAS_DRAW_SVG:
        raise RuntimeError(
            "could not import drawSvg, is it installed? (try pip install drawSvg==1.9.0)"
        )


def assert_file_is_creatable(session):
    try:
        with open(session.get_value(["settings", "export", "prefix"]), "w") as _:
            pass
        os.remove(session.get_value(["settings", "export", "prefix"]))
    except:
        raise RuntimeError("could not create file, does the given folder exist?")


def get_tsv(session, print_callback=lambda s: None):
    session = __conf_session(session)
    heatmap, tracks = "", ["", ""]

    heatmap = __get_header(session)
    heatmap += "#chr_x\tstart_x\tend_x\tchr_y\tstart_y\tend_y\tscore\n"
    for tup in session.get_heatmap_export(print_callback):
        heatmap += "\t".join([str(x) for x in tup]) + "\n"

    if session.get_value(["settings", "interface", "show_hide", "raw"]):
        for x_axis in [(True), (False)]:
            tracks[x_axis] = __get_header(session)
            tracks[x_axis] += "#chr_x\tstart_x\tend_x"
            for track in session.get_track_export_names(x_axis, print_callback):
                tracks[x_axis] += "\t" + track
            tracks[x_axis] += "\n"

            for tup in session.get_track_export(x_axis, print_callback):
                tracks[x_axis] += (
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
    return heatmap, tracks[0], tracks[1]


def export_tsv(session, print_callback=lambda s: None):
    assert_file_is_creatable(session)
    heatmap, track_x, track_y = get_tsv(session, print_callback)
    with open(
        session.get_value(["settings", "export", "prefix"]) + ".heatmap.tsv", "w"
    ) as out_file:
        out_file.write(heatmap)
    with open(
        session.get_value(["settings", "export", "prefix"]) + ".track.x.tsv", "w"
    ) as out_file:
        out_file.write(track_x)
    with open(
        session.get_value(["settings", "export", "prefix"]) + ".track.y.tsv", "w"
    ) as out_file:
        out_file.write(track_y)


def get_png(session):
    assert_has_draw_svg()
    return __draw(session).rasterize().pngData


def export_png(session):
    assert_file_is_creatable(session)
    data = get_png(session)
    with open(
        session.get_value(["settings", "export", "prefix"]) + ".png", "wb"
    ) as out_file:
        out_file.write(data)


def get_svg(session):
    assert_has_draw_svg()
    return __draw(session).asSvg()


def export_svg(session):
    assert_file_is_creatable(session)
    data = get_svg(session)

    with open(
        session.get_value(["settings", "export", "prefix"]) + ".svg",
        "w",
        encoding="utf-8",
    ) as out_file:
        out_file.write(data)
