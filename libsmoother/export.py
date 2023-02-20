from .quarry import Quarry
import drawSvg


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


def __conv_coords(p, is_bottom, idx, cds, w_cds, o_cds, w_plane, o_plane):
    return max(min(w_plane * (p - o_cds) / w_cds + o_plane, o_plane + w_plane), o_plane)


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
            [conv_coords(p, is_bottom, *transform) for p in ps]
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


def __draw_annotation(session, d, sizes, print_callback=lambda s: None):
    if sizes["show_anno"]:
        offset = 0
        if sizes["show_secondary"]:
            offset += sizes["secondary"] + sizes["margin"]

        offset_x = offset + sizes["annotation"] + sizes["margin"]

        d.append(
            drawSvg.Rectangle(
                offset_x, offset, sizes["heatmap"], sizes["annotation"], fill="#ffffff"
            )
        )
        d.append(
            drawSvg.Rectangle(
                offset, offset_x, sizes["annotation"], sizes["heatmap"], fill="#ffffff"
            )
        )

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


def __draw_secondary(session, d, sizes, print_callback=lambda s: None):
    if sizes["show_secondary"]:
        offset = 0

        anno_offset = 0
        if sizes["show_anno"]:
            anno_offset += sizes["annotation"] + sizes["margin"]

        offset_x = offset + anno_offset + sizes["secondary"] + sizes["margin"]

        d.append(
            drawSvg.Rectangle(
                offset_x, offset, sizes["heatmap"], sizes["secondary"], fill="#ffffff"
            )
        )
        d.append(
            drawSvg.Rectangle(
                offset, offset_x, sizes["secondary"], sizes["heatmap"], fill="#ffffff"
            )
        )

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


def __get_sizes(session):
    return {
        "show_heat": True,
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

    size -= sizes["margin"]

    return drawSvg.Drawing(size, size, displayInline=False)


def __draw(session):
    session = __conf_session(session)
    sizes = __get_sizes(session)
    d = __make_drawing(session, sizes)
    __draw_heatmap(session, d, sizes)
    __draw_annotation(session, d, sizes)
    __draw_secondary(session, d, sizes)
    return d


def export_png(session):
    __draw(session).savePng(
        session.get_value(["settings", "export", "prefix"]) + ".png"
    )


def export_svg(session):
    __draw(session).saveSvg(
        session.get_value(["settings", "export", "prefix"]) + ".svg"
    )
