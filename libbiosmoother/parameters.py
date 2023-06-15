try:
    import importlib.resources as pkg_resources
except ImportError:
    # Try backported to PY<37 `importlib_resources`.
    import importlib_resources as pkg_resources


def json_get(keys, d):
    t = d
    for k in keys:
        t = t[k]
    return t


def __exists(keys, d):
    t = d
    for k in keys:
        if not k in t:
            return False
        t = t[k]
    return True


def __is_checkbox(val):
    return isinstance(val, bool)


def is_spinner(val):
    return isinstance(val, dict) and sorted(list(val.keys())) == sorted(
        [
            "min",
            "max",
            "val",
            "step",
            "spinner_max_restricted",
            "spinner_min_restricted",
        ]
    )


def is_range_spinner(val):
    return isinstance(val, dict) and sorted(list(val.keys())) == sorted(
        [
            "min",
            "max",
            "val_max",
            "val_min",
            "step",
            "spinner_max_restricted",
            "spinner_min_restricted",
        ]
    )


def __is_multi_choice(key, val, valid):
    return (
        isinstance(val, str)
        and __exists(key, valid)
        and isinstance(json_get(key, valid), list)
    )


def spinner_is_int(d):
    return (
        isinstance(d["min"], int)
        and isinstance(d["max"], int)
        and isinstance(d["step"], int)
    )


def __is_parameter(key, conf, valid):
    val = json_get(key, conf)
    if is_spinner(val):
        return True
    if is_range_spinner(val):
        return True
    if __is_checkbox(val):
        return True
    if __is_multi_choice(key, val, valid):
        return True
    return False


def list_parameters(conf, valid, curr_name=[]):
    for key in json_get(curr_name, conf).keys():
        if key in ["previous", "next"]:
            continue
        curr_key = curr_name + [key]
        if __is_parameter(curr_key, conf, valid):
            yield curr_name + [key]
        elif isinstance(json_get(curr_key, conf), dict):
            yield from list_parameters(conf, valid, curr_name=curr_key)


def values_for_parameter(key, conf, valid):
    val = json_get(key, conf)
    if __is_checkbox(val):
        return [True, False]
    if is_spinner(val):
        c = (val["min"] + val["max"]) / 2
        if spinner_is_int(val):
            c = int(c)
        return [val["min"], c, val["max"]]
    if is_range_spinner(val):
        c = (val["min"] + val["max"]) / 2
        if spinner_is_int(val):
            c = int(c)
        return [
            (val["min"], val["min"]),  # zero range
            (val["min"], val["max"]),  # full range
            (val["min"], c),  # below max
            (c, val["max"]),  # above min
        ]
    if __is_multi_choice(key, val, valid):
        return json_get(key, valid)


def open_valid_json():
    return (pkg_resources.files("libbiosmoother") / "conf" / "valid.json").open("r")
