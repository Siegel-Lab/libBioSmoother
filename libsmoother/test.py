from .parameters import list_parameters, values_for_parameter, open_valid_json, json_get, is_spinner, is_range_spinner, spinner_is_int
import json
from .quarry import open_default_json
import traceback
import random
import sys

def __test_config(quarry, name, idx, skip_first):
    if idx >= skip_first:
        try:
            quarry.clear_cache()
            print("[" + str(idx) + "]", name)
            quarry.update_all(lambda s: None)
        except Exception:
            traceback.print_exc()
    else:
        print("[" + str(idx) + "]", name, "- skipped")

    return idx + 1


def __configure(quarry, conf):
    quarry.set_value(["settings"], conf)


def __list_configurations():
    with open_default_json() as default_file:
        default_json = json.load(default_file)
    return [
        (default_json["display_name"], default_json), 
    ]

def __vary_params_one_by_one(quarry, default_json, valid_json, idx, skip_first):
    print("- Varied parameters -")
    for p in list_parameters(default_json, valid_json):
        default_val = quarry.get_value(p)
        for v in values_for_parameter(p, default_json, valid_json):
            if v == default_val:
                continue
            if is_spinner(json_get(p, default_json)):
                quarry.set_value(["settings"] + p + ["val"], v)
            elif is_range_spinner(json_get(p, default_json)):
                quarry.set_value(["settings"] + p + ["val_min"], v[0])
                quarry.set_value(["settings"] + p + ["val_max"], v[1])
            else:
                quarry.set_value(["settings"] + p, v)
            idx = __test_config(quarry, " ".join(["set", ".".join(p), "to", str(v)]), idx, skip_first)
        quarry.set_value(p, default_val)

def __test_default_configs(quarry, default_json, valid_json, idx, skip_first):
    for name, conf in __list_configurations():
        print("-- Testing configuration:", name, "--")
        __configure(quarry, conf)
        idx = __test_config(quarry, "Default run", idx, skip_first)
        idx = __vary_params_one_by_one(quarry, default_json, valid_json, idx, skip_first)
    return idx


def __config_randomly(quarry, default_json, valid_json):
    __configure(quarry, default_json)
    for p in list_parameters(default_json, valid_json):
        d = json_get(p, default_json)
        if is_spinner(d):
            min_, max_ = json_get(p + ["min"], default_json), json_get(p + ["max"], default_json)
            if spinner_is_int(d):
                v = random.randint(min_, max_)
            else:
                v = random.uniform(min_, max_)
            quarry.set_value(["settings"] + p + ["val"], v)
        elif is_range_spinner(d):
            min_, max_ = json_get(p + ["min"], default_json), json_get(p + ["max"], default_json)
            if spinner_is_int(d):
                v = sorted([random.randint(min_, max_), random.randint(min_, max_)])
            else:
                v = sorted([random.uniform(min_, max_), random.uniform(min_, max_)])
            quarry.set_value(["settings"] + p + ["val_min"], v[0])
            quarry.set_value(["settings"] + p + ["val_max"], v[1])
        else:
            v = list(values_for_parameter(p, default_json, valid_json))
            quarry.set_value(["settings"] + p, random.choice(v))
    print(quarry.get_value(["settings"]))

def __test_random_configs(quarry, default_json, valid_json, idx, skip_first, seed=None):
    if seed is None:
        seed = random.randrange(sys.maxsize)
    print("Seed set to: ", seed)
    conf_seeds = [random.seed(seed) for _ in range(5)]

    for s in conf_seeds:
        random.seed(s)
        print("-- Random configuration --")
        __config_randomly(quarry)
        idx = __test_config(quarry, idx, skip_first)
        idx = __vary_params_one_by_one(quarry, default_json, valid_json, idx, skip_first)
    return idx



def test(quarry, seed, skip_first):
    with open_valid_json() as valid_file:
        valid_json = json.load(valid_file)
    with open_default_json() as default_file:
        default_json = json.load(default_file)
    idx = 0

    idx = __test_default_configs(quarry, default_json, valid_json, idx, skip_first)
    idx = __test_random_configs(quarry, default_json, valid_json, seed, idx, skip_first)

    print("Done.")
