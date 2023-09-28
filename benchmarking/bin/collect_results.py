import os
import pickle
import gzip
import pickle

RESO = ["50000", "10000", "5000"]
SUBSAMPLES = ["1", "0.2", "0.4", "0.6", "0.8"] # "0.1", "0.01", "0.001", "0.0001"
#ORGANISMS = ["drosophila_m", "e_coli", "s_cerevisiae", "t_brucei", "human"]
ORGANISMS = ["drosophila_m", "e_coli", "s_cerevisiae", "t_brucei"]
ORGA_SUBSAMPLE = ["t_brucei"]
#SUBSAMPLE_ORGA = ["even_1", "even_2"]
SUBSAMPLE_ORGA = ["even_1"]
NO_PARAMS = [""] #, "-c.-m", "-q.-c", "-q.-c.-m"]
PARAMS = ["", "-q", "-c", "-m", "-q.-c", "-q.-m", "-c.-m", "-q.-c.-m"]

def read_line_from_file(filename, line_num):
    with open(filename, 'r') as f:
        return list(f.readlines())[line_num].strip()
    
def read_line_col_from_file(filename, line_num, col_num):
    return read_line_from_file(filename, line_num).split()[col_num]

data = {}

def get_data(organisms, subsamples, params, resos):
    error_or_in_progress_count = 0
    for organism in organisms:
        print("organism:", organism)
        if not organism in data:
            data[organism] = {}
            
        unique_interactions = read_line_from_file("organisms/" + organism + "/out/unique_interactions_per_sq_kb.txt", 0)
        data[organism]["unique_interactions_per_sq_kb"] = unique_interactions
        print("\t", "number of unique interactions per square kb", unique_interactions)
        unique_interactions_post_subs = read_line_from_file("organisms/" + organism + "/out/unique_interactions_per_sq_kb_post_subsample_1.txt", 0)
        data[organism]["unique_interactions_post_subs_1"] = unique_interactions_post_subs
        print("\t", "number of unique interactions per square kb post subsampling even_1", unique_interactions_post_subs)
        unique_interactions_post_subs = read_line_from_file("organisms/" + organism + "/out/unique_interactions_per_sq_kb_post_subsample_2.txt", 0)
        data[organism]["unique_interactions_post_subs_2"] = unique_interactions_post_subs
        print("\t", "number of unique interactions per square kb post subsampling even_2", unique_interactions_post_subs)

        genome_size = 0
        for filename in os.listdir("organisms/" + organism + "/genome"):
            if filename.endswith(".sizes"):
                with open("organisms/" + organism + "/genome/" + filename, 'r') as f:
                    for line in f:
                        genome_size += int(line.split()[1])
                    data[organism]["genome_size"] = genome_size
                    print("\t", "genome size", genome_size)
        assert genome_size > 0

        for subsample in subsamples:
            print("\t", "subsample:", subsample)
            if not subsample in data[organism]:
                data[organism][subsample] = {}

            int_cnt_file = "organisms/" + organism + "/out/interaction_count." + subsample + ".txt"
            num_interactions = int(read_line_col_from_file(int_cnt_file, 0, 0))
            num_unique_interactions = int(read_line_col_from_file(int_cnt_file, 1, 0))
            data[organism][subsample]["num_interactions"] = num_interactions
            data[organism][subsample]["num_unique_interactions"] = num_unique_interactions
            print("\t", "\t", "number of interactions", num_interactions)
            print("\t", "\t", "number of unique interactions", num_unique_interactions)
            #assert num_interactions > 0
            for param in params:
                print("\t", "\t", "param:", param)
                if not param in data[organism][subsample]:
                    data[organism][subsample][param] = {}


                for reso in resos:
                    print("\t", "\t", "\t", "resolution:", reso)
                    data[organism][subsample][param][reso] = {}

                    build_status_file = "./organisms/" + organism + "/out/index_build_status." + \
                                                    reso + "." + subsample + param + ".txt"
                    if os.path.isfile(build_status_file):
                        build_status = read_line_from_file(build_status_file, 0)
                    else:
                        build_status = "not started"
                    data[organism][subsample][param][reso]["build_status"] = build_status
                    print("\t", "\t", "\t", "\t", "build status:", build_status)
                    if build_status == "OK":
                        index_size_file = "organisms/" + organism + "/out/index_size." + reso + "." + subsample + param + ".txt"
                        index_size = read_line_col_from_file(index_size_file, 0, 0)
                        index_size_h = read_line_col_from_file(index_size_file, 1, 0)
                        print("\t", "\t", "\t", "\t", "index size:", index_size, "( =", index_size_h, ")")
                        data[organism][subsample][param][reso]["index_size"] = index_size

                        index_build_time = read_line_col_from_file("organisms/" + organism + "/out/2_index_build_times." + 
                                                                reso + "." + subsample + param + ".txt", 6, 1)
                        print("\t", "\t", "\t", "\t", "index build time:", index_build_time)
                        data[organism][subsample][param][reso]["index_build_time"] = index_build_time

                        index_init_time = read_line_col_from_file("organisms/" + organism + "/out/2_index_build_times." + 
                                                                reso + "." + subsample + param + ".txt", 1, 1)
                        print("\t", "\t", "\t", "\t", "index init time:", index_init_time)
                        data[organism][subsample][param][reso]["index_init_time"] = index_init_time

                        data[organism][subsample][param][reso]["query_times"] = {}
                        all_rep_ok = True
                        for rep_int in range(100):
                            rep_str = str(rep_int)
                            query_status_file = "organisms/" + organism + "/out/index_query_times." + \
                                                        reso + "." + subsample + param + "." + rep_str + ".txt"
                            pickle_file = "organisms/" + organism + "/out/benchmark." + \
                                            reso + "." + subsample + param + "." + rep_str + ".pickle"
                            if os.path.isfile(query_status_file) and os.path.isfile(pickle_file):
                                query_status = read_line_from_file(query_status_file, 0)
                                #query_status = "broken"
                            else:
                                query_status = "not started"
                            #print("\t", "\t", "\t", "\t", "query_status", rep_int, query_status)
                            if query_status == "OK":
                                with open(pickle_file, "rb") as f:
                                    data[organism][subsample][param][reso]["query_times"][rep_int] = pickle.load(f)
                            else:
                                all_rep_ok = False
                        if all_rep_ok:
                            query_status = "OK"
                        else:
                            query_status = "Error or in progress"
                            error_or_in_progress_count += 1
                        data[organism][subsample][param][reso]["query_status"] = query_status
                        print("\t", "\t", "\t", "\t", "query_status", query_status)
                    else:
                        error_or_in_progress_count += 1

    print(error_or_in_progress_count, "tasks are unfinished")
    print()

get_data(ORGANISMS, SUBSAMPLE_ORGA, NO_PARAMS, RESO)
get_data(ORGA_SUBSAMPLE, SUBSAMPLES, NO_PARAMS, RESO)
get_data(ORGA_SUBSAMPLE, SUBSAMPLE_ORGA, PARAMS, RESO)

with open("index_sizes_data.pickle", "wb") as out_file:
    pickle.dump(data, out_file)