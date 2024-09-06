from mutate_model_class_latest import *

get_input = input("Please enter raw pdb and mutant name: ").split()
pdb = get_input[0]
make = "/".join(pdb.split("/")[:-1])
mut = get_input[1].split("/")[-1][:-8].split("_")
mut[0] = f"{mut[0]}/{mut[0]}"
try:
    new_model = MutateModel(*mut)
    MutateModel.mutate_model(new_model)
except Exception as error:
    record = f"{mut[0]}---{mut[1]}---{mut[2]}---{mut[3]}--{type(error).__name__}--{error}"
    newfile = open(f"{make}/error_msg_{mut[0]}---{mut[1]}---{mut[2]}---{mut[3]}--{type(error).__name__}--{error}.txt", 'w', newline="")
    newfile.seek(0, 0)
    newfile.write(record)
    newfile.close()
