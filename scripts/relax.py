import os
import yaml
import sys
import logging
def check_package(package_name,import_name):
    from importlib.util import find_spec

    has_lib = find_spec(package_name)
    if not has_lib:
        print("%s is not installed!!"%(import_name))
        return False

       

if __name__ == "__main__":
    package_list = {"biopython":"Bio","pandas":"pandas"}
    config_path = os.path.join(os.path.dirname(sys.argv[0]),"config.yaml")
    with open(config_path,"r") as f:
        config = yaml.load(f,Loader=yaml.FullLoader)
    #relax
    if config["input"]["step"] == "relax":
        #score_jd2
        if config["relax"]["need_prepare"]:
            cmd = os.path.join(os.path.dirname(sys.argv[0]),"jd2.run.sh") + " " +\
                config["input"]["path"] +" " +\
                os.path.join(config["input"]["rosetta_bin"],config["relax"]["score_jd2_exec"])
            ex = os.system(cmd)
            if ex:
                logging.error("run score_jd2 failed !")
                sys.exit(1)
            jd2_file = os.path.abspath(config["input"]["path"]) + "/jd2_pdb.list"
        else:
            jd2_file = config["input"]["path"]
        #if mpi_run relax
        if config["relax"]["mpi_run"]["use_mpi"]:
            mpi_cores = str(config["relax"]["mpi_run"]["mpi_cores"])
        else:
            mpi_cores = "0"
        #run relax
        cmd = os.path.join(os.path.dirname(sys.argv[0]),"relax.run.sh") + " " +\
            str(mpi_cores) + " " +\
            str(config["relax"]["relax_threads"])+ " " +\
            os.path.abspath(os.path.join(os.path.dirname(sys.argv[0]),config["relax"]["flag"])) + " " +\
            os.path.join(config["input"]["rosetta_bin"],config["relax"]["relax_exec"]) + " " +\
            jd2_file + " " +\
            config["input"]["path"]

        ex = os.system(cmd)
        if ex:
            logging.error("run relax failed !")
            sys.exit(1)
        relax_file = os.path.abspath(config["input"]["path"]) + "/relax_sc.list"
        could_flex = True
    elif config["input"]["step"] == "flex_ddg":
        relax_file = "None"
        could_flex = True
    else:
        logging.error("'step' not defined correctly !")
        sys.exit(2)
    
        

    