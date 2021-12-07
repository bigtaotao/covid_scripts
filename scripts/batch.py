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
    try:
        config_file = sys.argv[1]
        print("load config %s"%(config_file))
        config_path = (config_file)
    except:
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
            print("========= start running score_jd2 =========")
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
        print("========= start running relax =========")
        ex = os.system(cmd)
        if ex:
            logging.error("run relax failed !")
            sys.exit(1)
        relax_file = os.path.abspath(os.path.join(config["input"]["path"],"relax_sc.list"))
        could_flex = True
    elif config["input"]["step"] == "flex_ddg":
        relax_file = os.path.abspath(os.path.join(config["input"]["path"],"relax_sc.list"))
        if not os.path.exists(relax_file):
            relax_file = "None"
        could_flex = True
    else:
        logging.error("'step' not defined correctly !")
        sys.exit(2)
    #run flex_ddG
    if could_flex:
        import os
        from typing import SupportsAbs

        import sys
        from Bio.PDB import PDBParser
        import pandas as pd
        import copy
        import shutil
        import yaml
        import re
        print("========= start running make fles_ddG =========")
        relax_res = relax_file
        # relax_res = sys.argv[1]
        meta = pd.read_csv(config["input"]["path"] + "/meta.csv")

        to_csv = os.path.abspath(config["input"]["path"]) + "/flex_input.csv"
        input_path_root = "temp"
        mut_list = config["input"]["mutation"]
        def filter_name(name):
            r = re.search(r"(_[0-9]{4})+",name)
            try:
                name = name[0:r.start()]
            except:
                name = name.split(".")[0]
            return(name)
        if not relax_res == "None":
            with open(relax_res,"r") as f:
                file = f.read().splitlines()
            print("find %d score files"%(len(file)))


            result = []

            for each in file :
                with open(each,"r") as f:
                    res = f.read().splitlines()
                split_res = []
                for i in res:
                    split_res.append(i.split())
                split_res.pop(0)
                split_res.pop(0)
                split_pd  = pd.DataFrame(split_res)

                def add_name (x):
                    # return(x.iloc[21].split(".")[0]) # renum
                    
                    return(filter_name(x.iloc[21])) # no renum

                split_pd.loc[:,"name"] = split_pd.apply(add_name,axis = 1)

                name_list = list(set(split_pd.name.to_list()))
                print(str(name_list) + "ok")
                for each_name in name_list:
                    subset = split_pd.loc[split_pd.name == each_name]
                    min_idx = subset.iloc[:,1].astype("float").idxmin()
                    min_file = subset.at[min_idx,21]
                    # lig = min_file.split(".")[0]
                    # print(lig)
                    # lig = meta.loc[meta.Name == lig,"Column"].to_list()[0]

                    dirpath = os.path.dirname(each)
                    res_path = os.path.join(dirpath,min_file)
                    res_path = res_path + ".pdb"
                    

                    # add = res_path + "&" + lig

                    result.append(res_path)
        else:
            cmd = "find " + config["input"]["path"] + ' -name "*.pdb"'
            result = os.popen(cmd).read().splitlines()

        #In[]
        # print(str(result))
        # print(len(result))
        file_list = result
        content = []
        for each_pdb in file_list:
            each = each_pdb
            # name = os.path.basename(each).split(".")[0] # renum
            name = filter_name( os.path.basename(each)) #no renum
            print(name)
            chain = meta.loc[meta.Name == name,"rbd"].to_list()[0]
            content.append([each,name,chain])
            
        resi_dict = {"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E","GLY":"G","HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P","SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V"}

        # %%
        merge_res = pd.DataFrame()
        next_step = []
        ##find which have na res
        have_na_list = []
        for each_pdb in content:
            pdbparser=PDBParser()

            path = each_pdb[0]
            name = each_pdb[1]
            chain = each_pdb[2]
            pdb= pdbparser.get_structure(name,path)
            p = pdb[0][chain]

            seq = p.child_list[0].get_id()[1]

            need_mut = []

            for each in mut_list:
                try:

                    cache = resi_dict[p[mut_list[each][0]].get_resname()]
                    if cache == mut_list[each][1]:
                        res.append(each)
                    else:
                        need_mut.append(each)
                except:
                    pass

        
            merge_res = merge_res.append({"path":path,"name":name,"need_mut":need_mut},ignore_index = True)
            next_step.append([name,meta.loc[meta.Name == name,"comp"].to_list()[0],path,need_mut,meta.loc[meta.Name == name,"rbd"].to_list()[0]])
            #make init_input
            # next_step.append([name,meta.loc[meta.Name == name,"Column"].to_list()[0],each_pdb[3],need_mut,meta.loc[meta.Name == name,"Tags"].to_list()[0],shift])
        merge_res.to_csv(to_csv)
        # %%
        ##构建flex-ddG input文件夹
        os.makedirs(input_path_root,exist_ok= True)
        sub_dir = config["input"]["combine"]
        # sub_dir = ["single","multi","third"]

        def make_file (subpath,e_pdb,muts):
            mut_tail = ""
            for e_mut in muts:
                mut_tail = mut_tail + "_" + e_mut
            
            path = os.path.join(input_path_root,subpath,e_pdb[0]+mut_tail)
            os.makedirs(path,exist_ok= True)
            shutil.copy(e_pdb[2],path)
            with open(os.path.join(path,"chains_to_move.txt"),"w") as f:
                f.writelines(e_pdb[4])
            with open(os.path.join(path,"nataa_mutations.resfile"),"w") as f:
                f.writelines("NATAA \n")
                f.writelines("start \n")
                for e_mut in muts:
                    # print(e_mut)
                    f.writelines(str(mut_list[e_mut][0]) + " " + e_pdb[4] + " PIKAA " + mut_list[e_mut][1] + " \n") 

        def exhaust_pos(pos_list,num):
            res = []
            for i,pos in enumerate(pos_list):
                if num > 1:
                    cache = exhaust_pos(pos_list[i+1:],num-1)
                    for each in cache:
                        res.append([pos]+each)
                    
                else :
                    cache = []
                    for each in pos_list:
                        cache.append([each])
                    return cache
            return(res)
        # print(exhaust_pos(["a","b","c","d"],4))

        def filter_pos (pos_list):
            res = []
            for each_pos in pos_list:
                cache = []
                for each in each_pos:
                    cache.append(each[1:-1])
                if len(set(cache)) == len(cache):
                    res.append(each_pos)
            return(res)
        # print(filter_pos([['n501y', 'k417n'], ['n501y', 'k417t'], ['n501y', 'e484k'], ['k417n', 'k417t'], ['k417n', 'e484k'], ['k417t', 'e484k']]))
        for each in sub_dir:
            os.makedirs(os.path.join(input_path_root,each),exist_ok= True)
            for each_pdb in next_step:
                if sub_dir[each] == "all":
                    mut_posi = exhaust_pos(each_pdb[3],each_pdb[3])
                elif len(each_pdb[3]) >= sub_dir[each]:
                    mut_posi = exhaust_pos(each_pdb[3],sub_dir[each])
                    mut_posi = filter_pos(mut_posi)
                else:
                    continue
                    #add R  
                for pos in mut_posi:
                    mut_tail = ""
                    for e_mut in pos:
                        mut_tail = mut_tail + "_" + e_mut
                    with open(os.path.abspath(config["input"]["path"]) + "/mut.list","a") as f:
                        f.write(mut_tail[1:] + "\n")
                    make_file(each,each_pdb,pos)

        ## split input

        output_path = config["flex_ddg"]["root_path"] #最后输出的位置
        all_need_file = [
                        os.path.join(os.path.dirname(sys.argv[0]),"ddG-backrub.xml")
                        ]
        def cp_ddg(dst_path):
            global all_need_file
            for each_all_need in all_need_file:
                shutil.copy(each_all_need,dst_path)
            with open(os.path.join(dst_path,"run_flex.py"),"w") as f:
                content = "max_cpus = " + str(config["flex_ddg"]["max_cpus"]) + "\n"
                content += "nstruct = " + str(config["flex_ddg"]["nstruct"]) + "\n"
                content += "max_minimization_iter = " +  str(config["flex_ddg"]["max_minimization_iter"]) + "\n"
                content += "abs_score_convergence_thresh = " + str(config["flex_ddg"]["abs_score_convergence_thresh"]) + "\n"
                content += "number_backrub_trials = " + str(config["flex_ddg"]["number_backrub_trials"]) + "\n"
                content += "backrub_trajectory_stride = " + str(config["flex_ddg"]["backrub_trajectory_stride"]) + "\n"
                content += "rosetta_scripts_path = '" + os.path.join(config["input"]["rosetta_bin"],config["flex_ddg"]["rosetta_scripts_exec"]) + "'\n"
                content = content + r'''
import socket
import sys
import os
import subprocess

use_multiprocessing = True
if use_multiprocessing:
    import multiprocessing
path_to_script = 'ddG-backrub.xml'

if not os.path.isfile(rosetta_scripts_path):
    print('ERROR: "rosetta_scripts_path" variable must be set to the location of the "rosetta_scripts" binary executable')
    print('This file might look something like: "rosetta_scripts.linuxgccrelease"')
    raise Exception('Rosetta scripts missing')

def run_flex_ddg( name, input_path, input_pdb_path, chains_to_move, nstruct_i ):
    output_directory = os.path.join( 'output', os.path.join( name, '%02d' % nstruct_i ) )
    if not os.path.isdir(output_directory):
        os.makedirs(output_directory)

    flex_ddg_args = [
        os.path.abspath(rosetta_scripts_path),
        "-s %s" % os.path.abspath(input_pdb_path),
        '-parser:protocol', os.path.abspath(path_to_script),
        '-parser:script_vars',
        'chainstomove=' + chains_to_move,
        'mutate_resfile_relpath=' + os.path.abspath( os.path.join( input_path, 'nataa_mutations.resfile' ) ),
        'number_backrub_trials=%d' % number_backrub_trials,
        'max_minimization_iter=%d' % max_minimization_iter,
        'abs_score_convergence_thresh=%.1f' % abs_score_convergence_thresh,
        'backrub_trajectory_stride=%d' % backrub_trajectory_stride ,
        '-restore_talaris_behavior',
        '-in:file:fullatom',
        '-ignore_unrecognized_res',
        '-ignore_zero_occupancy false',
        '-ex1',
        '-ex2',
    ]

    log_path = os.path.join(output_directory, 'rosetta.out')

    print( 'Running Rosetta with args:' )
    print( ' '.join(flex_ddg_args) )
    print( 'Output logged to:', os.path.abspath(log_path) )
    print()

    outfile = open(log_path, 'w')
    process = subprocess.Popen(flex_ddg_args, stdout=outfile, stderr=subprocess.STDOUT, close_fds = True, cwd = output_directory)
    returncode = process.wait()
    outfile.close()

if __name__ == '__main__':
    cases = []
    for nstruct_i in range(1, nstruct + 1 ):
        for case_name in os.listdir('inputs'):
            case_path = os.path.join( 'inputs', case_name )
            for f in os.listdir(case_path):
                if f.endswith('.pdb'):
                    input_pdb_path = os.path.join( case_path, f )
                    break

            with open( os.path.join( case_path, 'chains_to_move.txt' ), 'r' ) as f:
                chains_to_move = f.readlines()[0].strip()

            cases.append( (case_name, case_path, input_pdb_path, chains_to_move, nstruct_i) )

    if use_multiprocessing:
        pool = multiprocessing.Pool( processes = min(max_cpus, multiprocessing.cpu_count()) )

    for args in cases:
        if use_multiprocessing:
            pool.apply_async( run_flex_ddg, args = args )
        else:
            run_flex_ddg( *args )

    if use_multiprocessing:
        pool.close()
        pool.join()

                                    '''
                f.write(content)
        sub_path = os.listdir(input_path_root)
        if config["flex_ddg"]["split_inputs"]["split"]:
            each_num = config["flex_ddg"]["split_inputs"]["each_have"]
        else:
            each_num = 9999999
        for each_sub in sub_path:
            split_dir_index = 0
            src_path = os.path.join(input_path_root,each_sub)
            file_list = os.listdir (src_path)
            len_file = len(file_list) 
            index = 0
            while 1 :
                dst_path = os.path.join(output_path,each_sub,str(split_dir_index))
                os.makedirs(dst_path,exist_ok=True)
                for i in range(each_num):
                    if index < len_file:
                        src_file = os.path.join(src_path,file_list[index])
                        dst_file = os.path.join(dst_path,"inputs",file_list[index])
                        shutil.copytree(src_file,dst_file)
                        cp_ddg(dst_path)
                        index += 1 
                    else:
                        break
                if index >= len_file:
                    break
                split_dir_index += 1
        
        logging.info("please run run_flex.py manually!")
            
        
        
        

    