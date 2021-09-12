# readme

All scripts are tested in python3.6.3, The whole process can be divided into two major steps: relax and flex ddG.

You can run ***batch.py*** directly:

```bash
python3 batch.py
```

 Or you can run ***relax.py*** to relax and then run ***make_flex_ddg.py*** to make flex_ddg input files:

If you ran relax first.py, then you can pass the path to relax_sc.list as an optional parameter to make_flex_ddg.py

```bash
python3 relax.py
python3 make_flex_ddg.py [PATH_TO relax_sc.list]
```

After running, you should get a folder named flex_ddg (You can modify the name of the folder through "root_path" in config.yaml), then you need to run run_flex.py in the folder:

```bash
python3 run_flex.py
```

Then you will get the running result in the subfolder "output", you need to run analyze_flex_ddG.py to analyze the running result:

```bash
python3 analyze_flex_ddG.py flex_ddg/output
```

If you want to get the wild type or mutant type structure, you can run the following command:

```bash
python3 extract_structures.py output
```
Our project is based on flex_ddg, their github address is https://github.com/Kortemme-Lab/flex_ddG_tutorial.git
