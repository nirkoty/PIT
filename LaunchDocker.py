import argparse
import json
import subprocess
import shutil
import os
from pathlib import Path
import re

def validate(config):

    if "output" not in config:
        raise Exception("Please specify an output path")
    else:
        if not isinstance(config["output"], str):
            raise Exception("Please specify a correct output path")
        path = Path(config["output"])
        if not path.parent.absolute().is_dir():
            raise Exception("Please specify a correct output path")

    if "conditions" not in config or not isinstance(config["conditions"], dict):
        raise Exception("Please specify conditions")

    if len(config["conditions"])==0:
        raise Exception("Please specify samples")

    if "reference_fasta" in config:
        if not isinstance(config["reference_fasta"], str) or not Path(config["reference_fasta"]).is_file():
            raise Exception("Please specify a correct reference fasta file")
        if not "reference_gff" in config:
            raise Exception("Please specify a reference gff file")

    if "reference_gff" in config:
        if not isinstance(config["reference_gff"], str) or not Path(config["reference_gff"]).is_file():
            raise Exception("Please specify a correct reference gff file")
        if not "reference_fasta" in config:
            raise Exception("Please specify a reference fasta file")

    if "species" in config and (not isinstance(config["species"], str) or re.search(r"[A-Za-z]+\s[A-Za-z]+", config["species"]) is None):
        raise Exception("Species name must be in binomial notation")


    for condition, conditionObj in config["conditions"].items():
        if not isinstance(conditionObj, dict):
            raise Exception("Expecting a dictionary for condition "+condition)

        for sample, sampleObj in conditionObj.items():
            if len(sampleObj) == 0 or (("left" not in sampleObj and "right" not in sampleObj) and "single" not in sampleObj):
                raise Exception(f"Expecting RNA-Seq reads files for {condition} {sample}")

            if "left" in sampleObj:
                if not isinstance(sampleObj["left"], str) or not Path(sampleObj["left"]).is_file():
                    raise Exception(f"RNA-Seq reads file not found for {condition} {sample}: {sampleObj['left']}")

                if "right" not in sampleObj:
                    raise Exception(f"Expecting a right hand RNA-Seq reads file for {condition} {sample}")

            if "right" in sampleObj:
                if not isinstance(sampleObj["right"], str) or not Path(sampleObj["right"]).is_file():
                    raise Exception(f"RNA-Seq reads file not found for {condition} {sample}: {sampleObj['right']}")

                if "left" not in sampleObj:
                    raise Exception(f"Expecting a left hand RNA-Seq reads file for {condition} {sample}")

    if "mutations" in config and not isinstance(config["mutations"], bool):
        raise Exception(f"Expecting true or false value for mutations")

    if not isinstance(config["ms"]["runs"], dict):
        raise Exception("Expecting dictionnary with run names as keys")

    if "ms" in config and "runs" in config["ms"]:
        for runName, runObj in config["ms"]["runs"].items():
            if not "files" in runObj:
                raise Exception(f"Please specify files for run {runName}")

            if isinstance(runObj["files"], str):
                if not Path(runObj["files"]).is_file():
                        raise Exception(f"File not found for run {runName}: {runObj['files']}")
            else:
                for file in runObj["files"]:
                    if not isinstance(file, str) or not Path(file).is_file():
                        raise Exception(f"File not found for run {runName}: {file}")
            if "SILAC" in runObj:
                if not isinstance(runObj["SILAC"], dict):
                    raise Exception(f"Expecting dictionnary of SILAC labels for conditions in {runName} run")
                for silacCondition, labels in runObj["SILAC"].items():
                    if silacCondition not in config["conditions"]:
                        raise Exception(f"Condition {silacCondition} used in SILAC labels for run {runName} is not defined as a condition in the RNA-Seq")
                    if not isinstance(labels["label"], list) or len([x for x in labels["label"] if not isinstance(x, str)])>0:
                        raise Exception(f"Expecting a list of SILAC labels for run {runName}")

                    if "samples" in labels:
                        for sample in labels["samples"]:
                            if not sample in config["conditions"][silacCondition]:
                                raise Exception(f"Sample {sample} used in {runName} run is not defined inside {silacCondition} condition")

            elif "TMT" in runObj:
                tmtLabels = ["126", "127C", "128C", "129C", "127N", "128N", "129N", "130N", "130C", "131"]
                for label in runObj["TMT"].values():
                    if label not in tmtLabels:
                        raise Exception(f"{label} is not a valid TMT label in run {runName}")
                for sample in list(runObj["TMT"].keys()):
                    if not "/" in sample and sample in config["conditions"] and len(config["conditions"][sample])==1:
                        runObj["TMT"][f"{sample}/1"] = runObj["TMT"][sample]
                        del runObj["TMT"][sample]

            elif "ITRAQ" in runObj:
                pass
            else: #label free
                if not "condition" in runObj:
                    raise Exception(f"Please select a condition to use as database for run {runName}")
                if not runObj["condition"] in config["conditions"]:
                    raise Exception(f"Condition {runObj['condition']} used for run {runName} is not defined")



        if "combine" in config["ms"]:
            for combinedRun, combinedObj in config["ms"]["combine"].items():
                for run in combinedObj["runs"]:
                    if not run in config["ms"]["runs"]:
                        raise Exception(f"Run {run} used in {combinedRun} combined run is not defined")

def moveStart(config):

    filesMoved = []



    output = config["output"]

    os.makedirs(output, exist_ok=True)
    with open(f"{output}/config_base.json", "w") as f:
        json.dump(config, f, indent=4)

    print("Copying files")


    config["output"] = f"/project/"
    configPath = f"{output}/config_docker.json"


    try:

        for condition, conditionObj in config["conditions"].items():
            os.makedirs(f"{output}/{condition}", exist_ok=True)
            for sample, sampleObj in conditionObj.items():
                os.makedirs(f"{output}/{condition}/{sample}",exist_ok=True)
                if "left" in sampleObj:
                    if sampleObj["left"]!=f"{output}/{condition}/{sample}/{os.path.basename(sampleObj['left'])}":
                        print(sampleObj["left"])
                        shutil.copy(sampleObj["left"], f"{output}/{condition}/{sample}/{os.path.basename(sampleObj['left'])}")
                        filesMoved.append([sampleObj["left"], f"{output}/{condition}/{sample}/{os.path.basename(sampleObj['left'])}"])
                    sampleObj["left"] = f"/project/{condition}/{sample}/{os.path.basename(sampleObj['left'])}"
                if "right" in sampleObj:
                    if sampleObj["right"]!=f"{output}/{condition}/{sample}/{os.path.basename(sampleObj['right'])}":
                        print(sampleObj["right"])
                        shutil.copy(sampleObj["right"], f"{output}/{condition}/{sample}/{os.path.basename(sampleObj['right'])}")
                        filesMoved.append([sampleObj["right"], f"{output}/{condition}/{sample}/{os.path.basename(sampleObj['right'])}"])
                    sampleObj["right"] = f"/project/{condition}/{sample}/{os.path.basename(sampleObj['right'])}"
                if "single" in sampleObj:
                    if sampleObj["single"]!=f"{output}/{condition}/{sample}/{os.path.basename(sampleObj['single'])}":
                        print(sampleObj["single"])
                        shutil.copy(sampleObj["single"], f"{output}/{condition}/{sample}/{os.path.basename(sampleObj['single'])}")
                        filesMoved.append([sampleObj["single"], f"{output}/{condition}/{sample}/{os.path.basename(sampleObj['single'])}"])
                    sampleObj["single"] = f"/project/{condition}/{sample}/{os.path.basename(sampleObj['single'])}"

        if "reference_fasta" in config:
            if config["reference_fasta"]!=f"{output}/{os.path.basename(config['reference_fasta'])}":
                print("reference_fasta")
                print(config["reference_fasta"])
                shutil.copy(config["reference_fasta"], f"{output}/{os.path.basename(config['reference_fasta'])}")
            config["reference_fasta"] = f"/project/{os.path.basename(config['reference_fasta'])}"
        if "reference_gff" in config:
            if config["reference_gff"]!=f"{output}/{os.path.basename(config['reference_gff'])}":
                print(config["reference_gff"])
                shutil.copy(config["reference_gff"], f"{output}/{os.path.basename(config['reference_gff'])}")
            config["reference_gff"] = f"/project/{os.path.basename(config['reference_gff'])}"

        os.makedirs(f"{output}/ms", exist_ok=True)
        for run, runObj in config["ms"]["runs"].items():
            os.makedirs(f"{output}/ms/{run}", exist_ok=True)
            newFiles = []
            if isinstance(runObj["files"], str):
                newFiles.append(f"/project/ms/{run}/{os.path.basename(runObj['files'])}")
                if runObj['files']!=f"{output}/ms/{run}/{os.path.basename(runObj['files'])}":
                    print(runObj['files'])
                    shutil.copy(runObj['files'], f"{output}/ms/{run}/{os.path.basename(runObj['files'])}")
                    filesMoved.append([runObj['files'], f"{output}/ms/{run}/{os.path.basename(runObj['files'])}"])
                runObj["files"] = newFiles
            else:
                for file in runObj["files"]:
                    newFiles.append(f"/project/ms/{run}/{os.path.basename(file)}")
                    if file!=f"{output}/ms/{run}/{os.path.basename(file)}":
                        print(file)
                        shutil.copy(file, f"{output}/ms/{run}/{os.path.basename(file)}")
                        filesMoved.append([file, f"{output}/ms/{run}/{os.path.basename(file)}"])
                runObj["files"] = newFiles

        with open(configPath, "w") as f:
            json.dump(config, f, indent=4)
    except Exception as e:
        print(e)
        # for files in filesMoved:
        #     shutil.copy(files[1], files[0])

    return output


def moveEnd(configPath):
    config = json.load(open(configPath))
    configBase = json.load(open(f"{os.path.dirname(os.path.realpath(configPath))}/config_base.json"))

    config["output"] = f"/project/"
    for condition, conditionObj in config["conditions"].items():
        for sample, sampleObj in conditionObj.items():
            if "left" in sampleObj:
                shutil.move(sampleObj["left"], configBase[condition][sample]["left"])
            if "right" in sampleObj:
                shutil.move(sampleObj["right"], configBase[condition][sample]["right"])
            if "single" in sampleObj:
                shutil.move(sampleObj["single"], configBase[condition][sample]["single"])

    for run, runObj in config["mzml"]["runs"].items():
        i=0
        for file in runObj["files"]:
            shutil.move(file, f"{configBase['ms'][run]['files'][i]}")
            i+=1

    shutil.move(f"{os.path.dirname(os.path.realpath(configPath))}/config_base", f"{os.path.dirname(os.path.realpath(configPath))}/config.json")




if __name__ == '__main__':



    parser = argparse.ArgumentParser(description='This script starts the analysis')

    parser.add_argument("-c", "--config", nargs=1, required=True, help="path to config file", metavar="PATH")
    parser.add_argument("-s", "--singularity", nargs=1, required=False, help="path to config file", metavar="PATH")
    parser.add_argument("-r", "--run", action="store_true")


    args = parser.parse_args()

    config = json.load(open(args.config[0]))

    try:
        validate(config)
    except Exception as e:
        print(e)


    outputDir = moveStart(config)

    if args.run:

        if args.singularity:
            if args.singularity[0].endswith(".sif"):
                os.system(f"singularity build -F --sandbox {outputDir}/pit_sandbox {args.singularity[0]}")
                process = subprocess.Popen(
                    ['singularity', 'run', "-e", "-B", f'{outputDir}:/project/', f'{outputDir}/pit_sandbox'])
            else:
                process = subprocess.Popen(['singularity', 'run', '--writable', "-e", "-B", f'{outputDir}:/project/', f"{args.singularity[0]}"])


        else:
            process = subprocess.Popen(['docker', 'run', '-it', '-v', f'{outputDir}:/project/', 'pit'])
        try:
            print('Running in process', process.pid)
            process.wait()
        except subprocess.TimeoutExpired:
            print('Timed out - killing', process.pid)
            process.kill()
        print("Done")


    # except Exception as e:
    #     print(e)

    # finally:
    #     moveEnd(args.config[0])
