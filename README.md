PIT

- Requirements
  - python3
  - docker


install docker image with 
```
docker pull nirkoty/pit
```

Running PIT requires writting a configuration file that contains the paths to the input files, defines the conditions and samples in the
experiment and mass spectrometry runs.
Below is a template of how to write the config file. Examples are provided in the "examples" folder.
```
{
	"output": "PATH TO OUTPUT FOLDER",
	"reference_fasta": "PATH TO REFERENCE GENOME FASTA FILE FOR REFERENCE GUIDED ASSEMBLY",
    "reference_gff": "PATH TO REFERENCE GENOME ANNOTATION (GTF) FILE FOR REFERENCE GUIDED ASSEMBLY",
    "threads": 12,


    "ms":{
        "runs": {

            "RUN_NAME":{
                "files": [PATHS TO RAW FILES],
                "modifications": {
                        "fixed": [
                               USE MAXQUANT NOTATION FOR PTM
                        ],
                        "variable": [
                                
                        ]

                    },
                "TMT": {
                    "Nsi/1": "126",
                    "Nsi/2": "127C",
                    "Nsi/3": "128C",
                    "Nsi/4": "129C",
                    "si/1": "127N",
                    "si/2": "128N",
                    "si/3": "129N",
                    "si/4": "130N",

                }
            }
            
        }
        
    },

    "mutations": true,
    "splicing": true


	"conditions": {
            "CONDITION NAME": {


                
                    "SAMPLE NAME":{

                        "left": "PATH TO FASTQ FORWARD FILE",
                        "right": "PATH TO FASTQ REVERSE FILE"
                    }
                
            }
        }

    
}
```
The analysis can be started with docker using the command
```
python3 LaunchDocker.py -c PATH_TO_CONFIG_FILE -r
```
or with singularity using the command 
```
python3 LaunchDocker.py -c PATH_TO_CONFIG_FILE -r -s PATH_TO_SINGULARITY_IMAGE
```