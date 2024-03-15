# SeaLion_container

## Create a Container 

To improve the usability of the workflow, an apptainer/singularity container was designed and build. 
This should allow a broader and user friendliy distribution of the software. Also this will improve the reproducability of analysis performed by the workflow.

## Build Container

To build the container locally inside the code directory we can run the following commands:

For Apptainer: 
```
cd ./code
apptainer build SeaLion_container.sif SeaLion_container.def
```

For Singularity:
```
cd ./code
singularity build SeaLion_container.sif SeaLion_container.def
```


The above command needs to be excecuted in the same directory as the perl scirpts, which will be imported in the container during the build process.


## Run software 

The workflow can be started with the following command from within the directory with the test data.


```
cd testdata
CONTAINER="../code/SeaLion_container.sif" # provide the full path to the location of your container.

apptainer exec ${CONTAINER} sealion_gold_beta4.pl -s
# when using singularity
singularity exec ${CONTAINER} sealion_gold_beta4.pl -s
```


For further information about Apptainer and Singularity Containers, please have a look at the official documentation:

https://apptainer.org/docs/user/main/quick_start.html

https://docs.sylabs.io/guides/3.5/user-guide/quick_start.html

