### 2019_1_22
* Pulled tools from perry tool file and directed to the correct directories and tools.
* Started running snpsift rule successfully, awaiting completion.
* Need to check where established part of pipeline ends to connect to this and what features if any are pulled.
* If no features processed or pulled, need to add to pipeline to pull those which could take an amount of additional work.

### 2019_1_14
* Pulled more scripts from previous project and updated constant file getting ready to test run

### 2019_1_11 2
* Seems as though the pipeline somewhat simultaneously predicts and creates hdf5 db while querying for the results
* Important for plugging this pipeline in is what type of vcf files are needed for matrix creation and from snpeff
* Normalized split vcf by sample and chr are needed for matrix creation while collapsed vcf files are needed from snpeff for vcf anno and snpsift
* Pulled over scripts and const, need to clean some stuff up with those

### 2019_1_11
* Beginning to compile the rules into one file as requested also adding comments to give some information to the rules
* In addition to the rule comments, adding sectional comments to explain next step in pipeline and the original file of the code
* Have gemini and query rules in file but unsure of next file/step after query, need to find out and lay out entire pipeline

### 2019_1_10 2
* So pipeline needs to be from running snpeff on to the end with decision tree prediction
* This will help with less of the pipeline to organize and pull over from one or two repos
* It seems as though steps are to organize snpeff output, build gemini, query gemini, make matrices, run prediciton
* Need to find if additional steps are between any of those steps and how much filtering/organization is involved

### 2019_1_10
* Initialize project repo
* Begin to pull code from repos nb_convergence and target_qc
* It is tricky that the project is in two repos and functions of the decision tree require all samples analyzed so it may be best to work backward
