# nl_repo Overview
This repository contains pipelines and scripts for processing chemical and biological data, including SMILES-based compound analysis and activity/planning pipelines. The primary scripts are located in query_pipelines/{cmp,act,pln}_pipeline, with additional utilities in scripts and specialized directories (20250326_literature_spider, etc.).

This README provides instructions to clone the repository, set up the environment, and run the scripts in these directories scripts.
Prerequisites

SSH Access: You need SSH access to the EC2 instance with the provided key file (<key.pem>) and public IP (<ec2-ip>).
GitHub Access: Ensure you can access https://github.com/sensoriumtx/nl_repo.git (no authentication required for public repo).
Local Setup (optional for VS Code):
Install Visual Studio Code with the Remote - SSH and Dev Containers extensions.



# Setup Instructions
1. SSH into EC2
On your local machine, connect to the EC2 instance:
bash```
ssh -i <key.pem> ubuntu@<ec2-ip>
```

2. Enter Docker Container
Access the running Docker container:
bash```
docker exec -it "container_name" bash
```
    *utilize docker ps for containernames

3. Clone the Repository
In the container, clone the repository:
bash ```
cd /home/rstudio/knowledge-graph
git clone https://github.com/sensoriumtx/nl_repo.git nl_repo
```

If the nl_repo directory already exists, update it:
bash```
cd /home/rstudio/knowledge-graph/nl_repo
git pull origin main
```

4. Fix Git Ownership
Git may report a “dubious ownership” error due to file ownership in the container. Fix it:

bash```
git config --global --add safe.directory /home/rstudio/knowledge-graph/nl_repo
```

Verify:

bash```
cd /home/rstudio/knowledge-graph/nl_repo
git status
```

Expected output: “On branch main, Your branch is up to date with ‘origin/main’.”


# Running Scripts
The repository contains R, Python, and shell scripts in query_pipelines, scripts, and other subdirectories. Key scripts are in query_pipelines/{cmp,act,pln}_pipeline.

Directory Structure

query_pipelines/cmp_pipeline:
smiles_pipeline.r: Processes SMILES data for compounds.
cmp_pipeline.r: Main compound pipeline.
wildcard_cmp.r: Wildcard compound processing.
    archive/: deprecated versions.


query_pipelines/act_pipeline:
act_pipeline.r: Activity pipeline.
act_for_cmp_and_pln_pipeline.r: Supports compound/planning pipelines.
archive/: deprecated versions.


query_pipelines/pln_pipeline:
pln_pipeline.r: Planning pipeline.
batch_pln_pipeline.r: Batch processing for planning.
archive/: deprecated versions.


scripts/: Utility scripts (R base, Python).
Other: These scripts are direct graph pulls or post-hoc procedures for data enrichment.

