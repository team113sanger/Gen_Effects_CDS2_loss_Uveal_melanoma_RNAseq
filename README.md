# r_project_template

## Getting started
This repo serves as minimal template for setting up R analysis projects in team113. 
Feel free to edit as you'd like. 

If you're interested in reproducing your R environment, the project contains an renv file, which records the 
dependencies 

```
renv::init() # To initialise revn
renv::install("glue") # To add a package
renv::snapshot() # To record any packages you've used
renv::restore() # To rebuild an environment from the renv.lockfile
````


```
git remote add origin https://gitlab.internal.sanger.ac.uk/team113_projects/training/templates/r_project_template.git
git branch -M main
git push -uf origin main
```