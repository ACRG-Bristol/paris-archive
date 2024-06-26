# Template repo for using with agage-archive to create netcdf data archives

To use this repo:

- Use this repository as a template for a new respository by clicking on the "Use this template" button on the top-right side of the Github page
- Clone your new repository
- Make sure that the ```agage-archive``` repo and its dependencies are in your environment
- Rename the code folder ```project_archive``` to ```{project_name}_archive``` (it's important that the folder name ends in ```_archive```)
- Rename the ```project``` directory in the ```data``` folder with the name of your project
- Edit ```data/project/attributes.json```
- Edit the readme that will be included in your data archive in ```data/project/README.md```
- Edit the repository readme
- Run the ```config.py``` script and edit the template config file in that is generated in ```project_archive/config.yaml``` with the correct paths to your data and output directories/zip archives
- Edit the spreadsheets ```data_combination.xlsx```, ```data_exclude.xlsx```, ```data_release_schedule.xlsx``` and ```scale_defaults.csv``` in the ```data/project``` folder as appropriate for your dataset. Remove or add sites, instruments and species as required
- Edit ```project_archive/run.py``` with the appropriate code to process your data files. Import and use functions from agage-archive as much as possible, for standardisation
- Run the run script
- Edit the project name in ```notebooks/visualisation.ipynb``` to view the archive contents
