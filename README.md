### This repository stores Genome Apps and code to help manage Genome Apps<p>

#### How Genome Apps are run in Guardiome software
1. Genome App is clicked on the UI
2. client passes the name of the path to the Genome App to the server
3. the server calls run_genome_app in genomeapp.py of this repository
4. run_gneome_app function checks if Genome App has already been run by looking for the 'date_run' key in result.json in the Genome App results folder
5. If 'date_run' isn't in result.json, or result.json doesnt exist yet, run_genome_app calls the run_app function in the run.py file in the tools folder of the Genome App
6. run_app runs the Genome App causing result.json to be saved in the results folder of the Genome App
7. run_genome_app adds the 'date_run' key to result.json and saves result.json
8. run_genome_app then reads result.json and passes the result.json dictionary to the client
9. client displays Genome App result
