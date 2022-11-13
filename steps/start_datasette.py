import os

# TODO something to create



db_remove_command = 'rm ../warehouse/datacompilations/core.db'

db_create_command = 'csvs-to-sqlite ../warehouse/datacompilations/core.csv ../warehouse/datacompilations/core.db'

start_command = 'datasette ../warehouse/datacompilations/core.db -o --setting max_returned_rows 10000'

os.system(start_command)