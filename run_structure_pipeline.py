from rocketry import Rocketry
from rocketry.args import Return, Task
from rocketry.conds import every, after_success, after_fail, after_finish

import argparse

from rich import print
from rich.console import Console


from shared.pipeline import load_config, read_lockfile, write_lockfile


from steps.structure import increment_localpdb_version, update_localpdb
from steps.structure import localpdb_query
from steps.structure import fetch_sabdab_data, fetch_stcrdab_data, download_coordinates
from steps.structure import classify_peptide



filename = None
config = None
mhc_class = None

console = Console()

parser = argparse.ArgumentParser()
parser.add_argument("--mhc_class")



app = Rocketry()



### Conditions ###

@app.cond()
def is_first_task(this_task=Task()):
    return this_task.name.startswith("do_before")


@app.cond()
def is_complete():
    if read_lockfile(filename) == 'completed':
        return True 
    else:
        return False

steps = [
    localpdb_query,
    fetch_sabdab_data,
    fetch_stcrdab_data,
    download_coordinates,
    classify_peptide
]


def output_step_title(n:int, step_name:str):
    step_title = step_name.replace('_',' ').capitalize()
    print(f'Step {n + 1} of {len(steps)} | {step_title}')
    print ()


def run_action(process, mhc_class):
    action_number = process['next_action']
    step_action = steps[action_number]
    step_name = step_action.__name__
    
    output_step_title(action_number, step_name)
    success, errors, error_count = step_action(mhc_class)

    process['next_action'] += 1
    process['actions'][step_name] = {
        'success':success,
        'errors':errors,
        'error_count':error_count
    }

    return process


### Tasks ###

@app.task(every('1 seconds') & is_first_task)
def do_before():
    if read_lockfile(filename) != 'initialised':
        print (f'To run the pipeline again remove the lock file "{filename}"')
        exit()
    else:
        process = {'mhc_class':mhc_class, 'next_action':0, 'actions':{}}
        increment_localpdb_version(config)
        #update_localpdb(config)
        print ('Starting pipeline...')
        write_lockfile(filename, 'running')
        return process


@app.task(after_success(do_before))
def do_first(arg=Return(do_before)):
    return run_action(arg, mhc_class)


@app.task(after_success(do_first))
def do_second(arg=Return(do_first)):
    return run_action(arg, mhc_class)


@app.task(after_success(do_second))
def do_third(arg=Return(do_second)):
    return run_action(arg, mhc_class)


@app.task(after_success(do_third))
def do_fourth(arg=Return(do_third)):
    return run_action(arg, mhc_class)


@app.task(after_success(do_fourth))
def do_fifth(arg=Return(do_fourth)):
    return run_action(arg, mhc_class)


@app.task(after_success(do_fifth))
def do_sixth(arg=Return(do_fifth)):
    return run_action(arg, mhc_class)


@app.task(after_success(do_sixth))
def do_seventh(arg=Return(do_sixth)):
    return run_action(arg, mhc_class)

# assign last task to the function name of the final action
last_task = do_fifth

@app.task(after_success(last_task))
def do_last(arg=Return(last_task)):
    print ('Last task - cleanup')
    print (arg)
    print ('Pipeline completed')
    write_lockfile(filename, 'completed')
    
    

if __name__ == "__main__":
    print()
    console.rule('[bold cyan]histo.fyi structure pipeline[/bold cyan]')
    config = load_config()
    args = parser.parse_args()
    mhc_class = args.mhc_class
    if not mhc_class:
        print ('[red]error![/red]MHC class is a required argument, please try again')
    else:
        filename = f'tmp/structure_{mhc_class}.lock'
        app.run()
    print()

