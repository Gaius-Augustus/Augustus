import git
import json
import os
import psutil

def commit_info(path_to_repo):
    # path_to_repo is the working tree or its .git folder. Open the working tree: git.Repo.init
    # on the .git folder makes a new empty repository in it, which has no commit.
    path = os.path.abspath(path_to_repo)
    if os.path.basename(path) == '.git':
        path = os.path.dirname(path)
    repo = git.Repo(path)
    return repo.commit().committed_datetime.isoformat(), str(repo.commit())


def store_additional_data(rev, cdate, resources, filename):
    data = {
        'revision': rev,
        'softmasking': True,
        'commit_date': cdate,
        'resources': resources
    }

    with open(filename, 'w') as file:
        json.dump(data, file, indent=4)


def check_memory():
    # startup memory check
    m = psutil.virtual_memory()
    print(f'Total memory: {m.total*10**(-9):.2f} GB')
    print(f'Currently used memory: {m.used*10**(-9):.2f} GB')
