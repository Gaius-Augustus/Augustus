import git
import json


def commit_info(path_to_repo):
    repo = git.Repo.init(path_to_repo)
    return repo.commit().committed_datetime.isoformat(), str(repo.commit())


def store_additional_data(rev, cdate, exec_minutes, filename):
    data = {
        'revision': rev,
        'softmasking': True,
        'commit_date': cdate,
        'execution_time': exec_minutes
    }

    with open(filename, 'w') as file:
        json.dump(data, file, indent=4)
