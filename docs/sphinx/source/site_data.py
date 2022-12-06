import datetime
import json
from pathlib import Path
import utils_docs as utils
from utils_docs import mp # same instance of mpmath's mp imported in ivcurves/utils


def load_scores_database():
    with open('scores_database.json', 'r') as file:
        return json.load(file)


def submissions():
    """
    Get the submissions in the scores database, excluding those marked broken.

    Returns
    -------
        A dict mapping pull request numbers to submission data.
    """
    database = load_scores_database()
    return {pr_number: submission_data
            for pr_number, submission_data in database.items()
            if 'broken' not in submission_data.keys()}


def to_ghuser(username):
    return f':ghuser:`{username}`'


def to_pull(pr_number):
    return f':pull:`{pr_number}`'


def link_to_submission_code(commit, username, submission_main):
    return f'https://github.com/cwhanse/ivcurves/tree/{commit}/submissions/{username}/{submission_main}'


def link_to_submission_docs(username, submission_main):
    return f'submissions/{username}/{Path(submission_main).stem}.html'


def datetime_from_github_datetime_str(ghdatetime_str):
    return datetime.datetime.strptime(ghdatetime_str, '%Y-%m-%dT%H:%M:%SZ')


def leaderboard_entry_list():
    entries = []

    for pr_number, submission_data in submissions().items():
        code_link = link_to_submission_code(
            submission_data['merge_commit'],
            submission_data['username'],
            submission_data['submission_main']
        )
        docs_link = link_to_submission_docs(
            submission_data['username'],
            submission_data['submission_main']
        )
        entries.append({
            'submission': f'{to_ghuser(submission_data["username"])} ({to_pull(pr_number)})',
            'overall_score': sum(mp.mpmathify(v) for v in submission_data['test_sets'].values()),
            'submission_datetime': datetime_from_github_datetime_str(submission_data['submission_datetime']),
            'links':  f'`Code <{code_link}>`__, `Docs <{docs_link}>`__'
        })

    # order entries from lowest score to highest, and then by submission datetime
    entries.sort(key=lambda l: l['submission_datetime'])
    entries.sort(key=lambda l: l['overall_score'])

    for idx, entry in enumerate(entries):
        entry['rank'] = f'#{idx + 1}'
        entry['overall_score'] = mp.nstr(entry['overall_score'])
        # use the submission date instead of datetime
        entry['submission_date'] = entry['submission_datetime'].strftime('%m/%d/%Y')

    return entries


def compare_submissions_table_rows():
    entries = []

    for pr_number, submission_data in submissions().items():
        code_link = link_to_submission_code(
            submission_data['merge_commit'],
            submission_data['username'],
            submission_data['submission_main']
        )
        docs_link = link_to_submission_docs(
            submission_data['username'],
            submission_data['submission_main']
        )

        entry = {'Submission': f'{to_ghuser(submission_data["username"])} ({to_pull(pr_number)})'}

        # per case scores
        for name, score in submission_data['test_sets'].items():
            entry[name] = mp.nstr(mp.mpmathify(score))

        entry['Links'] = f'`Code <{code_link}>`__, `Docs <{docs_link}>`__'

        entries.append(entry)

    # convert entries to a table with the keys of an entry as its columns
    table_rows = []
    table_rows.append(list(entries[0].keys()))

    for entry in entries:
        table_rows.append(list(entry.values()))

    return table_rows


def test_set_name_to_parameters_and_image():
    mapping = {}
    test_set_names = list(utils.get_filenames_in_directory(utils.TEST_SETS_DIR))
    test_set_names.sort()
    for name in test_set_names:
        info = {}
        for idx, parameters in utils.read_iv_curve_parameter_sets(f'{utils.TEST_SETS_DIR}/{name}').items():
            info[idx] = {'title': f'Case {idx}',
                         'parameters': list(map(mp.nstr, parameters[:-1])), # exclude cells_in_series
                         'image_path': f'./_images/test_cases/{utils.make_iv_curve_name(name, idx)}.png'}
        mapping[name] = info
    return mapping

