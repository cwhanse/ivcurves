import datetime
import json
import utils_docs as utils
from mpmath import mp


def load_scores_database():
    with open('scores_database.json', 'r') as file:
        return json.load(file)


def to_ghuser(username):
    return f':ghuser:`{username}`'


def to_pull(pr_number):
    return f':pull:`{pr_number}`'


def date_from_github_datetime_str(ghdatetime_str):
    ghdatetime = datetime.datetime.strptime(ghdatetime_str, '%Y-%m-%dT%H:%M:%SZ')
    return ghdatetime.strftime('%m/%d/%Y')


def leaderboard_entry_list():
    database = load_scores_database()
    leaderboard_entries = []

    for pr_number, submission_data in database.items():
        leaderboard_entries.append({
            'pr_number': to_pull(pr_number),
            'username': to_ghuser(submission_data['username']),
            'overall_score': sum(mp.mpmathify(v) for v in submission_data['test_sets'].values()),
            'submission_date': date_from_github_datetime_str(submission_data['submission_datetime'])
        })

    # order entries from lowest score to highest
    leaderboard_entries.sort(key=lambda l: l['overall_score'])

    for idx, entry in enumerate(leaderboard_entries):
        entry['rank'] = f'#{idx + 1}'
        entry['overall_score'] = mp.nstr(entry['overall_score'])

    return leaderboard_entries


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

