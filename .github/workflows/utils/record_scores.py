import argparse
import csv
import enum
import json
import pathlib


ROOT_DIR = pathlib.Path(f'{__file__}/../../../..').resolve()
TEST_SETS_DIR = pathlib.Path(f'{ROOT_DIR}/test_sets')


def load_json(filename):
    with open(filename, 'r') as file:
        return json.load(file)


def save_json(json_dict, filename):
    with open(filename, 'w') as file:
        return json.dump(json_dict, file, indent=2)


def load_overall_scores(filename):
    """
    Buildings a dictionary from test set filenames (excluding file extensions)
    to strings represented scores by reading from a CSV file with these
    columns: test_set, and score.

    Parameters
    ----------
    filename : str
        The path to the CSV file. The file extension must be included.

    Returns
    -------
    dict
    """
    overall_scores = {}
    with open(filename, newline='') as file:
        reader = csv.DictReader(file)
        for row in reader:
            overall_scores[row['test_set']] = row['score']
    return overall_scores


def test_set_filenames():
    """
    Returns a set of filenames in the directory ``directory_path``.
    The filenames do not have file extensions.

    Returns
    -------
    set
        A set of filenames without file extensions.
    """
    return {entry.stem for entry in TEST_SETS_DIR.iterdir() if entry.is_file()}


def validate_overall_scores(overall_scores):
    """
    Given a dictionary of overall scores (see :func:`load_overall_scores`),
    the following are validated:

    - There is a score for every test set.
    - Each score is parsable as a float.

    Otherwise, a ``ValueError`` is raised.

    Parameters
    ----------
    overall_scores : dict
        A dictionary from test set filenames (excluding file extensions) to
        strings representing scores.
    """
    overall_scores_keys = set(overall_scores.keys())

    # Proceed with normal validation
    valid_test_set_filenames = test_set_filenames()
    missing_test_set_filenames = valid_test_set_filenames - overall_scores_keys

    if missing_test_set_filenames:
        raise ValueError(f'Missing scores from these test sets: {", ".join(missing_test_set_filenames)}')

    for name, score_str in overall_scores.items():
        try:
            float(score_str) # validate is a number
        except ValueError:
            raise ValueError("The score of test set '{name}' must parse to a float: {score_str}")


def write_overall_scores_to_database(database, pr_number, pr_author, pr_closed_datetime, overall_scores):
    """
    Writes an entry in the JSON scores database.
    Entries are of this form:

    .. code-block:: json

       "<pr_number>": {
           "username": "<pr_author>",
           "submission_datetime": "<pr_closed_datetime>",
           "test_sets": {
               "<test_set_filename>": "<score>"
               ...
           }
       }

    As tables, the database has two tables: submissions, and test_set_scores.
    The submissions table has these columns: ``pr_number``, ``username``,
    ``submission_datetime``, and ``test_sets``.
    The primary key is ``pr_number``.
    The test_set_scores table has these columns: ``test_set``, and ``score``.
    The primary key is a compound key of ``pr_number`` and ``test_set``.

    Parameters
    ----------
    database : dict
        The scores database.

    pr_number : int
        The pull request number.

    pr_author : str
        The author of the pull request.

    overall_scores : dict
        A dictionary from test set filenames (excluding file extensions) to
        strings representing scores.
    """
    database[pr_number] = {'username': pr_author,
                           'submission_datetime': pr_closed_datetime,
                           'test_sets': overall_scores}


def get_argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--pr-author', dest='pr_author', type=str,
                        help='GitHub username of the pull request author')
    parser.add_argument('--pr-number', dest='pr_number', type=int,
                        help='GitHub pull request number')
    parser.add_argument('--pr-closed-at', dest='pr_closed_datetime', type=str,
                        help='Datetime when the GitHub pull request closed')
    parser.add_argument('--overall-scores-path', dest='overall_scores_path',
                        type=str, help='Path to the CSV of overall scores')
    parser.add_argument('--database-path', dest='database_path', type=str,
                        help='Path to the JSON scores database')
    return parser


if __name__ == '__main__':
    args = get_argparser().parse_args()

    overall_scores = load_overall_scores(args.overall_scores_path)
    scorer_code = validate_overall_scores(overall_scores)
    database = load_json(args.database_path)
    write_overall_scores_to_database(database, args.pr_number, args.pr_author, args.pr_closed_datetime, overall_scores)
    save_json(database, f'{ROOT_DIR}/{args.database_path}')

